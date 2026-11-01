"""A DataProvider over a RefSeq GFF3 annotation and an indexed genomic FASTA.

:class:`RefSeqDataProvider` builds transcript models from the GFF3's exon and CDS
features and reads every sequence it serves (genomic, transcript and translated
protein) from the FASTA, through :class:`SequenceProxy`. It needs ``pysam``,
installed by the ``refseq`` extra (``pip install hgvs-weaver[refseq]``).
"""

from __future__ import annotations

import collections
import gzip
import json
import logging
import os
import typing

import pysam.libcfaidx

import weaver

logger = logging.getLogger(__name__)


def _check_range(ac: str, start: int, end: int | None) -> None:
    """Refuses a range that starts before the sequence or ends before it starts.

    Raises:
        weaver.DataProviderError: Naming ``ac`` and the range.
    """
    if start < 0 or (end is not None and end < start):
        raise weaver.DataProviderError(f"{ac}: [{start}, {end}) is not a range")


class IdentifierKind:
    """Enum for sequence identifier types."""

    Genomic = "g"
    Transcript = "c"
    Protein = "p"
    Gene = "gene"


class SequenceProxy:
    """Proxy for accessing genomic sequences, supporting recording and replay."""

    def __init__(self, fasta_path: str, cache_path: str | None = None, mode: str | None = None) -> None:
        self.fasta_path = fasta_path
        self.cache_path = cache_path or os.environ.get("WEAVER_SEQ_CACHE")
        self.mode = mode or os.environ.get("WEAVER_SEQ_MODE", "live")
        self.cache: dict[str, str] = {}
        # Replay only: every recorded base by accession and 0-based position, so
        # that any requested range can be served, not just the ranges recorded.
        self._known: dict[str, dict[int, str]] = {}
        self.fasta: typing.Any = None
        self.references: list[str] = []

        if self.mode == "replay":
            if not self.cache_path:
                raise ValueError("replay mode needs a sequence cache: pass cache_path or set WEAVER_SEQ_CACHE")
            with open(self.cache_path) as f:
                self.cache = json.load(f)

            self._known = self._index_recorded_bases(self.cache)

            # Try to get references from manifest first
            manifest_data = self.cache.get("_manifest")
            if isinstance(manifest_data, dict):
                manifest = typing.cast("dict[str, list[str]]", manifest_data)
                fname = os.path.basename(self.fasta_path)
                if fname in manifest:
                    self.references = manifest[fname]
            else:
                # Fallback: extract unique accessions from keys like "AC:start-end"
                refs = set()
                for key in self.cache:
                    if ":" in key and not key.startswith("_"):
                        refs.add(key.split(":")[0])
                self.references = list(refs)
        else:
            self.fasta = pysam.libcfaidx.FastaFile(fasta_path)
            self.references = list(self.fasta.references)

    def fetch(self, ac: str, start: int, end: int | None = None) -> str:
        """The bases of ``ac`` in ``[start, end)``, to the end of the sequence where ``end`` is None.

        A range past the end of the sequence answers the bases there, none if it starts past the end.

        Raises:
            weaver.DataProviderError: If ``ac`` is not held, or ``[start, end)`` is not a range.
        """
        _check_range(ac, start, end)
        key = f"{ac}:{start}-{end}"

        if self.mode == "replay":
            if key in self.cache:
                return self.cache[key]
            return self._replay_range(ac, start, end)

        if ac not in self.fasta:
            raise weaver.DataProviderError(f"{ac}: not in the reference data")
        seq = str(self.fasta.fetch(ac, start, end).upper())
        if self.mode == "record":
            self.cache[key] = seq
        return seq

    @staticmethod
    def _index_recorded_bases(cache: dict[str, str]) -> dict[str, dict[int, str]]:
        """Spreads recorded "AC:start-end" fragments into per-position bases."""
        known: dict[str, dict[int, str]] = {}
        for key, seq in cache.items():
            if key.startswith("_") or ":" not in key or not isinstance(seq, str):
                continue
            ac, _, rng = key.rpartition(":")
            start_s, _, _ = rng.partition("-")
            try:
                start = int(start_s)
            except ValueError:
                continue
            bases = known.setdefault(ac, {})
            for i, base in enumerate(seq):
                bases[start + i] = base
        return known

    def _replay_range(self, ac: str, start: int, end: int | None) -> str:
        """Serves a range from recorded bases, with N for positions never recorded.

        A caller may page through a sequence in blocks far wider than any recorded
        fragment; returning the full block keeps it from concluding the sequence
        ends where the recording does. A base that was never recorded reads as N,
        so a comparison that depends on it fails rather than silently passing.
        """
        bases = self._known.get(ac)
        if not bases:
            raise weaver.DataProviderError(f"{ac}: not in the recorded sequences {self.cache_path}")
        if end is None:
            end = max(bases) + 1
        return "".join(bases.get(i, "N") for i in range(start, max(start, end)))

    def save_cache(self) -> None:
        if self.mode == "record" and self.cache_path:
            # Merge with existing cache if it exists
            existing = {}
            if os.path.exists(self.cache_path):
                with open(self.cache_path) as f:
                    existing = json.load(f)

            # Update sequences
            existing.update(self.cache)

            # Update manifest
            manifest = existing.setdefault("_manifest", {})
            fname = os.path.basename(self.fasta_path)
            manifest[fname] = self.references

            os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
            with open(self.cache_path, "w") as f:
                json.dump(existing, f, indent=2)
            logger.info("Saved cache with %d sequences and manifest for %s", len(existing) - 1, fname)


class RefSeqDataProvider:
    """DataProvider implementation using RefSeq GFF and FASTA files."""

    def __init__(self, gff_path: str, fasta_path: str) -> None:
        """Initializes the provider and loads the GFF file.

        Args:
          gff_path: Path to the RefSeq GFF3 file (can be gzipped).
          fasta_path: Path to the indexed RefSeq genomic FASTA file.
        """
        self.gff_path = gff_path
        self.fasta_path = fasta_path
        # (tx_ac, chrom) -> TranscriptData
        self.transcripts: dict[tuple[str, str], typing.Any] = {}
        # tx_ac -> list of ref_ac
        self.tx_to_refs: dict[str, list[str]] = collections.defaultdict(list)
        self.gene_to_transcripts: dict[str, list[str]] = collections.defaultdict(list)
        self.chrom_to_transcripts: dict[str, list[typing.Any]] = collections.defaultdict(list)
        self.accession_map: dict[str, tuple[str, str]] = {}  # protein_id -> tx_id
        self.mane_select: dict[str, str] = {}  # gene_name -> tx_ac (MANE Select)

        self._load_gff()
        self.fasta = SequenceProxy(fasta_path)

        self._transcript_cache: dict[tuple[str, str], str] = {}

    def _load_gff(self) -> None:
        """Parses the GFF file into internal transcript models."""
        logger.info("Loading RefSeq GFF from %s into memory...", self.gff_path)

        # (tx_ac, chrom) -> { 'exons': set(), 'cds': set(), 'info': {} }
        tx_data: dict[tuple[str, str], dict[str, typing.Any]] = collections.defaultdict(
            lambda: {"exons": set(), "cds": set(), "info": {}},
        )
        genes: dict[str, str] = {}
        id_to_tx_ac: dict[str, str] = {}

        opener = gzip.open if self.gff_path.endswith(".gz") else open

        with opener(self.gff_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                chrom, source, feature_type, start, end, _score, strand, _frame, attr_str = parts
                attrs: dict[str, str] = {}
                for item in attr_str.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs[k] = v

                feat_id = attrs.get("ID")
                parent = attrs.get("Parent")
                tx_ac = attrs.get("transcript_id")

                if not tx_ac:
                    if feat_id and (feat_id.startswith(("rna-NM_", "rna-NR_"))):
                        tx_ac = feat_id[4:]
                    elif parent and (parent.startswith(("rna-NM_", "rna-NR_"))):
                        tx_ac = parent[4:]
                    elif parent in id_to_tx_ac:
                        tx_ac = id_to_tx_ac[parent]

                if feature_type == "gene":
                    if feat_id:
                        genes[feat_id] = attrs.get("gene", attrs.get("Name", ""))

                elif feature_type in [
                    "mRNA",
                    "transcript",
                    "tRNA",
                    "ncRNA",
                    "lnc_RNA",
                    "rRNA",
                    "scRNA",
                    "snRNA",
                    "snoRNA",
                ]:
                    if tx_ac:
                        if feat_id:
                            id_to_tx_ac[feat_id] = tx_ac
                        key = (tx_ac, chrom)
                        if not tx_data[key]["info"] or source in ["RefSeq", "BestRefSeq"]:
                            tx_data[key]["info"] = {
                                "strand": strand,
                                "parent": parent,
                                "gene_name": attrs.get("gene", ""),
                            }
                        tag = attrs.get("tag", "")
                        if "MANE Select" in tag or "MANE_Select" in tag:
                            gene_name = attrs.get("gene", "")
                            if gene_name:
                                self.mane_select[gene_name] = tx_ac

                elif feature_type == "exon":
                    if tx_ac:
                        tx_data[(tx_ac, chrom)]["exons"].add((int(start), int(end)))

                elif feature_type == "CDS" and tx_ac:
                    tx_data[(tx_ac, chrom)]["cds"].add((int(start), int(end)))
                    prot_id = attrs.get("protein_id")
                    if prot_id:
                        tx_data[(tx_ac, chrom)]["info"]["protein_id"] = prot_id

        logger.info("Finalizing %d transcript-reference pairs...", len(tx_data))
        for (tx_id, chrom), data in tx_data.items():
            info = data["info"]
            if not info or "strand" not in info:
                continue

            exons_raw: list[tuple[int, int]] = list(data["exons"])
            if not exons_raw:
                continue

            cds: list[tuple[int, int]] = list(data["cds"])
            strand = info["strand"]
            protein_id = info.get("protein_id", "")

            py_exons: list[dict[str, typing.Any]] = []
            current_tx_pos = 0
            exons_for_tx = sorted(exons_raw, key=lambda x: x[0], reverse=(strand == "-"))

            g_min = min(e[0] for e in exons_raw)
            g_max = max(e[1] for e in exons_raw)

            for i, (e_start, e_end) in enumerate(exons_for_tx):
                exon_len = e_end - e_start + 1
                py_exons.append(
                    {
                        "ord": i,
                        "transcript_start": current_tx_pos,
                        "transcript_end": current_tx_pos + exon_len,
                        "reference_start": e_start - 1,
                        "reference_end": e_end - 1,
                        "alt_strand": 1 if strand == "+" else -1,
                        "cigar": f"{exon_len}=",
                    },
                )
                current_tx_pos += exon_len

            cds_start_idx = None
            cds_end_idx = None
            if cds:
                # Biological start codon base
                g_cds_start = min(c[0] for c in cds) if strand == "+" else max(c[1] for c in cds)
                # Biological last base of stop codon
                g_cds_end = max(c[1] for c in cds) if strand == "+" else min(c[0] for c in cds)

                cds_start_idx = self._genomic_to_tx(g_cds_start, py_exons, "+" if strand == "+" else "-")
                cds_end_idx = self._genomic_to_tx(g_cds_end, py_exons, "+" if strand == "+" else "-")

            gene_name = info["gene_name"]
            if not gene_name and info["parent"]:
                gene_name = genes.get(info["parent"], "")

            record = {
                "ac": tx_id,
                "gene": gene_name,
                "cds_start_index": cds_start_idx,
                "cds_end_index": cds_end_idx,
                "strand": 1 if strand == "+" else -1,
                "reference_accession": chrom,
                "exons": py_exons,
                "protein_id": protein_id,
                "start": g_min,
                "end": g_max,
            }
            self.transcripts[(tx_id, chrom)] = record
            self.tx_to_refs[tx_id].append(chrom)
            self.chrom_to_transcripts[chrom].append(record)

            if protein_id:
                self.accession_map[protein_id] = (tx_id, chrom)
            if gene_name:
                self.gene_to_transcripts[gene_name].append(tx_id)

        logger.info("GFF loading complete.")

    def _genomic_to_tx(self, g_pos: int, exons: list[dict[str, typing.Any]], strand: str) -> int | None:
        """Maps genomic position to transcript index.

        Args:
          g_pos: 1-based genomic coordinate.
          exons: List of exon dictionaries.
          strand: "+" or "-.

        Returns:
          0-based transcript index or None if not in exons.
        """
        g_0 = g_pos - 1
        for exon in exons:
            if exon["reference_start"] <= g_0 <= exon["reference_end"]:
                if strand == "+":
                    return exon["transcript_start"] + (g_0 - exon["reference_start"])
                return exon["transcript_start"] + (exon["reference_end"] - g_0)
        return None

    def get_seq(
        self,
        ac: str,
        start: int,
        end: int | None,
        kind: str | weaver.IdentifierType,
        force_plus: bool = False,
    ) -> str:
        """Retrieves the bases or residues of ``ac`` in the 0-based half-open range ``[start, end)``.

        A protein is the translation of its transcript's CDS; a transcript is spliced from the genome.

        Raises:
            weaver.DataProviderError: If the provider holds no sequence ``ac`` of that kind, or ``[start, end)`` is
                not a range. A range past the end of a sequence it holds is not an error: it answers the bases
                there, none if it starts past the end.
        """
        _check_range(ac, start, end)
        kind_name = str(kind).lower()
        if "proteinaccession" in kind_name or kind_name in {"p", "protein_accession"}:
            res = self.accession_map.get(ac)
            if not res:
                raise weaver.DataProviderError(f"{ac}: not in the reference data")
            tx_ac, chrom = res
            tx_info = self.transcripts[(tx_ac, chrom)]
            cds_start = tx_info["cds_start_index"]
            cds_end = tx_info["cds_end_index"]
            if cds_start is None or cds_end is None:
                raise weaver.DataProviderError(f"{ac}: the exons of {tx_ac} on {chrom} do not hold its CDS")

            tx_seq = self._get_tx_seq(tx_ac, chrom, 0, None, force_plus=force_plus).upper()
            translated = self._translate_cds(tx_seq[cds_start : cds_end + 1])
            return translated[start:end]

        if "transcript" in kind_name or kind_name == "c":
            return self._get_tx_seq(ac, self._default_reference(ac), start, end, force_plus=force_plus).upper()

        return self.fasta.fetch(ac, start, end)

    def _get_full_tx_seq(self, tx_ac: str, ref_ac: str) -> str:
        """Builds and caches the full transcript sequence (respecting strand)."""
        tx = self.transcripts[(tx_ac, ref_ac)]
        seq_parts = []
        # Exons are already ordered 5' -> 3' in tx["exons"]
        for exon in tx["exons"]:
            s = self.fasta.fetch(tx["reference_accession"], exon["reference_start"], exon["reference_end"] + 1)
            if tx["strand"] == -1:
                s = self.reverse_complement(s)
            seq_parts.append(s)
        return "".join(seq_parts)

    def _get_full_tx_seq_cached(self, tx_ac: str, ref_ac: str) -> str:
        try:
            return self._transcript_cache[(tx_ac, ref_ac)]
        except KeyError:
            self._transcript_cache[(tx_ac, ref_ac)] = self._get_full_tx_seq(tx_ac, ref_ac)
            return self._transcript_cache[(tx_ac, ref_ac)]

    def _get_tx_seq(self, tx_ac: str, ref_ac: str, start: int, end: int | None, force_plus: bool = False) -> str:
        """Builds a transcript sequence from genomic exons with optional transformations."""
        tx = self.transcripts[(tx_ac, ref_ac)]
        full_seq = self._get_full_tx_seq_cached(tx_ac, ref_ac)

        if force_plus and tx["strand"] == -1:
            # reverse_complement(natural_seq) gives genomic-plus orientation
            full_seq = self.reverse_complement(full_seq)

        return full_seq[start:end]

    def get_transcript(self, transcript_ac: str, reference_ac: str | None) -> dict[str, typing.Any]:
        """Returns the transcript model for the given accession."""
        if reference_ac:
            tx = self.transcripts.get((transcript_ac, reference_ac))
            if tx is None:
                raise weaver.DataProviderError(f"{transcript_ac}: not placed on {reference_ac}")
            return tx

        return self.transcripts[(transcript_ac, self._default_reference(transcript_ac))]

    def _default_reference(self, tx_ac: str) -> str:
        """The reference a transcript is read on when none is named: a primary chromosome in the FASTA if any.

        Raises:
            weaver.DataProviderError: If the GFF3 places no transcript under exactly ``tx_ac``, version included.
        """
        refs = self.tx_to_refs.get(tx_ac)
        if not refs:
            raise weaver.DataProviderError(f"{tx_ac}: not in the reference data")
        ref_ac = next((r for r in refs if r.startswith("NC_0000") and r in self.fasta.references), None)
        if not ref_ac:
            ref_ac = next((r for r in refs if r in self.fasta.references), refs[0])
        return ref_ac

    def get_transcripts_for_region(self, chrom: str, start: int, end: int) -> list[str]:
        """Finds transcripts overlapping a genomic region."""
        transcripts = self.chrom_to_transcripts.get(chrom, [])
        results = []
        for tx in transcripts:
            # Check overlap: (start1 <= end2) and (end1 >= start2)
            if start <= tx["end"] and end >= tx["start"]:
                results.append(tx["ac"])
        return list(set(results))

    def get_symbol_accessions(
        self,
        symbol: str,
        source_kind: str,
        target_kind: str,
    ) -> list[tuple[weaver.IdentifierType, str]]:
        """The accessions ``symbol`` maps to, by the kinds weaver asks for; any other pair is empty.

        ``gene -> c``: a gene's transcripts; ``gene -> p``: a gene's proteins; ``c -> p``: a transcript's protein;
        ``p -> c``: a protein's transcript.
        """
        if source_kind == IdentifierKind.Gene and target_kind == IdentifierKind.Transcript:
            return [
                (weaver.IdentifierType.TranscriptAccession, tx_ac)
                for tx_ac in dict.fromkeys(self.gene_to_transcripts.get(symbol, []))
            ]
        if source_kind == IdentifierKind.Gene and target_kind == IdentifierKind.Protein:
            proteins = (self._protein_of(tx_ac) for tx_ac in self.gene_to_transcripts.get(symbol, []))
            return [(weaver.IdentifierType.ProteinAccession, ac) for ac in dict.fromkeys(proteins) if ac]
        if source_kind == IdentifierKind.Transcript and target_kind == IdentifierKind.Protein:
            protein = self._protein_of(symbol) if symbol in self.tx_to_refs else ""
            return [(weaver.IdentifierType.ProteinAccession, protein)] if protein else []
        if source_kind == IdentifierKind.Protein and target_kind == IdentifierKind.Transcript:
            if symbol not in self.accession_map:
                return []
            tx_ac, _chrom = self.accession_map[symbol]
            return [(weaver.IdentifierType.TranscriptAccession, tx_ac)]
        return []

    def _protein_of(self, tx_ac: str) -> str:
        """The protein the GFF3 names for a transcript it annotates; empty for a non-coding one."""
        return self.transcripts[(tx_ac, self._default_reference(tx_ac))]["protein_id"]

    def get_identifier_type(self, identifier: str) -> str:
        """What the provider holds under ``identifier``, whatever its form; ``unknown`` if it holds nothing.

        A transcript or protein the GFF3 annotates, a sequence in the FASTA, or a gene the GFF3 names.
        """
        if identifier in self.tx_to_refs:
            return "transcript_accession"
        if identifier in self.accession_map:
            return "protein_accession"
        if identifier in self.fasta.references:
            return "genomic_accession"
        if identifier in self.gene_to_transcripts:
            return "gene_symbol"
        return "unknown"

    def reverse_complement(self, seq: str) -> str:
        """Returns the reverse complement of a sequence."""
        complement = str.maketrans("ATCGNautcgn", "TAGCNtagcgn")
        return seq.translate(complement)[::-1]

    def _translate_cds(self, seq: str) -> str:
        """Translates a CDS nucleotide sequence into an amino acid sequence (1-letter codes)."""
        codon_table = {
            "ATA": "I",
            "ATC": "I",
            "ATT": "I",
            "ATG": "M",
            "ACA": "T",
            "ACC": "T",
            "ACG": "T",
            "ACT": "T",
            "AAC": "N",
            "AAT": "N",
            "AAA": "K",
            "AAG": "K",
            "AGC": "S",
            "AGT": "S",
            "AGA": "R",
            "AGG": "R",
            "CTA": "L",
            "CTC": "L",
            "CTG": "L",
            "CTT": "L",
            "CCA": "P",
            "CCC": "P",
            "CCG": "P",
            "CCT": "P",
            "CAC": "H",
            "CAT": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGA": "R",
            "CGC": "R",
            "CGG": "R",
            "CGT": "R",
            "GTA": "V",
            "GTC": "V",
            "GTG": "V",
            "GTT": "V",
            "GCA": "A",
            "GCC": "A",
            "GCG": "A",
            "GCT": "A",
            "GAC": "D",
            "GAT": "D",
            "GAA": "E",
            "GAG": "E",
            "GGA": "G",
            "GGC": "G",
            "GGG": "G",
            "GGT": "G",
            "TCA": "S",
            "TCC": "S",
            "TCG": "S",
            "TCT": "S",
            "TTC": "F",
            "TTT": "F",
            "TTA": "L",
            "TTG": "L",
            "TAC": "Y",
            "TAT": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGC": "C",
            "TGT": "C",
            "TGA": "*",
            "TGG": "W",
        }
        protein = []
        seq = seq.upper()
        for i in range(0, len(seq) - len(seq) % 3, 3):
            codon = seq[i : i + 3]
            protein.append(codon_table.get(codon, "X"))
        return "".join(protein)

    def to_json(self) -> None:
        return None
