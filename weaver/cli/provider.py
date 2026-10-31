"""The biocommons ``hgvs`` adapter the validation compares weaver against."""

# The biocommons Interface this module subclasses declares its methods without annotations.
# pyright: reportIncompatibleMethodOverride=false

from __future__ import annotations

import typing

import hgvs.dataproviders.interface

import weaver
from weaver import refseq

# Importable from here for callers written against 1.0, where the class lived in this module.
RefSeqDataProvider = refseq.RefSeqDataProvider


class ReferenceHgvsDataProvider(hgvs.dataproviders.interface.Interface):
    """Bridge between weaver DataProvider and hgvs library Interface.

    The biocommons ``Interface`` declares its abstract methods without annotations, so every typed override reads
    as incompatible to a checker; the overrides follow the library's documented signatures, and the file-level
    pyright comment at the top turns that one report off for this module.
    """

    def __init__(self, refseq_provider: weaver.DataProvider) -> None:
        self.url = "local://refseq"
        self.required_version = "1.1"
        super().__init__()
        self.rp = refseq_provider

    def get_seq(self, ac: str, start: int | None = None, end: int | None = None) -> str:
        kind = "g" if ac.startswith("NC_") else "c"
        if ac.startswith("NP_"):
            kind = "p"
        return self.rp.get_seq(ac, start or 0, end, kind)

    def get_tx_info(
        self,
        tx_ac: str,
        alt_ac: str | None = None,
        _alt_aln_method: str | None = None,
    ) -> dict[str, typing.Any] | None:
        try:
            tx = self.rp.get_transcript(tx_ac, alt_ac)
            return {
                "hgnc": tx["gene"],
                "cds_start_i": tx["cds_start_index"],
                "cds_end_i": tx["cds_end_index"] + 1 if tx["cds_end_index"] is not None else None,
                "strand": tx["strand"],
                "alt_ac": tx["reference_accession"],
                "alt_aln_method": "transcript",
            }
        except Exception:
            return None

    def get_tx_exons(
        self,
        tx_ac: str,
        alt_ac: str | None = None,
        _alt_aln_method: str | None = None,
    ) -> list[dict[str, typing.Any]] | None:
        try:
            tx = self.rp.get_transcript(tx_ac, alt_ac)
            res: list[dict[str, typing.Any]] = []
            exons_transcript = sorted(tx["exons"], key=lambda x: x["transcript_start"])
            for i, e in enumerate(exons_transcript):
                res.append(
                    {
                        "tx_ac": tx_ac,
                        "alt_ac": tx["reference_accession"],
                        "alt_aln_method": "transcript",
                        "ord": i,
                        "tx_start_i": e["transcript_start"],
                        "tx_end_i": e["transcript_end"],
                        "alt_start_i": e["reference_start"],
                        "alt_end_i": e["reference_end"] + 1,
                        "alt_strand": tx["strand"],
                        "cigar": e["cigar"],
                    },
                )
            res.sort(key=lambda x: x["alt_start_i"])
            return res
        except Exception:
            return None

    def get_tx_identity_info(self, tx_ac: str) -> dict[str, typing.Any] | None:
        try:
            tx = self.rp.get_transcript(tx_ac, None)
            total_len = tx["exons"][-1]["transcript_end"]
            return {
                "hgnc": tx["gene"],
                "lengths": [total_len],
                "tx_ac": tx_ac,
                "alt_acs": [tx["reference_accession"]],
                "cds_start_i": tx["cds_start_index"],
                "cds_end_i": tx["cds_end_index"] + 1 if tx["cds_end_index"] is not None else None,
            }
        except Exception:
            return None

    def data_version(self) -> str:
        return "1.1"

    def schema_version(self) -> str:
        return "1.1"

    def get_assembly_map(self, _assembly_name: str) -> dict[str, typing.Any]:
        return {}

    def get_gene_info(self, _gene: str) -> dict[str, typing.Any]:
        return {}

    def get_pro_ac_for_tx_ac(self, tx_ac: str) -> str | None:
        try:
            tx = self.rp.get_transcript(tx_ac, None)
            return str(tx.get("protein_id"))
        except Exception:
            return None

    def get_tx_for_gene(self, _gene: str) -> list[str]:
        return []

    def get_tx_for_region(self, _alt_ac: str, _alt_aln_method: str, _start: int, _end: int) -> list[str]:
        return []

    def get_tx_mapping_options(self, tx_ac: str) -> list[dict[str, typing.Any]]:
        """Every reference the transcript is placed on where the provider lists them, else its default placement."""
        try:
            listed = getattr(self.rp, "tx_to_refs", None)
            refs = (
                listed.get(tx_ac, [])
                if listed is not None
                else [self.rp.get_transcript(tx_ac, None)["reference_accession"]]
            )
            return [{"tx_ac": tx_ac, "alt_ac": r, "alt_aln_method": "transcript"} for r in refs]
        except Exception:
            return []

    def get_similar_transcripts(self, _tx_ac: str) -> list[str]:
        return []

    def get_acs_for_protein_seq(self, _seq: str) -> list[str]:
        return []
