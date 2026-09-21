"""DigestTable: refget accessions computed once from local sequences."""

import pathlib
import typing

import pytest

import weaver
from weaver.refget import DigestTable, sequence_digest

GENOME = "ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA"


def test_digest_follows_the_refget_normalisation_rule() -> None:
    # The specification's own example.
    assert sequence_digest("ACGT") == "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
    # Soft-masking and line breaks do not change the accession.
    assert sequence_digest("acgt") == sequence_digest("ACGT")
    assert sequence_digest("AC\nGT\n") == sequence_digest("ACGT")


def test_table_round_trips_through_tsv_and_answers_both_ways(tmp_path: pathlib.Path) -> None:
    digest = sequence_digest(GENOME)
    table = DigestTable([("NC_TEST.1", len(GENOME), digest), ("chr1", len(GENOME), digest)])
    path = tmp_path / "refget.tsv"
    table.to_tsv(path)
    loaded = DigestTable.from_tsv(path)
    assert len(loaded) == 2
    assert loaded.get_refget_accession("NC_TEST.1") == digest
    assert loaded.get_refget_accession("NC_OTHER.1") is None
    assert loaded.get_accession_for_refget(digest) == "NC_TEST.1"
    assert loaded.get_accession_for_refget(f"ga4gh:{digest}") == "NC_TEST.1"
    assert loaded.accessions_for_refget(digest) == ["NC_TEST.1", "chr1"]
    assert loaded.length("chr1") == len(GENOME)
    assert path.read_text().splitlines()[0] == "sequence\tlength\tsha512t24u"


def test_a_table_without_the_columns_is_refused(tmp_path: pathlib.Path) -> None:
    path = tmp_path / "bad.tsv"
    path.write_text("name\tdigest\nx\tSQ.y\n")
    with pytest.raises(ValueError, match="lacks columns"):
        DigestTable.from_tsv(path)
    with pytest.raises(ValueError, match="not a refget accession"):
        DigestTable([("x", 1, "md5:abc")])


class Provider:
    """Serves the soft-masked genome, as a FASTA-backed provider would after uppercasing."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        raise weaver.DataProviderError(f"no transcript {ac}")

    def get_seq(self, _ac: str, start: int, end: int | None, _kind: str) -> str:
        return GENOME[start:end]

    def get_symbol_accessions(self, _symbol: str, _s: str, _t: str) -> list[tuple[weaver.IdentifierType, str]]:
        return []

    def get_identifier_type(self, _identifier: str) -> weaver.IdentifierType:
        return weaver.IdentifierType.GenomicAccession


def test_the_table_agrees_with_what_the_mapper_computes_and_names_the_sequence_back() -> None:
    provider = Provider()
    computed = weaver.VariantMapper(provider).to_vrs(weaver.parse("NC_TEST.1:g.7G>C"))
    table = DigestTable([("NC_TEST.1", len(GENOME), sequence_digest(GENOME.lower()))])
    mapper = weaver.VariantMapper(provider, refget=table)
    allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.7G>C"))
    assert (
        allele["location"]["sequenceReference"]["refgetAccession"]
        == computed["location"]["sequenceReference"]["refgetAccession"]
    )
    assert str(mapper.from_vrs(allele)) == "NC_TEST.1:g.7G>C"


def test_from_fasta_hashes_the_uppercased_sequence(tmp_path: pathlib.Path) -> None:
    pysam = pytest.importorskip("pysam")
    fasta = tmp_path / "g.fa"
    fasta.write_text(f">NC_TEST.1 a test\n{GENOME[:20].lower()}\n{GENOME[20:]}\n")
    pysam.faidx(str(fasta))
    table = DigestTable.from_fasta(fasta)
    assert table.get_refget_accession("NC_TEST.1") == sequence_digest(GENOME)
    assert table.length("NC_TEST.1") == len(GENOME)
