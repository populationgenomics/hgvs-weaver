"""Alleles in cis, c.[145C>T;147C>G], through the bindings: parsing, dicts, VRS."""

import base64
import hashlib
import json
import typing

import pytest

import weaver

# The transcript is genomic indices 10..50; c.1 is index 5 of it, so g. = c. + 15.
UTR5 = "GGGGG"
CDS = "ATGAAACTGGCCTATCGCTAA"  # M K L A Y R *
UTR3 = "CCGTATAAGTAAGG"
TRANSCRIPT = UTR5 + CDS + UTR3
GENOME = "T" * 10 + TRANSCRIPT + "C" * 10
PROTEIN = "MKLAYR"
SEQUENCES = {"NC_X.1": GENOME, "NM_X.1": TRANSCRIPT, "NP_X.1": PROTEIN}


def sha512t24u(seq: str) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(seq.encode()).digest()[:24]).decode()


class Provider:
    """One single-exon coding transcript, its genome and protein, with refget lookups."""

    def get_transcript(self, ac: str, _ref: str | None) -> dict[str, typing.Any]:
        if ac != "NM_X.1":
            raise weaver.DataProviderError(f"no transcript {ac}")
        return {
            "ac": ac,
            "gene": "TEST",
            "cds_start_index": len(UTR5),
            "cds_end_index": len(UTR5) + len(CDS) - 1,
            "strand": 1,
            "reference_accession": "NC_X.1",
            "exons": [
                {
                    "transcript_start": 0,
                    "transcript_end": len(TRANSCRIPT),
                    "reference_start": 10,
                    "reference_end": 10 + len(TRANSCRIPT) - 1,
                    "alt_strand": 1,
                    "cigar": f"{len(TRANSCRIPT)}=",
                },
            ],
        }

    def get_seq(self, ac: str, start: int, end: int | None, _kind: str) -> str:
        return SEQUENCES[ac][start:end]

    def get_symbol_accessions(self, symbol: str, _s: str, t: str) -> list[tuple[weaver.IdentifierType, str]]:
        if t == "p" and symbol == "NM_X.1":
            return [(weaver.IdentifierType.ProteinAccession, "NP_X.1")]
        return []

    def get_identifier_type(self, identifier: str) -> weaver.IdentifierType:
        if identifier.startswith("NC_"):
            return weaver.IdentifierType.GenomicAccession
        if identifier.startswith("NP_"):
            return weaver.IdentifierType.ProteinAccession
        return weaver.IdentifierType.TranscriptAccession

    def get_transcripts_for_region(self, _chrom: str, _start: int, _end: int) -> list[str]:
        return ["NM_X.1"]

    def get_refget_accession(self, _ac: str) -> str | None:
        return None  # let weaver compute it

    def get_accession_for_refget(self, refget: str) -> str | None:
        return next((ac for ac, seq in SEQUENCES.items() if f"SQ.{sha512t24u(seq)}" == refget), None)


@pytest.fixture
def mapper() -> weaver.VariantMapper:
    provider = Provider()
    return weaver.VariantMapper(provider, refget=provider)


@pytest.mark.parametrize(
    "hgvs",
    [
        "NC_X.1:g.[22C>T;28T>G]",
        "NC_012920.1:m.[8993T>G;9000del]",
        "NM_X.1:c.[7C>T;13T>G]",
        "NM_X.1(TEST):c.[122-6T>A;153C>T;200_201insA]",
        "NR_X.1:n.[7C>T;13T>G]",
        "NM_X.1:r.[7c>u;13u>g]",
        "NP_X.1:p.[Lys2Leu;Ala4del]",
        "NM_X.1:c.[7C>T]",
    ],
)
def test_cis_alleles_parse_and_format(hgvs: str) -> None:
    var = weaver.parse(hgvs)
    assert str(var) == var.format() == hgvs
    assert var.ac == hgvs.split("(")[0].split(":")[0]
    assert var.coordinate_type == hgvs.split(":")[1][0]
    members = var.members
    assert members is not None
    assert len(members) == hgvs.count(";") + 1
    for m in members:
        assert m.members is None
        assert m.ac == var.ac
        assert m.gene == var.gene
        assert m.coordinate_type == var.coordinate_type
    first_posedit = hgvs.split("[")[1].split(";")[0].rstrip("]")
    assert members[0].format() == f"{hgvs.split(':')[0]}:{var.coordinate_type}.{first_posedit}"


def test_plain_variants_have_no_members() -> None:
    assert weaver.parse("NM_X.1:c.7C>T").members is None


def test_alleles_in_trans_are_refused() -> None:
    with pytest.raises(weaver.UnsupportedOperationError, match="two molecules"):
        weaver.parse("NM_X.1:c.[7C>T];[13T>G]")
    with pytest.raises(weaver.ParseError):
        weaver.parse("NM_X.1:c.[7C>T;]")


def test_dict_round_trip() -> None:
    var = weaver.parse("NM_X.1(TEST):c.[7C>T;13T>G]")
    d = var.to_dict()
    assert d["variant_type"] == "CisPhased"
    assert d["ac"] == "NM_X.1"
    assert d["gene"] == "TEST"
    assert [m["variant_type"] for m in d["members"]] == ["Coding", "Coding"]
    assert str(weaver.Variant.from_dict(d)) == "NM_X.1(TEST):c.[7C>T;13T>G]"
    assert str(weaver.Variant.from_dict(json.loads(var.to_json()))) == str(var)


def test_mapping_methods_want_the_members(mapper: weaver.VariantMapper) -> None:
    cis = weaver.parse("NM_X.1:c.[7C>T;13T>G]")
    with pytest.raises(ValueError, match=r"allele in cis.*members"):
        mapper.c_to_g(cis)
    members = cis.members
    assert members is not None
    assert [str(mapper.c_to_g(m)) for m in members] == ["NC_X.1:g.22C>T", "NC_X.1:g.28T>G"]


def test_normalise_and_validate_go_member_by_member(mapper: weaver.VariantMapper) -> None:
    # c.4 is the first A of AAA (c.4_6): the deletion shifts to c.6.
    assert str(mapper.normalize_variant(weaver.parse("NM_X.1:c.[4del;13T>G]"))) == "NM_X.1:c.[6del;13T>G]"
    assert weaver.parse("NM_X.1:c.[7C>T;13T>G]").validate(Provider())
    assert not weaver.parse("NM_X.1:c.[7C>T;13A>G]").validate(Provider())


def test_to_vrs_is_a_cis_phased_block_of_the_members_alleles(mapper: weaver.VariantMapper) -> None:
    var = weaver.parse("NM_X.1:c.[7C>T;13T>G]")
    block = mapper.to_vrs(var)
    assert block["type"] == "CisPhasedBlock"
    assert block["id"].startswith("ga4gh:CPB.")
    assert block["id"] == f"ga4gh:CPB.{block['digest']}"
    assert mapper.vrs_id(var) == block["id"]
    members = var.members
    assert members is not None
    assert [m["id"] for m in block["members"]] == [mapper.vrs_id(m) for m in members]
    assert block["members"][0]["location"]["start"] == 21
    assert block["members"][0]["location"]["end"] == 22
    assert block["sequenceReference"] == block["members"][0]["location"]["sequenceReference"]
    assert block["sequenceReference"]["refgetAccession"] == f"SQ.{sha512t24u(GENOME)}"
    assert block["expressions"] == [{"syntax": "hgvs.c", "value": "NM_X.1:c.[7C>T;13T>G]"}]
    # The identifier does not depend on the order written; the members do.
    reversed_ = mapper.to_vrs(weaver.parse("NM_X.1:c.[13T>G;7C>T]"))
    assert reversed_["id"] == block["id"]
    assert reversed_["members"] == block["members"][::-1]
    assert mapper.vrs_id(weaver.parse("NC_X.1:g.[22C>T;28T>G]")) == block["id"]


@pytest.mark.parametrize(
    ("hgvs", "expected"),
    [
        ("NM_X.1:c.[7C>T;13T>G]", "NC_X.1:g.[22C>T;28T>G]"),
        ("NM_X.1:c.[4del;13T>G]", "NC_X.1:g.[21del;28T>G]"),
        ("NC_X.1:m.[22C>T;28T>G]", "NC_X.1:g.[22C>T;28T>G]"),
        ("NP_X.1:p.[Lys2Leu;Ala4del]", "NP_X.1:p.[Lys2Leu;Ala4del]"),
    ],
)
def test_from_vrs_round_trip(mapper: weaver.VariantMapper, hgvs: str, expected: str) -> None:
    block = mapper.to_vrs(weaver.parse(hgvs))
    back = mapper.from_vrs(block)
    assert str(back) == expected
    assert mapper.vrs_id(back) == block["id"]
    accession = expected.split(":")[0]
    assert str(mapper.from_vrs(json.dumps(block), accession)) == expected


def test_a_block_stating_its_sequence_once_reads_back(mapper: weaver.VariantMapper) -> None:
    def allele(start: int, alt: str) -> dict[str, typing.Any]:
        return {
            "type": "Allele",
            "location": {"type": "SequenceLocation", "start": start, "end": start + 1},
            "state": {"type": "LiteralSequenceExpression", "sequence": alt},
        }

    reference = {"type": "SequenceReference", "refgetAccession": f"SQ.{sha512t24u(GENOME)}"}
    block: dict[str, typing.Any] = {
        "type": "CisPhasedBlock",
        "members": [allele(21, "T"), allele(27, "G")],
        "sequenceReference": reference,
    }
    assert str(mapper.from_vrs(block)) == "NC_X.1:g.[22C>T;28T>G]"
    reference["refgetAccession"] = "SQ.x"
    with pytest.raises(weaver.HGVSError):
        mapper.from_vrs(block)


def test_equivalence_compares_the_sets_of_members(mapper: weaver.VariantMapper) -> None:
    searcher = Provider()
    a = weaver.parse("NM_X.1:c.[7C>T;13T>G]")

    def level(other: str) -> weaver.EquivalenceLevel:
        return mapper.equivalent_level(a, weaver.parse(other), searcher)

    assert level("NM_X.1:c.[13T>G;7C>T]") == weaver.EquivalenceLevel.Analogous
    assert level("NC_X.1:g.[22C>T;28T>G]") == weaver.EquivalenceLevel.Analogous
    assert level("NM_X.1:c.[7C>T;13T>A]") == weaver.EquivalenceLevel.Different
    assert level("NM_X.1:c.7C>T") == weaver.EquivalenceLevel.Different
    assert mapper.equivalent(weaver.parse("NM_X.1:c.[7C>T]"), weaver.parse("NC_X.1:g.22C>T"), searcher)
