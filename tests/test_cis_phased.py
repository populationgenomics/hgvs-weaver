"""Alleles in cis, c.[145C>T;147C>G], through the bindings: parsing, dicts, VRS."""

from __future__ import annotations

import base64
import hashlib
import json
import re
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


def test_an_allele_projects_to_the_genome_as_an_allele(mapper: weaver.VariantMapper) -> None:
    g = mapper.c_to_g(weaver.parse("NM_X.1:c.[7C>T;13T>G]"))
    members = g.members
    assert members is not None
    assert [str(m) for m in members] == ["NC_X.1:g.22C>T", "NC_X.1:g.28T>G"]


def test_an_allele_projects_from_the_genome(mapper: weaver.VariantMapper) -> None:
    c = mapper.g_to_c(weaver.parse("NC_X.1:g.[22C>T;28T>G]"), "NM_X.1")
    assert str(c) == "NM_X.1:c.[7C>T;13T>G]"


def test_every_transcript_gets_the_allele(mapper: weaver.VariantMapper) -> None:
    found = mapper.g_to_c_all(weaver.parse("NC_X.1:g.[22C>T;28T>G]"), Provider())
    assert [str(v) for v in found] == ["NM_X.1:c.[7C>T;13T>G]"]


def test_an_allele_is_respelled_in_rna(mapper: weaver.VariantMapper) -> None:
    r = mapper.c_to_r(weaver.parse("NM_X.1:c.[7C>T;13T>G]"))
    assert str(r) == "NM_X.1:r.[7c>u;13u>g]"


def test_an_rna_allele_is_respelled_in_dna(mapper: weaver.VariantMapper) -> None:
    c = mapper.r_to_c(weaver.parse("NM_X.1:r.[7c>u;13u>g]"))
    assert str(c) == "NM_X.1:c.[7C>T;13T>G]"


def test_an_rna_allele_projects_to_the_genome(mapper: weaver.VariantMapper) -> None:
    g = mapper.r_to_g(weaver.parse("NM_X.1:r.[7c>u;13u>g]"))
    assert str(g) == "NC_X.1:g.[22C>T;28T>G]"


def test_a_member_that_fails_raises_its_own_error_naming_it(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.ValidationError, match=r"Member 2 of NM_X\.1:c\.\[7C>T;\*20del\], NM_X\.1:c\.\*20del"):
        mapper.c_to_g(weaver.parse("NM_X.1:c.[7C>T;*20del]"))


def test_an_allele_of_another_system_is_refused(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(ValueError, match="allele in cis"):
        mapper.c_to_g(weaver.parse("NC_X.1:g.[22C>T;28T>G]"))


def test_protein_vrs_refuses_an_allele(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(ValueError, match="allele in cis"):
        mapper.protein_vrs(weaver.parse("NM_X.1:c.[7C>T;13T>G]"))


def test_c_to_p_predicts_the_protein_the_members_make_together(mapper: weaver.VariantMapper) -> None:
    predicted = mapper.c_to_p(weaver.parse("NM_X.1:c.[4A>T;6A>T]"))
    assert str(predicted) == "NP_X.1:p.(Lys2Tyr)"
    assert predicted.members is None


def test_c_to_p_writes_residues_stated_beside_a_change_as_a_protein_allele(mapper: weaver.VariantMapper) -> None:
    predicted = mapper.c_to_p(weaver.parse("NM_X.1:c.[13T>G;10_12=]"))
    assert str(predicted) == "NP_X.1:p.[(Ala4=;Tyr5Asp)]"
    assert predicted.members is not None
    assert [str(m) for m in predicted.members] == ["NP_X.1:p.(Ala4=)", "NP_X.1:p.(Tyr5Asp)"]


def test_c_to_p_of_an_allele_stating_no_whole_codon_raises_no_protein_statement(
    mapper: weaver.VariantMapper,
) -> None:
    with pytest.raises(weaver.NoProteinStatementError, match="no whole codon"):
        mapper.c_to_p(weaver.parse("NM_X.1:c.[10=;20=]"))


def test_c_to_p_refuses_an_allele_of_genomic_members(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(ValueError, match=r"allele in cis of c\. members"):
        mapper.c_to_p(weaver.parse("NC_X.1:g.[22C>T;28T>G]"))


def test_an_allele_normalised_to_one_change_comes_back_as_a_plain_variant(mapper: weaver.VariantMapper) -> None:
    normalised = mapper.normalize_variant(weaver.parse("NM_X.1:c.[7C>T;8T>G]"))
    assert normalised.members is None
    assert str(normalised) == "NM_X.1:c.7_8delinsTG"


def test_validate_reads_the_allele_whole() -> None:
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
    # Neither the identifier nor the members, in the normalised order, depend
    # on the order written.
    reversed_ = mapper.to_vrs(weaver.parse("NM_X.1:c.[13T>G;7C>T]"))
    assert reversed_["id"] == block["id"]
    assert reversed_["members"] == block["members"]
    assert mapper.vrs_id(weaver.parse("NC_X.1:g.[22C>T;28T>G]")) == block["id"]


def test_a_coding_allele_and_its_genomic_projection_share_a_vrs_identifier(mapper: weaver.VariantMapper) -> None:
    # c.4 and c.6 share a codon, which groups them in the c. text but not in VRS
    assert mapper.vrs_id(weaver.parse("NM_X.1:c.[4A>T;6A>T]")) == mapper.vrs_id(weaver.parse("NC_X.1:g.[19A>T;21A>T]"))


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


def test_judging_an_ambiguous_allele_raises_its_readings(mapper: weaver.VariantMapper) -> None:
    with pytest.raises(weaver.AmbiguousDescriptionError) as raised:
        mapper.equivalent_level(
            weaver.parse("NM_X.1:c.[7_8insA;7_8insG]"),
            weaver.parse("NM_X.1:c.7_8insAG"),
            Provider(),
        )
    assert sorted(str(r) for r in raised.value.readings) == ["NM_X.1:c.7_8insAG", "NM_X.1:c.7_8insGA"]


def test_an_ambiguity_is_an_ambiguous_description_error_with_its_readings() -> None:
    with pytest.raises(weaver.AmbiguousDescriptionError) as raised:
        weaver.parse("NM_X.1:c.[7_8insA;7_8insG]").validate(Provider())
    assert not isinstance(raised.value, weaver.ValidationError)
    assert sorted(str(r) for r in raised.value.readings) == ["NM_X.1:c.7_8insAG", "NM_X.1:c.7_8insGA"]


def test_validate_raises_on_a_contradictory_allele() -> None:
    with pytest.raises(weaver.ValidationError, match=re.escape("NM_X.1:c.7C>G")):
        weaver.parse("NM_X.1:c.[7C>T;7C>G]").validate(Provider())
