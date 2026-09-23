# The Choices Weaver Makes

HGVS leaves room for judgement, and the tools that read it have settled on different answers. This
page lists every deliberate choice weaver makes where another tool could reasonably answer
otherwise: what weaver does, why, the recommendation it rests on, and what other tools do.

Each choice ends with two lines. **Spec** links the [HGVS nomenclature](https://hgvs-nomenclature.org/stable/)
page (or the VRS or refget specification) the choice rests on, or says there is none. **Tests** names
the tests that assert it, as `file::function`; files are under `hgvs-weaver/tests/` unless written
as `src/…` (unit tests) or `tests/…py` (Python), and property tests are in
`hgvs-weaver/tests/properties/main.rs`. "Agrees" and "differs" are stated only where the other
tool's behaviour was checked, by running it (VariantValidator's REST API, the biocommons `hgvs`
grammar and regression tables, ClinVar's 100,000-variant sample) or by reading its code; anything
else is marked not checked.

The tools referred to: **HGVS** is the nomenclature; **biocommons** is the Python `hgvs` package
weaver was originally a port of; **VariantValidator** is the web service; **ClinVar** is the
descriptions in its variant summary; **VRS** is the GA4GH Variation Representation Specification
2.0.1.

## Projection between sequences

**A projection states the target's bases.** A transcript record and its genome can differ at a
base (RefSeq transcripts are curated against submitted mRNAs). When a variant is projected, the
edit is resolved on the source to the bases it removes and the bases it puts there; where the
target holds different bases over the projected range, weaver writes the minimal edit that turns
the target's bases into that alternate. MUC2 `NM_002457.5:c.12468C>A`, over a record C and a
genomic G, is `g.1099802G>A`; the record's `c.1568=` on SHANK3, over a genomic T, is `g.…T>C`;
and the genome's `g.…T>C` projected onto that record is `c.1568=`. A duplication or inversion
carries the source's bases, because that is what the molecule becomes: `c.12468dup` is
`g.1099802delinsCC`, not a duplication of the genome's G. A deletion is positional. Where record
and genome agree, which is nearly everywhere, the edit is carried across unchanged.
*Biocommons agrees* (`replace_reference`). *VariantValidator agrees* on substitution, identity,
duplication and deletion at the MUC2 base.
Spec: [general](https://hgvs-nomenclature.org/stable/recommendations/general/), "descriptions on
RNA/protein level should describe the changes observed on that level"; a description is of the
sequence it is on.
Tests: `projection_reference_test::a_projection_states_the_genomes_bases_not_the_records`,
`::a_projection_states_the_transcripts_bases_not_the_genomes`,
`::the_re_read_is_in_the_targets_orientation_on_the_minus_strand`;
`projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
(MUC2 and SHANK3 with NCBI's alignments); `a_variant_agrees_with_its_projection` (property, where
record and genome agree).

**A change across a splice junction has no genomic form.** `r.44_47del` spanning an exon boundary
describes the spliced RNA. weaver refuses to project it to the genome and says why, while its
`c.` spelling still projects (to a deletion that includes the intron) and its protein is still
predicted. *VariantValidator differs*: it projects the `c.` reading with the warning "spans at
least one intron". weaver refuses because a warning is easy to miss and the two molecules are
not the same.
Spec: [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/)
and [RNA splicing](https://hgvs-nomenclature.org/stable/recommendations/RNA/splicing/), where exon
skipping is written as an `r.` deletion across the junction (`r.(3277_3432del)`), an RNA event.
Tests: `rna_test::r_projects_to_the_genome_within_one_exon_only`, `::r_predicts_the_protein_like_c`;
`rna_across_a_junction_has_a_protein_but_no_genomic_form` (property).

**No transcript-side normalisation before projecting.** weaver projects the variant where it is
written. *VariantValidator differs*: it 3'-shifts on the transcript first (`c.1568dup` becomes
`c.1569dup`) and reports the automapping. Normalisation is available separately
(`normalize_variant`) and is the caller's decision.
Spec: none; the [3' rule](https://hgvs-nomenclature.org/stable/recommendations/general/) says how
to write a variant, not that a tool must rewrite what it is given.
Tests: `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
(positions come back as given).

**Gapped exons are projected locally.** Where a transcript exon aligns to the genome with
insertions and deletions, weaver projects through the exon's CIGAR and writes a local edit.
*VariantValidator differs* on heavily gapped exons: for SHANK3's `c.1568dup` it falls back to a
2.5 kb whole-region delins from `g.50695049` to `g.50697558`. weaver's `g.50697558delinsCC` is
the deliberate answer.
Spec: none.
Tests: `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
(the SHANK3 dup case); `src/transcript_mapper.rs::test_g_to_n_cigar`.

**A repeat is projected as its whole run.** `c.10AC[3]` on the minus strand is widened to the full
run of the unit on the transcript before projecting, so the other strand reads it from the right
end. Without this the genomic repeat started mid-run. *Not checked against other tools.*
Spec: [repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/):
the count is the total number of units in the run.
Tests: `transcript_coordinates_test::repeat_on_the_minus_strand_projects_to_its_whole_run`.

**A position outside every exon is an error, not an extrapolation.** Transcript positions resolve
through the exon structure only; a `c.` position that falls in no exon and has no intronic offset
is a typed error. *Biocommons* extrapolates in some paths. weaver does not, because a silently
extrapolated coordinate is wrong without saying so.
Spec: [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "it is not allowed
to describe variants in nucleotides beyond the boundaries of a reference sequence".
Tests: none directly; `transcript_coordinates_round_trip` (property) covers every in-exon position.

**Intronic positions are carried as given.** An intronic edit projected onto a transcript has no
transcript base to compare against, so the edit is kept as written; projected onto the genome,
the bases it states are checked against the genome like any other. Normalisation leaves intronic
edits where they are.
Spec: [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): RNA and coding
reference sequences "do not contain intron sequences and can therefore not be used to describe
variants affecting these sequences".
Tests: `projection_reference_test::an_intronic_position_has_no_transcript_base_to_re_read`;
`transcript_coordinates_test::normalize_leaves_intronic_coding_variants_as_written`;
`intron_offsets_agree_from_either_exon` (property).

**`m.` is `g.` on the mitochondrial reference; `r.` is `c.` or `n.` in RNA letters.** `m.` and
`g.` differ only in the letter they are written with; normalisation, SPDI, validation and
equivalence are shared. Every `r.` operation is a conversion to the `c.` (coding transcript) or
`n.` (non-coding) spelling and back. *HGVS agrees* on the numbering. *VariantValidator* accepts
`r.` in exons but rejects `r.*10`, intronic `r.`, `r.spl` and `r.0`, all of which HGVS allows and
weaver accepts.
Spec: [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "nucleotide `r.123`
relates to `c.123` or `n.123`"; [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/)
for `r.0`, `r.spl`, `r.?`, `r.=`.
Tests: `transcript_coordinates_test::mitochondrial_variants_share_the_genomic_implementation`;
`rna_test::r_is_c_in_rna_letters`, `::r_on_a_non_coding_transcript_is_n`,
`::statements_about_the_transcript_round_trip`, `::r_normalises_validates_and_has_alleles_like_c`;
`rna_is_the_transcript_in_other_letters` (property).

## Normalisation

**3' shift, cyclic over repeats.** A deletion, duplication or insertion is shifted as far 3' as
the reference keeps matching it cyclically, on the strand of the sequence the variant is written
on. An insertion that repeats the bases immediately before it is written as a duplication, and an
insertion that slides has its bases rotated by how far it moved. *HGVS, biocommons,
VariantValidator and ClinVar agree.*
Spec: [general, 3' rule](https://hgvs-nomenclature.org/stable/recommendations/general/);
[insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/): "tandem
duplications are described as a duplication, not an insertion".
Tests: `src/normalize.rs::deletion_shifts_3_prime_and_states_no_bases`,
`::insertion_that_repeats_preceding_bases_becomes_duplication`, `::a_slid_insertion_is_rotated`,
`::shift_5_mirrors_shift_3`; `transcript_coordinates_test::normalize_converts_genomic_and_noncoding_insertions_to_duplications`;
`normalising_preserves_the_edited_sequence` and `ambiguous_range_is_sound_and_complete` (properties).

**A delins is not shifted.** Only pure deletions, duplications and insertions slide; a delins with
bases on both sides is left where it was written. Shifting it was wrong on 72 ClinVar rows.
*Biocommons agrees.*
Spec: [deletion-insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) gives
no shift rule for delins; the 3' rule is stated for the equivalent placements of one change.
Tests: `src/normalize.rs::delins_is_not_shifted`, `::substitution_is_left_alone`.

**`del` and `dup` are written bare.** Normalised output writes `c.306del`, not `c.306delC`. Bases
the input stated survive when the edit does not move and are dropped when it does, since they
would then be stale; nothing is filled in. *HGVS and VariantValidator agree; biocommons differs*
(it fills the deleted bases in). There is no option: this is the one spelling.
Spec: [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/): "the
recommendation is not to describe the variant as `g.33344591delA`".
Tests: `src/normalize.rs::deletion_shifts_3_prime_and_states_no_bases`;
`transcript_coordinates_test::mitochondrial_variants_share_the_genomic_implementation` (stated
bases kept when the edit does not move); `tests/test_api.py::test_normalization`.

**A repeat resolves to its whole run.** `c.10AC[3]` means every existing copy of the unit from
`c.10` becomes three copies, so its SPDI and protein reading cover the run. Before this, 0 of 394
ClinVar repeat rows matched on SPDI; after, 386. *ClinVar agrees* on those rows.
Spec: [repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/):
the bracketed number "represents the total count of repeat units".
Tests: `transcript_coordinates_test::repeat_on_the_minus_strand_projects_to_its_whole_run`;
`protein_allele_test::substitution_insertion_delins_dup_and_repeat_resolve_on_the_protein`;
`equivalence_cases_test::test_analogous_repeat_equivalence`, `::test_multi_unit_repeat_equivalence`.

## Protein consequences

Protein consequences are read from codons, not from a diff of two protein strings. The rules that
follow all bear on where the stop is.

**The declared CDS end is the stop.** The reference protein ends where the transcript record
declares the CDS ends, when that codon really is a stop; only a loosely marked CDS end falls back
to the first stop in the translation. A selenocysteine `TGA` inside the CDS is therefore not a
stop: SEPN1 `NM_020451.2:c.943G>A` is `p.(Gly315Ser)`. *Biocommons differs* (`p.?`, since it sees
two in-frame stops). *Not checked* against VariantValidator.
Spec: none; the nomenclature assumes the reference protein is known.
Tests: `src/protein.rs::a_selenocysteine_codon_is_not_the_stop`;
`real_transcripts_test::biocommons_real_transcript_cases` (MULTISTOP01, a recorded difference).

**A frameshift whose first changed residue is the stop is an extension.** `Ter902ArgextTer92`,
not `Ter902ArgfsTer93`. *HGVS and biocommons agree*; five ClinVar rows moved from mismatch to
Analogous.
Spec: [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/): "the
variant extends the amino acid sequence at the C-terminal end and is therefore by definition an
extension".
Tests: `src/protein.rs::a_frameshift_in_the_stop_codon_is_an_extension`,
`::stop_loss_is_an_extension_to_the_next_stop`; `real_transcripts_test::biocommons_real_transcript_cases`
(EXT01 to EXT06).

**An extension needs the stop codon itself to change.** An in-frame insertion just before the stop
that leaves the stop intact is an insertion, not an extension, even though the stop's position is
the first residue that differs.
Spec: [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/), the
same definition: the stop codon is what changes.
Tests: `src/protein.rs::an_insertion_before_the_stop_that_repeats_the_last_residue_is_not_an_extension`,
`::an_in_frame_deletion_reaching_the_stop_keeps_it_original`.

**A stop formed inside inserted bases is a delins ending in Ter**, not a frameshift: translation
ends before the new frame reads a single reference base.
Spec: [frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/):
"variants which introduce an immediate translation termination (stop) codon are described as
nonsense variant", not a frameshift.
Tests: `src/protein.rs::a_stop_inside_the_inserted_bases_is_a_delins_not_a_frameshift`;
`equivalence_cases_test::test_immediate_stop_normalization`.

**In-frame changes are written 3'-most.** The shared tail is trimmed, so an in-frame deletion or
insertion is written at its 3'-most equivalent residues. *HGVS agrees.* *ClinVar differs* in
spelling: it writes `Xxx_Yyyins…` at the 3' end of a run where weaver writes `dup`, and one-letter
repeat forms such as `p.490PRS[1]`; 170 rows in 100,000 are the same protein in another spelling,
and equivalence judges them Analogous.
Spec: [protein deletion](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/):
"the most C-terminal position possible of the reference sequence is arbitrarily assigned to have
been changed".
Tests: `src/protein.rs::in_frame_deletion_and_duplication_take_the_3_prime_position`;
`protein_insertions_are_written_3_prime_most` and `in_frame_indels_describe_the_translated_protein`
(properties).

**A start-codon change is written specifically.** `c.1A>G` is `p.(Met1Val)`. *HGVS and
VariantValidator prefer `p.Met1?`*; weaver keeps the specific prediction because it carries more
information, and `VariantTransformSettings(start_codon=HgvsQuestion)` rewrites it on request.
Spec: [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/):
`p.Met1?` when "the consequence, on the protein level, of a variant affecting the translation
initiation codon can not be predicted". A recorded difference.
Tests: `mapping_test::test_mapper_c_to_p_start_codon_subst`;
`real_transcripts_test::biocommons_real_transcript_cases` (INITMET01, recorded);
`src/transform.rs` unit tests for the rewrite.

**An edit across the CDS start.** A deletion that removes the whole CDS is `p.0?`; an edit that
starts in the 5'UTR and reaches into the CDS disrupts the start codon and is `p.Met1?`; an edit
entirely upstream is `p.?`. *ClinVar differs*: it writes `Met1fs`, `Met1_Glu2insGly…`, committing
to a consequence; weaver does not commit because the initiation site is not predictable.
*Biocommons agrees* on `p.Met1?` (its own test table). For a substitution entirely in the 5'UTR
*VariantValidator writes `p.(=)`* where weaver writes `p.?`; both are defensible, and weaver's is
the more honest since a 5'UTR change can create an upstream start.
Spec: [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/):
`p.0?` "when you predict that no protein is produced"; `p.Met1?` as above; `p.Met1ext-5` shows a
5'UTR change can activate an upstream initiation site.
Tests: `real_transcripts_test::biocommons_real_transcript_cases` (WHOLEGENE01/02, INITMET02/03);
`protein_allele_from_coding_test::statements_have_no_allele_and_a_wrong_protein_is_an_error`,
`::a_coding_variant_agrees_with_every_spelling_of_its_consequence` (`c.-3_2del`). The pure-5'UTR
`p.?` has no direct test.

**Frameshift length.** `fsTer N` counts to the first stop the new frame reaches; an immediate stop
is a plain nonsense substitution. *HGVS agrees.*
Spec: [frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/): the
stop position is counted "starting from the first changed amino acid as codon 1"; "the shortest
frameshift variant possible contains `fsTer2`".
Tests: `src/protein.rs::frameshift_reports_first_changed_residue_and_distance_to_stop`,
`::substitutions_nonsense_and_silence`; `snv_effect_follows_the_codon_table` (property).

**Protein alleles.** From a `p.` description, only edits that name a sequence have an allele:
frameshift, extension, `p.?`, `p.0?` and `p.Met1?` are unsupported. From a coding variant, every
consequence has one (`protein_allele`, `protein_vrs`), because the edited transcript's translation
is known: a frameshift's allele is the residues from the first change to the new stop. The
translated CDS must be the protein the provider serves; a difference is an annotation error and
is reported as one. A silent change is the reference allele over the codons the edit touched; a
deleted CDS is the deletion of the whole protein.
Spec: VRS 2.0.1 `Allele` on a protein `SequenceReference`; none in HGVS.
Tests: `protein_allele_test::consequences_have_no_allele`;
`protein_allele_from_coding_test::consequences_without_a_p_sequence_still_have_an_allele`,
`::definite_changes_give_the_same_allele_by_either_route`,
`::statements_have_no_allele_and_a_wrong_protein_is_an_error`;
`the_protein_allele_is_the_edited_translation` (property).

## Validation and stated bases

**Stated bases are checked by `validate`, and nowhere else.** A variant's stated reference
(`c.123A>G`'s A) is compared with the sequence by `validate`, and by `c_to_p`, which errors on a
mismatch. Canonical alleles, SPDI and VRS use the sequence's own bases and ignore what was stated,
for nucleotides and for protein residues alike, because the allele describes the sequence, not
the description. *VRS agrees.*
Spec: VRS 2.0.1: an `Allele` is a `SequenceLocation` and a state; the reference is the sequence's.
Tests: `src/allele.rs::a_stated_reference_that_disagrees_with_the_sequence_is_ignored`;
`transcript_coordinates_test::validate_checks_stated_reference_through_transcript_coordinates`;
`protein_allele_test::validation_checks_the_named_and_stated_residues`;
`projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
(the `validate` cases).

## Equivalence

**Judged by allele and by the protein left behind.** Two nucleotide variants are the same change
exactly when their canonical alleles are equal. A coding variant and a protein description agree
when the protein the variant leaves is the protein the description leaves; two descriptions agree
when they leave the same protein. `p.Tyr165Ter`, `p.Ala164_Tyr165insTer` and the `c.` deletion
behind them are one change. `Identity` is the same text, or a prediction written exactly as given.
[How it decides](equivalence_logic.md) has the full rules. *Biocommons differs*: it compares
normalised text.
Spec: none; HGVS describes, it does not compare.
Tests: `transcript_coordinates_test::canonical_alleles_make_spdi_vrs_and_equivalence_one_value`;
`equivalence_cases_test::test_clinvar_regression_tyr165ter`, `::test_analogous_fs_wildcard_unification`;
`protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`;
`a_variant_agrees_with_its_projection` (property).

**A description that says nothing matches nothing.** `p.?` matches no other description.
`p.Met1?` says only where the change starts and matches a description whose change starts there;
it does not match ClinVar's `Met1fs`, a commitment weaver did not make. Two ClinVar rows judge
Different for this reason, on purpose.
Spec: [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/)
defines `p.?` and `p.Met1?` as statements of what is not known.
Tests: `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`
(`c.-3_2del` against `p.Leu3fs` is Different).

**Judging with no protein sequence is an error, not a guess.** Earlier versions reconciled the
residues two descriptions happened to name; now the sequence is required.
Spec: none.
Tests: `equivalence_test::equivalence_needs_the_sequence`.

**Versions of one protein accession are compared on ours.** ClinVar often names an older `NP_`
version; the descriptions are read against the protein the provider serves.
Spec: none.
Tests: `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`
(`NP_X.1` against `NP_X.0`).

**A cis allele compares as a set.** `c.[a;b]` and `c.[b;a]` are the same change; a single-member
`c.[a]` compares as `a`; a cis allele never compares equal to a single variant of two or more
members.
Spec: [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/): variants
"should be listed in genomic order", which is a writing convention, not a different allele.
Tests: `cis_phased_test::cis_alleles_compare_as_sets_of_their_members_alleles`.

## Alleles, SPDI and VRS

**The canonical allele is fully justified.** A change is trimmed to what it alters and then
widened over the whole region in which it could equally be written, so every spelling of one
change is one allele with one identifier. *VRS and VOCA agree*; this is their normalisation.
`to_spdi_unambiguous` renders it. The plain `to_spdi` renders the 3'-normalised variant instead,
with an insertion placed at its second flank as SPDI counts.
Spec: VRS 2.0.1 "Normalization" (fully justified alleles); refget/SPDI interbase coordinates.
Tests: `src/allele.rs::deletion_anywhere_in_a_run_is_the_same_allele`,
`::duplication_and_insertion_of_the_unit_are_the_same_allele`,
`::literal_insertion_stays_literal_and_substitution_is_trimmed`;
`transcript_coordinates_test::plain_spdi_places_an_insertion_at_its_second_flank`;
`ambiguous_range_is_sound_and_complete` and `vrs_and_spdi_round_trip` (properties).

**Refget accessions are computed over the normalised sequence.** Letters uppercased, everything
else dropped, as the refget specification defines; chr19 of the NCBI FASTA, uppercased, hashes to
the specification's published `SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl`. Without a `Refget` lookup the
accession is computed from the whole sequence; `from_vrs` then needs the accession passed. The
digest is checked against the sequence either way.
Spec: [refget](https://samtools.github.io/hts-specs/refget.html), "Sequence Normalization" and
the `ga4gh` identifier.
Tests: `src/vrs.rs::refget_accession_is_over_the_normalised_sequence`,
`::digest_and_identifier_match_the_spec_example`;
`tests/test_digest_table.py::test_digest_follows_the_refget_normalisation_rule`;
`vrs_round_trip_test::the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked`.

**Uncertain breakpoints.** A deletion with uncertain breakpoints, `g.(?_100)_(200_?)del`, renders
as an Allele with Range bounds and an empty state, unnormalised; a duplication with uncertain
breakpoints has no Allele (its bases are unknown) and renders as a `CopyNumberChange` gain. An
unknown bound prints as `?`; a parenthesised exact position is exact.
Spec: [uncertain positions](https://hgvs-nomenclature.org/stable/recommendations/uncertain/):
`(A_B)_(C_D)` "where `B_C` describes the minimal extent and `A_D` the maximal", `?` for unknown;
VRS 2.0.1 `Range` and `CopyNumberChange`.
Tests: `vrs_range_test::uncertain_breakpoints_become_ranges`,
`::an_imprecise_deletion_has_an_empty_literal_state_and_json_ranges`,
`::only_deletions_may_have_uncertain_breakpoints`, `::unknown_bounds_round_trip_as_question_marks`;
`vrs_copy_change_test::an_imprecise_duplication_is_a_gain_over_range_bounds`,
`::to_vrs_variation_gives_a_change_for_an_imprecise_dup_only`;
`spec_summary_test::test_uncertain_intervals`.

**`copyChange` is a label.** VRS 2.0.1 defines `copyChange` as a string enum (`gain`, `loss`, …);
2.0.0 used EFO codes in a `MappableConcept`. weaver writes labels and reads labels, bare EFO
CURIEs and the 2.0.0 object form. The specification's published `CopyNumberChange` example digest
only reproduces when `copyChange` is digested as the pre-release string `EFO:0030071`; both digests
are pinned in tests so the machinery is anchored to a published value while output follows the
current schema. EFO's high-level loss is `EFO:0020073`, not in the `00300xx` block.
Spec: VRS 2.0.1 `vrs-source.yaml`, `CopyNumberChange.properties.copyChange.enum`.
Tests: `src/vrs.rs::copy_number_change_digests_the_label_over_the_location`,
`::copy_change_terms_are_read_by_label_or_efo_code`;
`vrs_copy_change_test::the_whole_gain_family_is_a_dup_and_the_loss_family_a_del`.

**`CisPhasedBlock` members are sorted before digesting**, as the digest-serialisation rule requires,
so the identifier does not depend on the order written; the `members` array keeps that order. The
schema's `minItems: 2` is relaxed to one because HGVS allows the degenerate `c.[145C>T]`. The
in-trans form `[..];[..]` is rejected: it describes two molecules.
Spec: VRS 2.0.1 computed identifiers, "order arrays of digests and ids by Unicode Character Set
values"; [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) for cis
`[a;b]` against trans `[a];[b]`.
Tests: `src/vrs.rs::cis_phased_block_digest_matches_the_spec_validation_data`,
`::a_cis_phased_block_carries_its_members_in_the_order_given`;
`cis_phased_test::the_identifier_is_the_same_whichever_way_the_members_are_written`,
`::alleles_in_trans_are_two_molecules_and_are_refused`.

**Reading back.** `from_vrs` and `from_spdi` give the 3'-normalised variant on the allele's own
sequence. An insertion trimmed to before the first base becomes a delins of that base, since HGVS
has no insertion before base 1; a replacement that is the reverse complement of what it replaces
reads back as `inv`; a `LengthExpression` comes back as `insN[n]`; a Range of copies is refused,
since HGVS has no spelling for it.
Spec: [numbering](https://hgvs-nomenclature.org/stable/background/numbering/) (no position before
1); [inversion](https://hgvs-nomenclature.org/stable/recommendations/DNA/inversion/);
[insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/) for `insN[n]`.
Tests: `vrs_round_trip_test::nucleotide_variants_round_trip_to_their_normalised_form`,
`::alleles_from_other_producers_parse`;
`vrs_length_test::length_expressions_read_back_as_the_recommended_spelling`;
`vrs_copy_number_test::from_vrs_checks_the_count_it_is_given`; `vrs_and_spdi_round_trip` (property).

**What is not rendered.** `Adjacency`, `Terminus` and `DerivativeMolecule` describe breakends and
fusions this grammar has no input for.
Spec: VRS 2.0.1; no HGVS input in this grammar.
Tests: none.

## Parsing

**The grammar is checked against biocommons's table**, 580 inputs over the 92 rules the two
grammars share, and agrees on all but one: weaver accepts a terminator inside an amino acid
sequence (`insTerGlu`), which *biocommons rejects* and *ClinVar writes*.
Spec: the nomenclature pages above; the table is biocommons's reading of them.
Tests: `grammar_test::biocommons_grammar_table` (the difference is listed in `KNOWN_DIFFERENCES`);
`parser_round_trips_canonical_hgvs` and `parser_never_panics` (properties).

**Forms accepted beyond biocommons.** Alleles in cis for every coordinate system; `r.(123a>g)`
predicted RNA changes and the statements `r.0`, `r.?`, `r.spl`, `r.=`; `p.0`, `p.0?`, `p.Met1?`;
uncertain genomic breakpoints with `?`; insertions of a stated length, `insN[20]`, `insN[(20_30)]`,
`delinsN[12]` and the older `ins(20)`; `copyN`; `r.` positions with `*`.
Spec: [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/);
[RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/);
[protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/);
[uncertain positions](https://hgvs-nomenclature.org/stable/recommendations/uncertain/);
[insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/).
Tests: `cis_phased_test::cis_alleles_parse_and_print_back_in_every_coordinate_system`;
`rna_test::predicted_changes_keep_their_parentheses`, `::statements_about_the_transcript_round_trip`;
`spec_summary_test::test_uncertain_intervals`, `::test_spec_summary_variants`;
`vrs_length_test::the_older_spelling_gives_the_same_allele`, `::an_uncertain_length_is_a_range`.

**Recommended spellings on output.** `insN[20]` rather than `ins(20)`; `extTer8` rather than `ext*8`
(*VariantValidator writes `ext*`*; the two are compared as equal); bare `del` and `dup`.
Spec: [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/);
[extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/), where
"both notations are acceptable"; [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/).
Tests: `vrs_length_test::length_expressions_read_back_as_the_recommended_spelling`;
`real_transcripts_test::biocommons_real_transcript_cases` (`ext*` compared equal to `extTer`).

## The data contract

**A range past the end of a sequence returns the bases that exist.** `get_seq` must return what is
there rather than error or return nothing, because the core pages through sequences in fixed
blocks and a short final block is how it learns where the sequence ends.
Spec: none; weaver's own interface.
Tests: `reference_store_agrees_with_direct_slicing` (property);
`tests/test_refget.py::test_sequences_come_from_the_server_with_ranges_clamped`.

**Interval methods are half-open and 0-based**; HGVS positions are 1-based and inclusive, with
`c.` skipping the non-existent position 0. Exon `reference_end` is inclusive.
Spec: [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "there is no
nucleotide `c.0`".
Tests: `transcript_coordinates_round_trip` and `intron_offsets_agree_from_either_exon` (properties);
`transcript_coordinates_test::spdi_interval_honours_cds_end_anchor`,
`::spdi_interval_honours_non_zero_cds_start`, `::spdi_interval_on_minus_strand`.

**A mapper keeps its cache.** A `VariantMapper` holds sequence blocks and refget accessions for as
long as it lives. Build one and reuse it.
Spec: none.
Tests: `tests/test_mapper_cache.py::test_sequences_are_fetched_once_per_mapper`,
`::test_a_refget_lookup_is_asked_before_the_sequence_is_hashed`.

**Refget is its own seam.** A sequence source need know nothing about digests; the `Refget` lookup
is a separate object, and a provider without one still works for everything but naming the
sequence behind a digest.
Spec: [refget](https://samtools.github.io/hts-specs/refget.html).
Tests: `vrs_round_trip_test::the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked`;
`tests/test_refget.py::test_refget_lookups_both_ways`;
`tests/test_digest_table.py::test_the_table_agrees_with_what_the_mapper_computes_and_names_the_sequence_back`.
