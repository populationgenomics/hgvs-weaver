# The Choices Weaver Makes

HGVS leaves room for judgement, and the tools that read it have settled on different answers. This
page lists every deliberate choice weaver makes where another tool could reasonably answer
otherwise: what weaver does, why, a worked example, the recommendation it rests on, and what other
tools do.

## How to read this page

Every choice has the same shape.

| Part | What it holds |
| :--- | :--- |
| **Example** | Input on the left, weaver's output on the right, taken from the test that asserts it. `→` reads "becomes". |
| **Agrees** / **Differs** | Only where the other tool was actually checked: by running it (VariantValidator's REST API, the biocommons `hgvs` grammar and regression tables, ClinVar's 100,000-variant sample) or by reading its code. Anything else is marked *not checked*. |
| **Spec** | The [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) page (or the VRS or refget specification) the choice rests on, with a verdict: *agrees* when the specification says or clearly implies the same, *differs* when weaver departs from what it recommends, *silent* when it says nothing. Departures are few and each says why. |
| **Tests** | The tests that assert it, as `file::function`. Files are under `hgvs-weaver/tests/` unless written as `src/…` (unit tests) or `tests/…py` (Python); property tests live in `hgvs-weaver/tests/properties/main.rs`. |

The tools referred to: **biocommons** is the Python `hgvs` package weaver was originally a port
of; **VariantValidator** is the web service; **ClinVar** is the descriptions in its variant
summary; **VRS** is the GA4GH Variation Representation Specification 2.0.1.

## At a glance

| Choice | Spec | biocommons | VariantValidator | ClinVar |
| :--- | :---: | :---: | :---: | :---: |
| [A projection states the bases of the target](#a-projection-states-the-bases-of-the-target) | agrees | agrees | agrees | – |
| [A change across a splice junction has no genomic form](#a-change-across-a-splice-junction-has-no-genomic-form) | agrees | – | differs | – |
| [No normalisation before projecting](#no-normalisation-before-projecting) | – | – | differs | – |
| [Gapped exons are projected locally](#gapped-exons-are-projected-locally) | – | – | differs | – |
| [A repeat is projected as its whole run](#a-repeat-is-projected-as-its-whole-run) | agrees | – | – | – |
| [A position outside every exon is an error](#a-position-outside-every-exon-is-an-error) | agrees | differs | – | – |
| [Intronic positions are carried as given](#intronic-positions-are-carried-as-given) | agrees | – | – | – |
| [Mitochondrial and RNA descriptions share the DNA machinery](#mitochondrial-and-rna-descriptions-share-the-dna-machinery) | agrees | – | differs | – |
| [Shift 3', cyclically over repeats](#shift-3-cyclically-over-repeats) | agrees | agrees | agrees | agrees |
| [A delins is not shifted](#a-delins-is-not-shifted) | – | agrees | – | – |
| [Deletions and duplications are written bare](#deletions-and-duplications-are-written-bare) | agrees | differs | agrees | – |
| [A repeat resolves to its whole run](#a-repeat-resolves-to-its-whole-run) | agrees | – | – | agrees |
| [The declared CDS end is the stop](#the-declared-cds-end-is-the-stop) | – | differs | – | – |
| [A frameshift that starts at the stop is an extension](#a-frameshift-that-starts-at-the-stop-is-an-extension) | agrees | agrees | – | agrees |
| [An extension needs the stop codon itself to change](#an-extension-needs-the-stop-codon-itself-to-change) | agrees | – | – | – |
| [A stop formed inside inserted bases is a delins ending in Ter](#a-stop-formed-inside-inserted-bases-is-a-delins-ending-in-ter) | agrees | – | – | – |
| [In-frame changes are written 3'-most](#in-frame-changes-are-written-3-most) | agrees | – | – | spelling differs |
| [A start-codon change is written specifically](#a-start-codon-change-is-written-specifically) | **differs** | differs | differs | – |
| [Edits outside the CDS are statements](#edits-outside-the-cds-are-statements) | agrees | agrees | differs | differs |
| [Frameshift length counts to the first new stop](#frameshift-length-counts-to-the-first-new-stop) | agrees | – | – | – |
| [Protein alleles come from the coding change](#protein-alleles-come-from-the-coding-change) | agrees | – | – | – |
| [Stated bases are checked by validate and nowhere else](#stated-bases-are-checked-by-validate-and-nowhere-else) | agrees | – | – | – |
| [Judged by allele and by the protein left behind](#judged-by-allele-and-by-the-protein-left-behind) | – | differs | – | – |
| [A description that says nothing matches nothing](#a-description-that-says-nothing-matches-nothing) | agrees | – | – | differs |
| [Judging with no protein sequence is an error](#judging-with-no-protein-sequence-is-an-error) | – | – | – | – |
| [Versions of one protein accession compare on ours](#versions-of-one-protein-accession-compare-on-ours) | – | – | – | – |
| [A cis allele compares as a set](#a-cis-allele-compares-as-a-set) | agrees | – | – | – |
| [The canonical allele is fully justified](#the-canonical-allele-is-fully-justified) | agrees | – | – | – |
| [Refget accessions are computed over the normalised sequence](#refget-accessions-are-computed-over-the-normalised-sequence) | agrees | – | – | – |
| [Uncertain breakpoints become Range bounds](#uncertain-breakpoints-become-range-bounds) | agrees | – | – | – |
| [copyChange is a label](#copychange-is-a-label) | agrees | – | – | – |
| [CisPhasedBlock members are sorted before digesting](#cisphasedblock-members-are-sorted-before-digesting) | **differs** | – | – | – |
| [Reading back gives the normalised variant](#reading-back-gives-the-normalised-variant) | agrees | – | – | – |
| [Breakends and fusions are not rendered](#breakends-and-fusions-are-not-rendered) | – | – | – | – |
| [The grammar is checked against the biocommons table](#the-grammar-is-checked-against-the-biocommons-table) | – | one difference | – | agrees |
| [Forms accepted beyond biocommons](#forms-accepted-beyond-biocommons) | agrees | differs | – | – |
| [Recommended spellings on output](#recommended-spellings-on-output) | agrees | – | differs | – |
| [A range written backwards is refused when parsed](#a-range-written-backwards-is-refused-when-parsed) | agrees | differs | agrees | – |
| [A range past the end of a sequence returns the bases that exist](#a-range-past-the-end-of-a-sequence-returns-the-bases-that-exist) | – | – | – | – |
| [Interval methods are half-open and 0-based](#interval-methods-are-half-open-and-0-based) | agrees | – | – | – |
| [A mapper keeps its cache](#a-mapper-keeps-its-cache) | – | – | – | – |
| [Refget is its own seam](#refget-is-its-own-seam) | agrees | – | – | – |

**Spec** is the HGVS nomenclature, or VRS or refget where the choice is about them; a dash means
it is silent. For a tool, a dash means it was not checked or the choice does not apply to it.
The two departures from a specification are in bold.

## Projection between sequences

### A projection states the bases of the target

**Where record and genome differ, the projected edit is written against the target's own
bases, as the minimal edit that produces the same alternate.**

A transcript record and its genome can differ at a base: RefSeq transcripts are curated against
submitted mRNAs. When a variant is projected, the edit is resolved on the source to the bases it
removes and the bases it puts there. Where the target holds different bases over the projected
range, weaver writes the edit that turns the target's bases into that alternate, collapsing to `=`
when they already coincide. A duplication or inversion carries the source's bases, because that is
what the molecule becomes. A deletion is positional. Where record and genome agree, which is nearly
everywhere, the edit is carried across unchanged.

**Example** (MUC2 has a record C over a genomic G; SHANK3 a record C over a genomic T):

```text
NM_002457.5:c.12468C>A       →  NC_000011.10:g.1099802G>A        the genome's G, not the record's C
NM_002457.5:c.12468dup       →  NC_000011.10:g.1099802delinsCC   a dup carries the record's C, twice
NM_002457.5:c.12468del       →  NC_000011.10:g.1099802del        a deletion is positional

NM_001372044.2:c.1568=       →  NC_000022.11:g.50697558T>C       no change on the record is a change on the genome
NC_000022.11:g.50697558T>C   →  NM_001372044.2:c.1568=           and the other way round
NM_001372044.2:c.1568C>T     →  NC_000022.11:g.50697558=
```

- **Agrees:** biocommons (`replace_reference`); VariantValidator on substitution, identity,
  duplication and deletion at the MUC2 base.
- **Spec (agrees):** [general](https://hgvs-nomenclature.org/stable/recommendations/general/):
  "descriptions on RNA/protein level should describe the changes observed on that level"; a
  description is of the sequence it is on.
- **Tests:** `projection_reference_test::a_projection_states_the_genomes_bases_not_the_records`,
  `::a_projection_states_the_transcripts_bases_not_the_genomes`,
  `::the_re_read_is_in_the_targets_orientation_on_the_minus_strand`;
  `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
  (MUC2 and SHANK3 with NCBI's alignments); `a_variant_agrees_with_its_projection` (property, where
  record and genome agree).

### A change across a splice junction has no genomic form

**An `r.` change spanning an exon boundary describes the spliced RNA, and weaver refuses to
project it to the genome.**

The `c.` spelling of the same change still projects, to a deletion that includes the intron, and
its protein consequence is still predicted. weaver refuses the `r.` projection because a warning is
easy to miss and the two molecules are not the same.

**Example** (exon 1 of `NM_R.1` ends at r.45):

```text
NM_R.1:r.44_47del   → genome    UnsupportedOperation: "spans a splice junction: the spliced RNA has no single genomic equivalent"
NM_R.1:r.44_47del   → protein   the same frameshift as NM_R.1:c.44_47del
NM_R.1:c.44_47del   → genome    projects, over the intron
NM_R.1:r.45del      → genome    projects: one exon
```

- **Differs:** VariantValidator projects the `c.` reading with the warning "spans at least one
  intron".
- **Spec (agrees):** [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/)
  and [RNA splicing](https://hgvs-nomenclature.org/stable/recommendations/RNA/splicing/), where exon
  skipping is written as an `r.` deletion across the junction (`r.(3277_3432del)`), an RNA event.
- **Tests:** `rna_test::r_projects_to_the_genome_within_one_exon_only`,
  `::r_predicts_the_protein_like_c`; `rna_across_a_junction_has_a_protein_but_no_genomic_form`
  (property).

### No normalisation before projecting

**weaver projects the variant where it is written; shifting it is the caller's decision.**

Normalisation is available separately, as `normalize_variant`.

**Example** (SHANK3, record C over a genomic T):

```text
NM_001372044.2:c.1568dup  →  NC_000022.11:g.50697558delinsCC      weaver, at the position given
                          →  c.1569dup first, then the genome      VariantValidator, reporting the automap
```

- **Differs:** VariantValidator 3'-shifts on the transcript first and reports the automapping.
- **Spec (silent):** none; the [3' rule](https://hgvs-nomenclature.org/stable/recommendations/general/)
  says how to write a variant, not that a tool must rewrite what it is given.
- **Tests:** `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
  (positions come back as given).

### Gapped exons are projected locally

**Where an exon aligns to the genome with insertions and deletions, weaver projects through the
exon's CIGAR and writes a local edit.**

**Example** (SHANK3, whose exon alignment is heavily gapped):

```text
NM_001372044.2:c.1568dup  →  NC_000022.11:g.50697558delinsCC                     weaver
                          →  NC_000022.11:g.50695049_50697558delins… (2.5 kb)    VariantValidator
```

- **Differs:** VariantValidator falls back to a whole-region delins on heavily gapped exons.
- **Spec (silent):** none.
- **Tests:** `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
  (the SHANK3 dup case); `src/transcript_mapper.rs::test_g_to_n_cigar`.

### A repeat is projected as its whole run

**A repeat is widened to the full run of its unit on the source before projecting, so the other
strand reads it from the right end.**

Without this the genomic repeat started mid-run.

**Example** (`NC_REP.1` is `TTTTT GCCATT GCCATT GCCATT AAAAA…`; `NM_REP_MINUS.1` reads it on the
minus strand, where the run is `AATGGC` from c.78):

```text
NM_REP_MINUS.1:c.78AATGGC[4]  →  NC_REP.1:g.6_23GCCATT[4]           the whole run, in plus-strand orientation
NC_REP.1:g.6GCCATT[4]         →  NM_REP_MINUS.1:c.78_95AATGGC[4]
```

- **Others:** not checked.
- **Spec (agrees):** [repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/):
  the count is the total number of units in the run.
- **Tests:** `transcript_coordinates_test::repeat_on_the_minus_strand_projects_to_its_whole_run`.

### A position outside every exon is an error

**Transcript positions resolve through the exon structure only; a position past the last exon is
a `ValidationError`, not an extrapolated coordinate.**

A silently extrapolated coordinate is wrong without saying so.

**Example** (`NM_PLUS10.1` is one exon of 100 bases; c.*60 is its last base):

```text
NM_PLUS10.1:c.*60T>A   →  projects
NM_PLUS10.1:c.*61A>G   →  ValidationError, from c_to_g, validate and to_spdi_unambiguous alike
```

- **Differs:** biocommons extrapolates in some paths.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "it is not
  allowed to describe variants in nucleotides beyond the boundaries of a reference sequence".
- **Tests:** `transcript_coordinates_test::a_position_in_no_exon_is_an_error_not_an_extrapolation`;
  `transcript_coordinates_round_trip` (property) covers every in-exon position.

### Intronic positions are carried as given

**An intronic edit has no transcript base to compare against, so it is kept as written on the
transcript; on the genome its bases are checked like any other. Normalisation leaves it where it
is.**

**Example** (`NM_SPLICED.1`, whose intron is genome[50..60) on `NC_D.1`, where g.53 is A):

```text
NC_D.1:g.53A>G          →  NM_SPLICED.1:c.40+3A>G      carried as given
NM_SPLICED.1:c.40+3C>G  →  NC_D.1:g.53A>G              re-read against the genome, which has A

normalize NM_PLUS10.1:c.30_30+1insA   →  NM_PLUS10.1:c.30_30+1insA    unchanged, not an error
normalize NM_MINUS10.1:c.1+5del       →  NM_MINUS10.1:c.1+5del
```

- **Others:** not checked.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): RNA and coding
  reference sequences "do not contain intron sequences and can therefore not be used to describe
  variants affecting these sequences".
- **Tests:** `projection_reference_test::an_intronic_position_has_no_transcript_base_to_re_read`;
  `transcript_coordinates_test::normalize_leaves_intronic_coding_variants_as_written`;
  `intron_offsets_agree_from_either_exon` (property).

### Mitochondrial and RNA descriptions share the DNA machinery

**`m.` is `g.` on the mitochondrial reference; `r.` is `c.` or `n.` in RNA letters.**

`m.` and `g.` differ only in the letter they are written with: normalisation, SPDI, validation and
equivalence are shared. Every `r.` operation is a conversion to the `c.` (coding transcript) or
`n.` (non-coding) spelling and back; the statements `r.0`, `r.spl`, `r.?` and `r.(=)` map to
`p.0`, `p.?`, `p.?` and `p.(=)`.

**Example:**

```text
normalize NC_TEST.1:m.1012_1013insT   →  NC_TEST.1:m.1012dup      as g. would, written back as m.

NM_R.1:r.-3_5delinsuu   ⇄  NM_R.1:c.-3_5delinsTT
NR_R.1:r.10c>g          ⇄  NR_R.1:n.10C>G
NM_R.1:r.(10c>g)        →  NM_R.1:c.10C>G                        the prediction flag is dropped on the way
NM_R.1:r.spl            →  NP_R.1:p.?
NM_R.1:r.0              →  NP_R.1:p.0
```

- **Agrees:** HGVS, on the numbering.
- **Differs:** VariantValidator accepts `r.` in exons but rejects `r.*10`, intronic `r.`, `r.spl`
  and `r.0`, all of which HGVS allows.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "nucleotide
  `r.123` relates to `c.123` or `n.123`";
  [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/) for
  `r.0`, `r.spl`, `r.?`, `r.=`.
- **Tests:** `transcript_coordinates_test::mitochondrial_variants_share_the_genomic_implementation`;
  `rna_test::r_is_c_in_rna_letters`, `::r_on_a_non_coding_transcript_is_n`,
  `::statements_about_the_transcript_round_trip`, `::r_predicts_the_protein_like_c`,
  `::r_normalises_validates_and_has_alleles_like_c`; `rna_is_the_transcript_in_other_letters`
  (property).

## Normalisation

### Shift 3', cyclically over repeats

**A deletion, duplication or insertion is shifted as far 3' as the reference keeps matching it
cyclically, on the strand of the sequence the variant is written on.**

An insertion that repeats the bases immediately before it becomes a duplication, and an insertion
that slides has its bases rotated by how far it moved.

**Example** (on `ACGT` repeated, and on the ten bases `TTCAGCAGTT`):

```text
NC_TEST.1:g.1012_1013insT      →  NC_TEST.1:g.1012dup            repeats the T before it
NC_TEST.1:g.1012_1013insACGT   →  NC_TEST.1:g.1997_2000dup       slides to the end of the repeat
NC_TEST.1:g.1012_1013insAA     →  NC_TEST.1:g.1013_1014insAA     slides one base past an A, still an insertion
NM_PLUS0.1:n.12_13insT         →  NM_PLUS0.1:n.12dup             n. and g. alike

on TTCAGCAGTT:  g.3_5del       →  g.6_8del                       the first CAG shifts onto the second
on GGGGGGGGG:   g.2_3insGA     →  g.3_4insAG                     slid one base, so rotated by one
```

- **Agrees:** HGVS, biocommons, VariantValidator and ClinVar.
- **Spec (agrees):** [general, 3' rule](https://hgvs-nomenclature.org/stable/recommendations/general/);
  [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/): "tandem
  duplications are described as a duplication, not an insertion".
- **Tests:** `src/normalize.rs::deletion_shifts_3_prime_and_states_no_bases`,
  `::insertion_that_repeats_preceding_bases_becomes_duplication`, `::a_slid_insertion_is_rotated`,
  `::shift_5_mirrors_shift_3` (the last two rows above, read from the 0-based ranges those tests
  assert); `transcript_coordinates_test::normalize_converts_genomic_and_noncoding_insertions_to_duplications`;
  `normalising_preserves_the_edited_sequence` and `ambiguous_range_is_sound_and_complete`
  (properties).

### A delins is not shifted

**Only pure deletions, duplications and insertions slide; a delins with bases on both sides stays
where it was written.**

Shifting it was wrong on 72 ClinVar rows.

**Example** (on `TTCAGCAGTT`; the base after the range equals the first deleted base, but sliding
would put the new bases after it, a different sequence):

```text
g.3_5delinsTT   →  g.3_5delinsTT
g.3C>T          →  g.3C>T
```

- **Agrees:** biocommons.
- **Spec (silent):** [deletion-insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/)
  gives no shift rule for delins; the 3' rule is stated for the equivalent placements of one
  change.
- **Tests:** `src/normalize.rs::delins_is_not_shifted`, `::substitution_is_left_alone` (read from
  the 0-based ranges they assert).

### Deletions and duplications are written bare

**Normalised output writes `c.306del`, not `c.306delC`. Stated bases survive when the edit does
not move and are dropped when it does; nothing is filled in.**

Dropped bases would be stale, and the reference has them anyway. There is no option: this is the
one spelling.

**Example:**

```text
NM_TEST:c.4_5del               →  NM_TEST:c.5_6del               bare, shifted
NC_TEST.1:m.1008_1010delTAC    →  NC_TEST.1:m.1008_1010delTAC    unmoved, so the stated bases stay
on TTCAGCAGTT:  g.3_5delCAG    →  g.6_8del                       moved, so they go
```

- **Agrees:** HGVS and VariantValidator.
- **Differs:** biocommons fills the deleted bases in.
- **Spec (agrees):** [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/): "the
  recommendation is not to describe the variant as `g.33344591delA`".
- **Tests:** `src/normalize.rs::deletion_shifts_3_prime_and_states_no_bases`;
  `transcript_coordinates_test::mitochondrial_variants_share_the_genomic_implementation`;
  `tests/test_api.py::test_normalization`.

### A repeat resolves to its whole run

**`c.10AC[3]` means every existing copy of the unit from c.10 becomes three copies, so its allele
and protein reading cover the run.**

Before this, 0 of 394 ClinVar repeat rows matched on SPDI; after, 386.

**Example** (`NC_REP.1` is `TTTTT GCCATT GCCATT GCCATT AAAAA…`; all four are one allele):

```text
NM_REP_MINUS.1:c.78AATGGC[4]
NC_REP.1:g.6GCCATT[4]           →  the same canonical allele: the run of three becomes four
NC_REP.1:g.18_23dup
NC_REP.1:g.23_24insGCCATT
```

- **Agrees:** ClinVar, on those rows.
- **Spec (agrees):** [repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/):
  the bracketed number "represents the total count of repeat units".
- **Tests:** `transcript_coordinates_test::repeat_on_the_minus_strand_projects_to_its_whole_run`;
  `protein_allele_test::substitution_insertion_delins_dup_and_repeat_resolve_on_the_protein`;
  `equivalence_cases_test::test_analogous_repeat_equivalence`,
  `::test_multi_unit_repeat_equivalence`.

## Protein consequences

Protein consequences are read from codons, not from a diff of two protein strings. The rules that
follow all bear on where the stop is. The unit-test examples below use the coding sequence
`ATG AAA CTG GCC TAT CGC TAA` (Met Lys Leu Ala Tyr Arg Ter) followed by a 3'UTR, and their `c.`
spellings are read from the 0-based ranges the tests assert.

### The declared CDS end is the stop

**The reference protein ends where the transcript record declares the CDS ends, when that codon
really is a stop; only a loosely marked CDS end falls back to the first stop in the translation.**

A selenocysteine `TGA` inside the CDS is therefore not a stop.

**Example** (SEPN1 has selenocysteine codons at 127 and 462; a toy CDS `ATG AAA TGA CTG TAA`
reads Met Lys Sec Leu Ter):

```text
NM_020451.2:c.943G>A   →  NP_065184.2:p.(Gly315Ser)     weaver
                       →  NP_065184.2:p.?                biocommons, seeing two in-frame stops

toy c.10C>G            →  p.Leu4Val                      codon 4 is after the TGA and still residue 4
```

- **Differs:** biocommons.
- **Not checked:** VariantValidator.
- **Spec (silent):** none; the nomenclature assumes the reference protein is known.
- **Tests:** `src/protein.rs::a_selenocysteine_codon_is_not_the_stop`;
  `real_transcripts_test::biocommons_real_transcript_cases` (MULTISTOP01, a recorded difference).

### A frameshift that starts at the stop is an extension

**When the first changed residue is the stop, the protein is extended, not shifted.**

Five ClinVar rows moved from mismatch to Analogous.

**Example** (ATM; and the toy CDS with 3'UTR `CCG TAT AA`):

```text
NM_000051.3:c.9170_9171delGA   →  NP_000042.3:p.(Ter3057PheextTer4)

toy c.20_21del                 →  p.Ter7SerextTer2       the AA of the stop TAA goes: extension
toy c.17_18del                 →  p.Arg6LeufsTer4        two bases earlier, Arg6 changes first: frameshift
```

- **Agrees:** HGVS, biocommons, ClinVar.
- **Spec (agrees):** [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/):
  "the variant extends the amino acid sequence at the C-terminal end and is therefore by
  definition an extension".
- **Tests:** `src/protein.rs::a_frameshift_in_the_stop_codon_is_an_extension`,
  `::stop_loss_is_an_extension_to_the_next_stop`;
  `real_transcripts_test::biocommons_real_transcript_cases` (EXT01 to EXT06).

### An extension needs the stop codon itself to change

**An in-frame change just before the stop that leaves the stop intact is an insertion or
deletion, even though the stop's position is the first residue that differs.**

**Example** (CDS `ATG CTG TTT GTA TTG TGT CGT CTT TAA`, Met Leu Phe Val Leu Cys Arg Leu Ter):

```text
c.21_22insTTGTCT   →  p.Leu8_Ter9insSerLeu    the stop is intact two codons on
toy c.16_18del     →  p.Arg6del               the stop shifts left but is the same stop
```

- **Others:** not checked.
- **Spec (agrees):** [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/),
  the same definition: the stop codon is what changes.
- **Tests:** `src/protein.rs::an_insertion_before_the_stop_that_repeats_the_last_residue_is_not_an_extension`,
  `::an_in_frame_deletion_reaching_the_stop_keeps_it_original`.

### A stop formed inside inserted bases is a delins ending in Ter

**Translation ends before the new frame reads a single reference base, so the change is a delins
ending in Ter, not a frameshift.**

**Example** (the toy CDS; and a two-codon CDS `AAA GGG`):

```text
toy c.4_12delinsCTGTAA     →  p.Lys2_Ala4delinsLeuTer    in frame
toy c.4_12delinsCTGTAAGG   →  p.Lys2_Ala4delinsLeuTer    9 bases out, 8 in: the frame breaks, but after the stop
NM_1.1:c.1_2delinsTA       →  NP_1.1:p.(Lys1Ter)         not p.Lys1fsTer1
```

- **Others:** not checked.
- **Spec (agrees):** [frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/):
  "variants which introduce an immediate translation termination (stop) codon are described as
  nonsense variant", not a frameshift.
- **Tests:** `src/protein.rs::a_stop_inside_the_inserted_bases_is_a_delins_not_a_frameshift`;
  `equivalence_cases_test::test_immediate_stop_normalization`.

### In-frame changes are written 3'-most

**The shared tail is trimmed, so an in-frame deletion or insertion is written at its 3'-most
equivalent residues.**

**Example** (CDS `ATG AAA AAA CTG TAA`, Met Lys Lys Leu Ter):

```text
c.4_6del   →  p.Lys3del    the first Lys codon is deleted; the second is named
c.4_6dup   →  p.Lys3dup
```

- **Agrees:** HGVS.
- **Differs:** ClinVar, in spelling: it writes `Xxx_Yyyins…` at the 3' end of a run where weaver
  writes `dup`, and one-letter repeat forms such as `p.490PRS[1]`; 170 rows in 100,000 are the
  same protein in another spelling, and equivalence judges them Analogous.
- **Spec (agrees):** [protein deletion](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/):
  "the most C-terminal position possible of the reference sequence is arbitrarily assigned to have
  been changed".
- **Tests:** `src/protein.rs::in_frame_deletion_and_duplication_take_the_3_prime_position`;
  `protein_insertions_are_written_3_prime_most` and `in_frame_indels_describe_the_translated_protein`
  (properties).

### A start-codon change is written specifically

**`c.1A>G` is `p.(Met1Val)`, not `p.Met1?`, because the specific prediction carries more
information.**

This is the one choice on this page where weaver's default output departs from an explicit
recommendation, and it is a contentious one. The nomenclature's reasoning is that translation may
initiate elsewhere, so the consequence cannot be predicted; weaver's is that `p.(Met1Val)` says
which codon changed and how, is written as a prediction, and can be reduced to `p.Met1?` but not
recovered from it. `VariantTransformSettings` with `StartCodonConvention::HgvsQuestion` rewrites
it on request.

**Example:**

```text
NM_0001.3:c.1A>T          →  NP_0001.1:p.(Met1Leu)
NM_007199.2:c.1A>G        →  NP_009130.2:p.(Met1Val)     weaver
                          →  NP_009130.2:p.Met1?         biocommons and VariantValidator

transform NP_000051.2:p.(Met1Val)  →  NP_000051.2:p.Met1?    with start_codon = HgvsQuestion
```

- **Differs:** biocommons and VariantValidator write `p.Met1?`, as the nomenclature recommends. A
  recorded difference.
- **Spec (differs):** [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/):
  `p.Met1?` when "the consequence, on the protein level, of a variant affecting the translation
  initiation codon can not be predicted".
- **Tests:** `mapping_test::test_mapper_c_to_p_start_codon_subst`;
  `real_transcripts_test::biocommons_real_transcript_cases` (INITMET01, recorded);
  `src/transform.rs::test_transform_met1_to_question`.

### Edits outside the CDS are statements

**A deletion of the whole CDS is `p.0?`; an edit that starts in the 5'UTR and reaches into the CDS
is `p.Met1?`; an edit entirely upstream is `p.?`; an edit entirely in the 3'UTR is `p.(=)`.**

weaver does not commit to a consequence upstream because the initiation site is not predictable: a
5'UTR change can create an upstream start. Downstream of the stop codon nothing can change the
protein, so the prediction is no change, and no residue past the end of the protein is named.

**Example** (`NM_X.1` has a five-base 5'UTR and a fourteen-base 3'UTR):

```text
NM_X.1:c.-3G>A             →  NP_X.1:p.?         entirely upstream
NM_X.1:c.-3_2del           →  NP_X.1:p.Met1?     reaches into the start codon
NM_X.1:c.-5_*14del         →  NP_X.1:p.0?        the whole CDS
NM_X.1:c.*3A>G             →  NP_X.1:p.(=)       entirely downstream
NM_X.1:c.21A>G             →  NP_X.1:p.(Ter7=)   the stop codon itself, silently

NM_000249.3:c.-7_*46del    →  NP_000240.1:p.0?
NM_022051.2:c.-1_1insGCC   →  NP_071334.1:p.Met1?
```

- **Agrees:** biocommons, on `p.Met1?` and `p.0?` (its own test table).
- **Agrees:** VariantValidator, on `p.(=)` for the 3'UTR (older versions wrote `p.?`); ClinVar
  writes a 3'UTR change as silent at the stop, `p.Ter1648=`, which leaves the same protein.
- **Differs:** ClinVar writes `Met1fs`, `Met1_Glu2insGly…`, committing to a consequence. For a
  substitution entirely in the 5'UTR VariantValidator writes `p.(=)` where weaver writes `p.?`;
  both are defensible.
- **Spec (agrees):** [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/):
  `p.0?` "when you predict that no protein is produced"; `p.Met1?` as above; `p.Met1ext-5` shows
  a 5'UTR change can activate an upstream initiation site.
- **Tests:** `protein_allele_from_coding_test::edits_around_the_cds_start_are_statements_not_predictions`,
  `::a_change_entirely_in_the_three_prime_utr_leaves_the_protein_unchanged`,
  `::statements_have_no_allele_and_a_wrong_protein_is_an_error`,
  `::a_coding_variant_agrees_with_every_spelling_of_its_consequence`;
  `real_transcripts_test::biocommons_real_transcript_cases` (WHOLEGENE01/02, INITMET02/03).

### Frameshift length counts to the first new stop

**`fsTer N` counts residues from the first changed one to the first stop the new frame reaches;
an immediate stop is a plain nonsense substitution.**

**Example** (the toy CDS):

```text
toy c.17_18del   →  p.Arg6LeufsTer4    Leu Ala Val Ter: four residues counting the stop
toy c.4A>T       →  p.Lys2Ter          AAA → TAA, immediate
toy c.4A>C       →  p.Lys2Gln
toy c.6A>G       →  p.Lys2=            AAA → AAG
```

- **Agrees:** HGVS.
- **Spec (agrees):** [frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/):
  the stop position is counted "starting from the first changed amino acid as codon 1"; "the
  shortest frameshift variant possible contains `fsTer2`".
- **Tests:** `src/protein.rs::a_frameshift_in_the_stop_codon_is_an_extension`,
  `::frameshift_reports_first_changed_residue_and_distance_to_stop`,
  `::substitutions_nonsense_and_silence`; `snv_effect_follows_the_codon_table` (property).

### Protein alleles come from the coding change

**From a `p.` description, only edits that name a sequence have an allele. From a coding variant,
every consequence has one, because the edited transcript's translation is known.**

A frameshift's allele is the residues from the first change to the new stop. The translated CDS
must be the protein the provider serves; a difference is an annotation error and is reported as
one. A silent change is the reference allele over the codons the edit touched; a deleted CDS is
the deletion of the whole protein.

**Example** (`NP_X.1` is `MKLAYR`; SPDI is `sequence:position:deleted:inserted`):

```text
from p.:   NP_TEST.1:p.Lys2Ter      →  NP_TEST.1:1:KLAAAYRQ:
           NP_TEST.1:p.Leu3fs       →  UnsupportedOperation: describes a consequence, not a sequence

from c.:   NM_X.1:c.4A>C            →  NP_X.1:1:K:Q             the same allele as its p.Lys2Gln
           NM_X.1:c.5del            →  a frameshift's allele, from Lys2 to the end of the protein
           NM_X.1:c.19T>C           →  NP_X.1:6::QPYK           stop loss, read into the 3'UTR
           NM_X.1:c.-5_*14del       →  NP_X.1:0:MKLAYR:         the whole protein
           NM_X.1:c.-3_2del         →  UnsupportedOperation: p.Met1? names no protein sequence
           with NP_X.1 served as MKLAYQ  →  ValidationError: the CDS does not translate to NP_X.1, differing from residue 6
```

- **Spec (agrees):** VRS 2.0.1 `Allele` on a protein `SequenceReference`; none in HGVS.
- **Tests:** `protein_allele_test::consequences_have_no_allele`;
  `protein_allele_from_coding_test::consequences_without_a_p_sequence_still_have_an_allele`,
  `::definite_changes_give_the_same_allele_by_either_route`,
  `::statements_have_no_allele_and_a_wrong_protein_is_an_error`;
  `the_protein_allele_is_the_edited_translation` (property).

## Validation and stated bases

### Stated bases are checked by validate and nowhere else

**A variant's stated reference, whether a substitution's base or the bases written after `del`,
`dup` or the deleted part of a `delins`, is compared with the sequence by `validate`, and by
`c_to_p`; canonical alleles, SPDI and VRS use the sequence's own bases and ignore what was stated.**

The allele describes the sequence, not the description. This holds for nucleotides and for
protein residues alike.

**Example** (`NM_PLUS10.1` has G at c.1; `NP_TEST.1` is `MKLAAAYRQ`; and on `TTCAGCAGTT`):

```text
validate NM_PLUS10.1:c.1G>A   →  true
validate NM_PLUS10.1:c.1A>G   →  false                     the stated A is not there
validate NM_PLUS10.1:c.1delC  →  false                     so is a base stated after del or dup
validate NM_PLUS10.1:c.1_3del3insAA →  true                a count states no bases
validate NM_MINUS10.1:c.1+5T>A →  true                     intronic: accepted unchecked
validate NP_TEST.1:p.Arg2Leu  →  false                     residue 2 is Lys
validate NP_TEST.1:p.Trp3fs   →  false                     a frameshift's named residue is still checked

allele of g.4G>A on TTCAGCAGTT →  X:3:A:A                   the sequence has A: nothing changes, and the allele says so
```

- **Agrees:** VRS; VariantValidator, which reports a reference mismatch for `delC` on a record
  without the C.
- **Spec (agrees):** VRS 2.0.1: an `Allele` is a `SequenceLocation` and a state; the reference is the
  sequence's.
- **Tests:** `src/allele.rs::a_stated_reference_that_disagrees_with_the_sequence_is_ignored`;
  `transcript_coordinates_test::validate_checks_stated_reference_through_transcript_coordinates`;
  `protein_allele_test::validation_checks_the_named_and_stated_residues`;
  `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
  (the `validate` cases).

## Equivalence

### Judged by allele and by the protein left behind

**Two nucleotide variants are the same change exactly when their canonical alleles are equal; a
coding variant and a protein description agree when they leave the same protein.**

`Identity` is the same text, or a prediction written exactly as given. [How it
decides](equivalence_logic.md) has the full rules.

**Example:**

```text
NP_001337263.1:p.Tyr165Ter   vs  NP_001337263.1:p.Ala164_Tyr165insTer      Analogous: both end the protein after 164
NM_X.1:c.4A>C                vs  NP_X.1:p.(Lys2Gln)                        Identity: the prediction as weaver writes it
NM_X.1:c.4A>C                vs  NP_X.1:p.Lys2Gln                          Analogous
NM_X.1:c.4A>C                vs  NP_X.1:p.K2Q                              Analogous
NP_001.1:p.Arg97ProfsTer4    vs  NP_001.1:p.Arg97delinsProAlaValTer        equivalent: Pro Ala Val Ter is four residues
NP_001.1:p.Arg97ProfsTer4    vs  NP_001.1:p.Arg97delinsProAlaValLeuTer     not: one residue too many
```

- **Differs:** biocommons compares normalised text.
- **Spec (silent):** none; HGVS describes, it does not compare.
- **Tests:** `transcript_coordinates_test::canonical_alleles_make_spdi_vrs_and_equivalence_one_value`;
  `equivalence_cases_test::test_clinvar_regression_tyr165ter`,
  `::test_analogous_fs_wildcard_unification`;
  `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`;
  `a_variant_agrees_with_its_projection` (property).

### A description that says nothing matches nothing

**`p.?` matches no other description; `p.Met1?` says only where the change starts and matches
only a description whose change starts there.**

It does not match ClinVar's `Met1fs`, a commitment weaver did not make. Two ClinVar rows judge
Different for this reason, on purpose.

**Example:**

```text
NM_X.1:c.-3_2del     vs  NP_X.1:p.Leu3fs    not equivalent: the deletion is p.Met1?, which promises nothing about Leu3
NM_X.1:c.-5_*14del   vs  NP_X.1:p.0?        equivalent
NM_X.1:c.-5_*14del   vs  NP_X.1:p.0         equivalent
```

- **Differs:** ClinVar commits.
- **Spec (agrees):** [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/)
  defines `p.?` and `p.Met1?` as statements of what is not known.
- **Tests:** `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`.

### Judging with no protein sequence is an error

**Two protein descriptions are compared on the protein they leave, which needs the sequence;
without it the judgement is an error, not a guess.**

Earlier versions reconciled the residues two descriptions happened to name.

**Example:**

```text
NP_0001.1:p.Ala201_Val202insGlyProGlyAla  vs  NP_0001.1:p.Gly198_Ala201dup
    → Err, not Unknown: residues 199 and 200 are named by neither, and may well make these one change
```

- **Spec (silent):** none.
- **Tests:** `equivalence_test::equivalence_needs_the_sequence`.

### Versions of one protein accession compare on ours

**ClinVar often names an older `NP_` version; both descriptions are read against the protein the
provider serves. A different accession is not the same protein.**

**Example:**

```text
NP_X.1:p.(Lys2Gln)   vs  NP_X.0:p.Lys2Gln   equivalent
NP_X.1:p.Lys2Gln     vs  NP_Y.1:p.Lys2Gln   not equivalent
```

- **Spec (silent):** none.
- **Tests:** `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`.

### A cis allele compares as a set

**`c.[a;b]` and `c.[b;a]` are the same change; a single-member `c.[a]` compares as `a`; a cis
allele never equals a single variant of two or more members.**

**Example** (`NM_X.1` sits at genomic index 10 of `NC_X.1`, so c.7 is g.22):

```text
NM_X.1:c.[7C>T;13T>G]   vs  NM_X.1:c.[13T>G;7C>T]        Analogous
NM_X.1:c.[7C>T;13T>G]   vs  NC_X.1:g.[22C>T;28T>G]       Analogous: the same members on the genome
NM_X.1:c.[7C>T;13T>G]   vs  NM_X.1:c.[7C>T]              Different: a member missing
NM_X.1:c.[7C>T]         vs  NM_X.1:c.7C>T                Identity
NM_X.1:c.[7C>T;13T>G]   vs  NM_X.1:c.7C>T                Different
```

- **Spec (agrees):** [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/): variants
  "should be listed in genomic order", which is a writing convention, not a different allele.
- **Tests:** `cis_phased_test::cis_alleles_compare_as_sets_of_their_members_alleles`.

## Alleles, SPDI and VRS

### The canonical allele is fully justified

**A change is trimmed to what it alters and then widened over the whole region in which it could
equally be written, so every spelling of one change is one allele with one identifier.**

`to_spdi_unambiguous` renders it. The plain `to_spdi` renders the 3'-normalised variant instead,
with an insertion placed at its second flank as SPDI counts.

**Example** (on `TTCAGCAGTT`, SPDI `sequence:position:deleted:inserted`; and on `ACGT` repeated):

```text
X:g.3_5del  and  X:g.6_8del          →  X:2:CAGCAG:CAG           one allele: the run, one unit shorter
X:g.3_5dup  and  X:g.5_6insCAG       →  X:2:CAGCAG:CAGCAGCAG     one allele
X:g.2_3insGG                         →  X:2::GG                  nothing to widen over
X:g.4_5delAGinsAT                    →  X:4:G:T                  trimmed to the base that changes

plain to_spdi NC_TEST.1:g.1012_1013insCC   →  NC_TEST.1:1012::CC   3'-normalised, at the second flank
plain to_spdi NC_TEST.1:g.1013A>G          →  NC_TEST.1:1012:A:G
```

- **Agrees:** VRS and VOCA; this is their normalisation.
- **Spec (agrees):** VRS 2.0.1 "Normalization" (fully justified alleles); refget/SPDI interbase
  coordinates.
- **Tests:** `src/allele.rs::deletion_anywhere_in_a_run_is_the_same_allele`,
  `::duplication_and_insertion_of_the_unit_are_the_same_allele`,
  `::literal_insertion_stays_literal_and_substitution_is_trimmed` (read from the 0-based ranges
  they assert); `transcript_coordinates_test::plain_spdi_places_an_insertion_at_its_second_flank`;
  `ambiguous_range_is_sound_and_complete` and `vrs_and_spdi_round_trip` (properties).

### Refget accessions are computed over the normalised sequence

**Letters uppercased, everything else dropped, as the refget specification defines.**

Without a `Refget` lookup the accession is computed from the whole sequence, and `from_vrs` then
needs the accession passed. The digest is checked against the sequence either way.

**Example:**

```text
ACGT             →  SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2     the specification's own example
acgt             →  the same
AC⏎G T⏎          →  the same
chr19 of the NCBI GRCh38 FASTA  →  SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl, the specification's published value
```

- **Spec (agrees):** [refget](https://samtools.github.io/hts-specs/refget.html), "Sequence Normalization"
  and the `ga4gh` identifier.
- **Tests:** `src/vrs.rs::refget_accession_is_over_the_normalised_sequence`,
  `::digest_and_identifier_match_the_spec_example`;
  `tests/test_digest_table.py::test_digest_follows_the_refget_normalisation_rule`;
  `vrs_round_trip_test::the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked`.

### Uncertain breakpoints become Range bounds

**A deletion with uncertain breakpoints renders as an Allele with Range bounds and an empty
state, unnormalised; a duplication with uncertain breakpoints has no Allele, since its bases are
unknown, and renders as a CopyNumberChange gain.**

An unknown bound prints as `?`; a parenthesised exact position is exact.

**Example:**

```text
NC_TEST.1:g.(?_5)_(10_?)del      →  Allele  "start":[null,4]  "end":[10,null]  "state":{"type":"LiteralSequenceExpression","sequence":""}
NC_TEST.1:g.(3_5)_(10_12)del     →  Allele  start [2,4], end [10,12]
NC_TEST.1:g.(5)_(12)del          →  the same Allele as g.5_12del: exact, so normalised
NC_TEST.1:g.(3_5)_(10_12)dup     →  CopyNumberChange  "copyChange":"gain"  start [2,4], end [10,12]
NC_TEST.1:g.(3_5)_(10_12)inv     →  UnsupportedOperation
```

- **Spec (agrees):** [uncertain positions](https://hgvs-nomenclature.org/stable/recommendations/uncertain/):
  `(A_B)_(C_D)` "where `B_C` describes the minimal extent and `A_D` the maximal", `?` for
  unknown; VRS 2.0.1 `Range` and `CopyNumberChange`.
- **Tests:** `vrs_range_test::uncertain_breakpoints_become_ranges`,
  `::an_imprecise_deletion_has_an_empty_literal_state_and_json_ranges`,
  `::only_deletions_may_have_uncertain_breakpoints`, `::unknown_bounds_round_trip_as_question_marks`;
  `vrs_copy_change_test::an_imprecise_duplication_is_a_gain_over_range_bounds`,
  `::to_vrs_variation_gives_a_change_for_an_imprecise_dup_only`;
  `spec_summary_test::test_uncertain_intervals`.

### copyChange is a label

**VRS 2.0.1 defines `copyChange` as a string enum; weaver writes labels, and reads labels, bare EFO
CURIEs and the 2.0.0 `MappableConcept` form.**

The specification's published `CopyNumberChange` example digest only reproduces when `copyChange`
is digested as the pre-release string `EFO:0030071`; both digests are pinned in tests so the
machinery is anchored to a published value while output follows the current schema. EFO's
high-level loss is `EFO:0020073`, not in the `00300xx` block.

**Example:**

```text
NC_TEST.1:g.(3_5)_(10_12)dup   →  "copyChange":"gain"

read back:
  "gain"  "low-level gain"  "EFO:0030070"
  {"primaryCoding":{"code":"EFO:0030070","system":"https://www.ebi.ac.uk/efo/"}}   →  NC_TEST.1:g.(3_5)_(10_12)dup
  "loss"  "complete genomic loss"  "EFO:0030067"  "EFO:0020073"                     →  NC_TEST.1:g.(3_5)_(10_12)del
  "regional base ploidy"  "amplification"                                          →  UnsupportedOperation
```

- **Agrees:** VRS 2.0.1; 2.0.0 used EFO codes in a `MappableConcept`.
- **Spec (agrees):** VRS 2.0.1 `vrs-source.yaml`, `CopyNumberChange.properties.copyChange.enum`.
- **Tests:** `src/vrs.rs::copy_number_change_digests_the_label_over_the_location`,
  `::copy_change_terms_are_read_by_label_or_efo_code`;
  `vrs_copy_change_test::the_whole_gain_family_is_a_dup_and_the_loss_family_a_del`.

### CisPhasedBlock members are sorted before digesting

**The identifier does not depend on the order the members are written; the `members` array keeps
that order.**

The schema's `minItems: 2` is relaxed to one, a departure from VRS, because HGVS allows the
degenerate `c.[145C>T]` and a parser that accepts it needs somewhere to put it. The in-trans form
`[..];[..]` is rejected: it describes two molecules.

**Example:**

```text
NM_X.1:c.[7C>T;13T>G]   and  NM_X.1:c.[13T>G;7C>T]    →  the same ga4gh:CPB. id, members in the order written
NM_X.1:c.[7C>T;13T>G]   and  NC_X.1:g.[22C>T;28T>G]   →  the same id
NM_X.1:c.[7C>T]                                        →  a block of one
NM_X.1:c.[7C>T];[13T>G] →  UnsupportedOperation: describes two molecules (alleles in trans); parse each [..] as a cis allele instead
```

- **Spec (differs on minItems):** VRS 2.0.1 computed identifiers, "order arrays of digests and ids by Unicode Character
  Set values"; [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) for
  cis `[a;b]` against trans `[a];[b]`.
- **Tests:** `src/vrs.rs::cis_phased_block_digest_matches_the_spec_validation_data`,
  `::a_cis_phased_block_carries_its_members_in_the_order_given`;
  `cis_phased_test::the_identifier_is_the_same_whichever_way_the_members_are_written`,
  `::alleles_in_trans_are_two_molecules_and_are_refused`.

### Reading back gives the normalised variant

**`from_vrs` and `from_spdi` give the 3'-normalised variant on the allele's own sequence.**

An insertion trimmed to before the first base becomes a delins of that base, since HGVS has no
insertion before base 1; a replacement that is the reverse complement of what it replaces reads
back as `inv`; a `LengthExpression` comes back as `insN[n]`; a Range of copies is refused, since
HGVS has no spelling for it.

**Example** (`NC_TEST.1` is `ACGTTTGCAAGGCTAGCTAGCTTTTAACGGGATCGATCGA`):

```text
NC_TEST.1:g.4del        → VRS → NC_TEST.1:g.6del              the 3'-most spelling
NC_TEST.1:g.3_4insT     → VRS → NC_TEST.1:g.6dup
NC_TEST.1:g.13_16del    → VRS → NC_TEST.1:g.19_22del          CTAG rolls over CTAG CT
NC_TEST.1:g.6_7inv      → VRS → NC_TEST.1:g.6_7inv            TG → CA is its own reverse complement
insert T before base 1  → VRS → NC_TEST.1:g.1delinsTA
NC_TEST.1:g.10_11ins(20)        → VRS → NC_TEST.1:g.10_11insN[20]
from_spdi NC_TEST.1:6:1:C       →  NC_TEST.1:g.7G>C
from_spdi NP_TEST.1:1:K:L       →  NP_TEST.1:p.Lys2Leu
"copies":[3,null]               →  UnsupportedOperation
```

- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/) (no position
  before 1); [inversion](https://hgvs-nomenclature.org/stable/recommendations/DNA/inversion/);
  [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/) for `insN[n]`.
- **Tests:** `vrs_round_trip_test::nucleotide_variants_round_trip_to_their_normalised_form`,
  `::alleles_from_other_producers_parse`,
  `::an_insertion_before_the_first_base_reads_back_as_a_delins_of_that_base`;
  `vrs_length_test::length_expressions_read_back_as_the_recommended_spelling`;
  `vrs_copy_number_test::from_vrs_checks_the_count_it_is_given`; `vrs_and_spdi_round_trip`
  (property).

### Breakends and fusions are not rendered

**`Adjacency`, `Terminus` and `DerivativeMolecule` describe breakends and fusions this grammar
has no input for.**

**Example:** none; there is no HGVS input that would produce them.

- **Spec (silent):** VRS 2.0.1.
- **Tests:** none.

## Parsing

### The grammar is checked against the biocommons table

**580 inputs over the 92 rules the two grammars share, agreeing on all but one.**

weaver accepts a terminator inside an amino acid sequence, which biocommons rejects and ClinVar
writes.

**Example:**

```text
rule aat13_seq, input TerGly        →  accepted by weaver, rejected by biocommons
NP_…:p.Glu26_Glu27insTerGlu         →  parses; ClinVar writes it
```

- **Differs:** biocommons, on that one input.
- **Agrees:** ClinVar.
- **Spec (silent):** the nomenclature pages above; the table is biocommons's reading of them.
- **Tests:** `grammar_test::biocommons_grammar_table` (the difference is listed in
  `KNOWN_DIFFERENCES`); `parser_round_trips_canonical_hgvs` and `parser_never_panics`
  (properties).

### Forms accepted beyond biocommons

**Alleles in cis for every coordinate system; predicted RNA changes and RNA statements; protein
statements; uncertain genomic breakpoints; insertions of a stated length; copy number; `r.`
positions with `*`.**

**Example** (each parses and prints back as written):

```text
NM_004006.2:c.[145C>T;147C>G]      NP_X.1:p.[Lys2Leu;Ala4del]      NM_X.1:r.[7c>u;13u>g]
NM_R.1:r.(10c>g)                   NM_R.1:r.0   NM_R.1:r.?   NM_R.1:r.spl   NM_R.1:r.=
NP_X.1:p.0   NP_X.1:p.0?   NP_X.1:p.Met1?
NC_000009.11:g.(?_108337304)_(108337428_?)del
NC_TEST.1:g.10_11insN[20]   NC_TEST.1:g.10_11insN[(20_30)]   NC_TEST.1:g.10_12delinsN[5]   NC_TEST.1:g.10_11ins(20)
NC_000014.8:g.88401076_88459508copy4
NM_R.1:r.*5u>a
```

- **Differs:** biocommons rejects these.
- **Spec (agrees):** [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/);
  [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/);
  [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/);
  [uncertain positions](https://hgvs-nomenclature.org/stable/recommendations/uncertain/);
  [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/).
- **Tests:** `cis_phased_test::cis_alleles_parse_and_print_back_in_every_coordinate_system`;
  `rna_test::predicted_changes_keep_their_parentheses`, `::statements_about_the_transcript_round_trip`;
  `spec_summary_test::test_uncertain_intervals`, `::test_spec_summary_variants`;
  `vrs_length_test::the_older_spelling_gives_the_same_allele`, `::an_uncertain_length_is_a_range`.

### Recommended spellings on output

**`insN[20]` rather than `ins(20)`; `extTer8` rather than `ext*8`; bare `del` and `dup`.**

**Example:**

```text
NC_TEST.1:g.10_11ins(20)       →  NC_TEST.1:g.10_11insN[20]
NC_TEST.1:g.10_12del3insN[5]   →  NC_TEST.1:g.10_12delinsN[5]
NM_000051.3:c.9170_9171delGA   →  NP_000042.3:p.(Ter3057PheextTer4)     biocommons and VariantValidator write ext*4; compared as equal
```

- **Differs:** VariantValidator and biocommons write `ext*`.
- **Spec (agrees):** [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/);
  [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/), where
  "both notations are acceptable";
  [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/).
- **Tests:** `vrs_length_test::length_expressions_read_back_as_the_recommended_spelling`;
  `real_transcripts_test::biocommons_real_transcript_cases` (`ext*` compared equal to `extTer`).

### A range written backwards is refused when parsed

**A range whose start is after its end is a `ParseError`, not a variant that projects to a
plausible-looking genomic range and panics when its bases are sliced.**

Transcript positions order by region (5'UTR and CDS before the 3'UTR), then base, then intronic
offset. A range built in code that runs backwards is a `ValidationError` where it is resolved.
Published variants do contain them, as typos.

**Example:**

```text
NM_025137.4:c.100_50del    →  ParseError: runs backwards       and g.1100_1050del, p.Lys10_Leu5del, c.88-1_87+1del
NM_206933.4:c.8559_2A>G    →  ParseError                       as published
NM_025137.4:c.6331_6232insG →  ParseError                      as published
NM_X.1:c.87+1_88-1del      →  parses                           an intron, in order
NM_X.1:c.-5_10del          →  parses
```

- **Differs:** biocommons parses such a range and refuses it only in its intrinsic validator.
- **Agrees:** VariantValidator refuses it.
- **Spec (agrees):** [general recommendations](https://hgvs-nomenclature.org/stable/recommendations/general/),
  positions in a range are given 5' to 3'.
- **Tests:** `inverted_range_test::a_range_written_backwards_is_refused_when_parsed`,
  `::a_range_built_backwards_is_an_error_wherever_it_is_resolved`.

## The data contract

### A range past the end of a sequence returns the bases that exist

**`get_seq` must return what is there rather than error or return nothing.**

The core pages through sequences in fixed blocks, and a short final block is how it learns where
the sequence ends.

**Example** (a 40-base sequence):

```text
get_seq(ac, 35, 100)   →  the last five bases
get_seq(ac, 100, 110)  →  ""
```

- **Spec (silent):** none; weaver's own interface.
- **Tests:** `reference_store_agrees_with_direct_slicing` (property);
  `tests/test_refget.py::test_sequences_come_from_the_server_with_ranges_clamped`.

### Interval methods are half-open and 0-based

**HGVS positions are 1-based and inclusive, with `c.` skipping the non-existent position 0;
interval methods return 0-based half-open ranges. Exon `reference_end` is inclusive.**

**Example** (`NM_PLUS10.1` and `NM_MINUS10.1` are one exon at genomic 1000..1099 with the CDS at
transcript index 10):

```text
NM_PLUS10.1:c.1A>G      →  (1010, 1011)
NM_PLUS10.1:c.*1A>G     →  (1040, 1041)     the base after the stop
NM_MINUS10.1:c.1A>G     →  (1089, 1090)
NM_MINUS10.1:c.1_3del   →  (1087, 1090)     low to high on the genome
```

- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "there is no
  nucleotide `c.0`".
- **Tests:** `transcript_coordinates_round_trip` and `intron_offsets_agree_from_either_exon`
  (properties); `transcript_coordinates_test::spdi_interval_honours_cds_end_anchor`,
  `::spdi_interval_honours_non_zero_cds_start`, `::spdi_interval_on_minus_strand`.

### A mapper keeps its cache

**A `VariantMapper` holds sequence blocks and refget accessions for as long as it lives. Build one
and reuse it.**

**Example:**

```python
mapper = weaver.VariantMapper(provider)
mapper.c_to_p(weaver.parse("NM_TEST.1:c.4G>A"))   # fetches
mapper.c_to_p(weaver.parse("NM_TEST.1:c.4G>A"))   # no new fetch
weaver.VariantMapper(provider).c_to_p(...)         # a new mapper starts cold
```

- **Spec (silent):** none.
- **Tests:** `tests/test_mapper_cache.py::test_sequences_are_fetched_once_per_mapper`,
  `::test_a_refget_lookup_is_asked_before_the_sequence_is_hashed`.

### Refget is its own seam

**A sequence source need know nothing about digests; the `Refget` lookup is a separate object,
and a provider without one still works for everything but naming the sequence behind a digest.**

**Example:**

```python
table = DigestTable.from_fasta("genome.fa")            # or a RefgetProvider, or any object with the two methods
mapper = weaver.VariantMapper(provider, refget=table)
allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.7G>C"))
str(mapper.from_vrs(allele))                           # "NC_TEST.1:g.7G>C", the accession looked up from the digest

weaver.VariantMapper(provider).from_vrs(allele)                        # DataProviderError: no lookup
weaver.VariantMapper(provider).from_vrs(allele, "NC_TEST.1")           # works: the caller names the sequence
```

- **Spec (agrees):** [refget](https://samtools.github.io/hts-specs/refget.html).
- **Tests:** `vrs_round_trip_test::the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked`;
  `tests/test_refget.py::test_refget_lookups_both_ways`;
  `tests/test_digest_table.py::test_the_table_agrees_with_what_the_mapper_computes_and_names_the_sequence_back`.
