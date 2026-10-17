# The Choices Weaver Makes

HGVS leaves room for judgement, and the tools that read it have settled on different answers. This page lists every
deliberate choice weaver makes where another tool could reasonably answer otherwise: what weaver does, why, a worked
example, the recommendation it rests on, and what other tools do.

## The principles that decide a choice

[`PRODUCT.md`](../PRODUCT.md) sets weaver's direction: positions are types, consequences are read from the sequence,
every choice is written down, and the answer is an answer or an error, never a guess. The principles below sit one level
beneath those. They are how an entry on this page is decided, stated once so that a reader can see why the entries say
what they do, and so that the next choice can be checked against them. Each closes with the entries that apply it.

### The specification first, then the plain reading of the notation

**Where the HGVS nomenclature, VRS or refget says what a description means, weaver does that; where it is silent, weaver
takes the mechanical reading of what the notation says; and every departure, and every silence weaver fills, is recorded
here with its reason.**

A tool that fills a silence by taste gives answers nobody can predict. Reading the notation mechanically gives the
answer a careful reader would derive from the text, and recording each such reading lets someone disagree with one line
rather than with a behaviour. The **Spec** verdict on every entry is the ledger: *agrees*, *silent*, or *differs* and
why.

The specification writes a `?` after a residue only at the start codon; weaver reads `p.Arg12?` as it reads `p.Met1?`, a
change from that residue on of which nothing more is said, and records that reading as its own.

Applied in: [a description that says nothing matches nothing](#a-description-that-says-nothing-matches-nothing),
[a start-codon change is written specifically](#a-start-codon-change-is-written-specifically) (a departure),
[CisPhasedBlock members are sorted before digesting](#cisphasedblock-members-are-sorted-before-digesting) (a departure),
[a transcript base the genome lacks projects as an insertion](#a-transcript-base-the-genome-lacks-projects-as-an-insertion),
[mitochondrial and RNA descriptions share the DNA machinery](#mitochondrial-and-rna-descriptions-share-the-dna-machinery),
[the whole transcript unchanged is a statement](#the-whole-transcript-unchanged-is-a-statement-not-a-change).

### One reading is an answer; none, or several, is an error

**A description with one mechanical reading gets an answer, however unusual it looks. One with no reading is
contradictory, and is an error saying why. One with several readings that make different sequences, or that state
different bases, is ambiguous, and is an error that lists them: weaver never picks.**

`PRODUCT.md` says weaver answers or refuses; this says where the line falls. Refusing a description that has one reading
throws away a right answer, and answering one that has several is a guess. Keeping the two errors apart lets a caller
who knows which reading was meant take it.

`c.[=;10C>G]` says the whole transcript was analysed and no change found, and names a change: it has no reading and is
refused when parsed. `c.[10=;12C>G]`, a position unchanged beside a change, has one, and parses.

Applied in:
[the whole transcript unchanged is a statement](#the-whole-transcript-unchanged-is-a-statement-not-a-change),
[a position outside every exon is an error](#a-position-outside-every-exon-is-an-error),
[an open CDS end numbers no position](#an-open-cds-end-numbers-no-position),
[a range written backwards is refused when parsed](#a-range-written-backwards-is-refused-when-parsed),
[a stated base the reference lacks is an error](#a-stated-base-the-reference-lacks-is-an-error),
[judging with no protein sequence is an error](#judging-with-no-protein-sequence-is-an-error),
[overlapping members of a cis allele need one reading](#overlapping-members-of-a-cis-allele-need-one-reading),
[an ambiguous description is an error that lists its readings](#an-ambiguous-description-is-an-error-that-lists-its-readings).

### A description claims what it says, and no more

**Where the text of a description leaves something open, weaver's answer leaves it open too, rather than commit to the
likeliest value.**

`p.?` and `p.Met1?` record what is not known. Read as a specific change, they would match descriptions they never
claimed; and a weak reading can be strengthened by someone who knows more, where a strong one cannot be withdrawn. The
same holds for what weaver derives: it does not invent a genomic form, an allele or a protein for a statement that names
none.

`NM_X.1:c.-3_2del` loses the start codon and predicts `p.(Met1?)`; it is not equivalent to `p.Leu3fs`, because it says
nothing about Leu3.

Applied in: [a description that says nothing matches nothing](#a-description-that-says-nothing-matches-nothing),
[edits outside the CDS are statements](#edits-outside-the-cds-are-statements),
[the whole transcript unchanged is a statement](#the-whole-transcript-unchanged-is-a-statement-not-a-change),
[protein alleles come from the coding change](#protein-alleles-come-from-the-coding-change),
[uncertain breakpoints become Range bounds](#uncertain-breakpoints-become-range-bounds).

### Equality is the allele

**Two descriptions are the same change when they make the same sequence and state the same bases unchanged. Their
spelling, their coordinate system, and anything they record that is not part of the allele, such as a stated reference
base, do not enter. A description states bases unchanged through its identity members, `c.10=`, and through members that
cancel each other; one that only changes the sequence states none. Nucleotide descriptions are compared on the genome,
because a change written on a transcript is a change to the genome.**

Where an identity or a set of cancelling members is written is what it states, not a spelling: on a run of A,
`c.[4_5insA;5del]` states c.4 and c.5, and `c.[5_6insA;6del]` states c.5 and c.6. One change has many spellings: shifted
or not, `dup` or `ins`, `c.` or `g.`, with its bases stated or bare. Text comparison sees the spelling. Comparing on the
genome rather than the transcript keeps equivalence transitive: two transcript changes that leave one transcript can
still be different changes to the genome when they sit on either side of an exon junction.

On `TTCAGCAGTT`, `g.3_5del` and `g.6_8del` are one canonical allele, `X:2:CAGCAG:CAG`. `c.10=` and `c.20=` both make the
reference, and are different statements; `c.[7C>T;10_12=]` makes what `c.7C>T` makes, and says more.

Applied in: [judged by allele and by the protein left behind](#judged-by-allele-and-by-the-protein-left-behind),
[the canonical allele is fully justified](#the-canonical-allele-is-fully-justified),
[a repeat resolves to its whole run](#a-repeat-resolves-to-its-whole-run),
[a transcript base the genome lacks projects as an insertion](#a-transcript-base-the-genome-lacks-projects-as-an-insertion),
[a cis allele is the sequence its members make](#a-cis-allele-is-the-sequence-its-members-make).

### Only normalize rewrites a description

**Parsing, validating, projecting and comparing keep a change where and as it was written. `normalize` alone rewrites
it, into the preferred spelling for interchange; it keeps the sequence the description makes, normalising twice is
normalising once, and it never changes what the description means on the genome, so nothing it moves crosses an exon
junction.**

Where a change is written can carry what its author meant, and a result moved somewhere else is hard to relate to the
question asked. The preferred spelling serves the readers that need one: tools that match strings, and tools that read a
description's parts one at a time. A normalisation that changed the sequence, or the change on the genome, would hand
them a different variant. What is kept is the change's placement and shape; spellings the specification treats as one
value, `ins(20)` and `insN[20]`, `M` and `Met`, are one value when parsed.

`NM_001372044.2:c.1568dup` projects at c.1568, where VariantValidator shifts it to c.1569 first. When it does rewrite,
`normalize` writes `NM_J.1:c.4del` as `c.5del`, on the last base of its exon, not as `c.6del`, which the transcript
cannot tell from it but which deletes a base of the next exon on the genome.

Applied in: [no normalisation before projecting](#no-normalisation-before-projecting),
[an allele in cis projects member by member](#an-allele-in-cis-projects-member-by-member),
[intronic positions are carried as given](#intronic-positions-are-carried-as-given),
[shift 3', cyclically over repeats](#shift-3-cyclically-over-repeats),
[a shift stops at an exon/exon junction](#a-shift-stops-at-an-exonexon-junction),
[deletions and duplications are written bare](#deletions-and-duplications-are-written-bare),
[a prediction is bracketed as HGVS writes it](#a-prediction-is-bracketed-as-hgvs-writes-it),
[reading back gives the normalised variant](#reading-back-gives-the-normalised-variant),
[members of a cis allele that meet are written as one change](#members-of-a-cis-allele-that-meet-are-written-as-one-change).

### What weaver writes reads back as itself, and as the specification reads it

**A description weaver writes parses back to the value it was written from, and a reader following the specification
takes it to mean what weaver means.**

What weaver writes is stored and read later, by weaver and by other tools. A description that parses to another value
changes meaning in storage, and one the specification reads differently means different things to different readers. So
each form is written as the specification spells it, the parentheses of each statement as the specification gives that
statement, and no form is written that the specification gives another meaning, such as `X`, which has been used for a
stop. The property tests hold every projection, normalisation and reading back from SPDI and VRS to the first half.

`p.M1?` is stored and printed as `p.Met1?`, and the `p.(Met1?)` that `c_to_p` writes for a lost start codon parses to
the value it built.

Applied in: [recommended spellings on output](#recommended-spellings-on-output),
[a prediction is bracketed as HGVS writes it](#a-prediction-is-bracketed-as-hgvs-writes-it),
[a start-codon change is written specifically](#a-start-codon-change-is-written-specifically),
[positions number from 1, and an unknown one is not 0](#positions-number-from-1-and-an-unknown-one-is-not-0).

### Bounds before bases

**Where positions alone decide a question, weaver decides it from the positions, without building or reading sequence.**

A decision taken from bounds costs the same whatever the distance, and needs no sequence that may be unavailable or, for
an unknown bound, unlimited. One that builds sequence to answer a positional question costs in proportion to the
distance, and has nothing to read where a bound is unknown.

A range written backwards is refused when it is parsed, from its positions, rather than projected and found wrong when
its bases are sliced.

Applied in: [a range written backwards is refused when parsed](#a-range-written-backwards-is-refused-when-parsed),
[positions number from 1, and an unknown one is not 0](#positions-number-from-1-and-an-unknown-one-is-not-0),
[uncertain breakpoints become Range bounds](#uncertain-breakpoints-become-range-bounds).

## How to read this page

Every choice has the same shape.

| Part                     | What it holds                                                                                                                                                                                                                                                                                                                               |
| :----------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Example**              | Input on the left, weaver's output on the right, taken from the test that asserts it. `→` reads "becomes".                                                                                                                                                                                                                                  |
| **Agrees** / **Differs** | Only where the other tool was actually checked: by running it (VariantValidator's REST API, the biocommons `hgvs` grammar and regression tables, ClinVar's 100,000-variant sample) or by reading its code. Anything else is marked *not checked*.                                                                                           |
| **Spec**                 | The [HGVS nomenclature](https://hgvs-nomenclature.org/stable/) page (or the VRS or refget specification) the choice rests on, with a verdict: *agrees* when the specification says or clearly implies the same, *differs* when weaver departs from what it recommends, *silent* when it says nothing. Departures are few and each says why. |
| **Tests**                | The tests that assert it, as `file::function`. Files are under `hgvs-weaver/tests/` unless written as `src/…` (unit tests) or `tests/…py` (Python); property tests live in `hgvs-weaver/tests/properties/main.rs`.                                                                                                                          |

The tools referred to: **biocommons** is the Python `hgvs` package weaver was originally a port of; **VariantValidator**
is the web service; **ClinVar** is the descriptions in its variant summary; **VRS** is the GA4GH Variation
Representation Specification 2.0.1.

## At a glance

| Choice                                                                                                                              |                 Spec                 |    biocommons    |   VariantValidator   |     ClinVar      |
| :---------------------------------------------------------------------------------------------------------------------------------- | :----------------------------------: | :--------------: | :------------------: | :--------------: |
| [A projection states the bases of the target](#a-projection-states-the-bases-of-the-target)                                         |                agrees                |      agrees      |        agrees        |        –         |
| [A change across a splice junction has no genomic form](#a-change-across-a-splice-junction-has-no-genomic-form)                     |                agrees                |        –         |       differs        |        –         |
| [No normalisation before projecting](#no-normalisation-before-projecting)                                                           |                  –                   |        –         |       differs        |        –         |
| [Gapped exons are projected locally](#gapped-exons-are-projected-locally)                                                           |                  –                   |        –         |       differs        |        –         |
| [A transcript base the genome lacks projects as an insertion](#a-transcript-base-the-genome-lacks-projects-as-an-insertion)         |                  –                   |     differs      | differs inside a run |        –         |
| [A repeat is projected as its whole run](#a-repeat-is-projected-as-its-whole-run)                                                   |                agrees                |        –         |          –           |        –         |
| [A position outside every exon is an error](#a-position-outside-every-exon-is-an-error)                                             |                agrees                |     differs      |          –           |        –         |
| [An open CDS end numbers no position](#an-open-cds-end-numbers-no-position)                                                         |                agrees                |        –         |          –           |        –         |
| [Intronic positions are carried as given](#intronic-positions-are-carried-as-given)                                                 |                agrees                |        –         |          –           |        –         |
| [Mitochondrial and RNA descriptions share the DNA machinery](#mitochondrial-and-rna-descriptions-share-the-dna-machinery)           |      agrees; silent on `r.(?)`       |        –         |       differs        |        –         |
| [The whole transcript unchanged is a statement, not a change](#the-whole-transcript-unchanged-is-a-statement-not-a-change)          | agrees; silent on what `c.=` implies |     differs      |          –           |        –         |
| [An allele in cis projects member by member](#an-allele-in-cis-projects-member-by-member)                                           |  agrees; silent on the duplication   |        –         |          –           |        –         |
| [Shift 3', cyclically over repeats](#shift-3-cyclically-over-repeats)                                                               |                agrees                |      agrees      |        agrees        |      agrees      |
| [A shift stops at an exon/exon junction](#a-shift-stops-at-an-exonexon-junction)                                                    |                agrees                |     differs      |          –           |      agrees      |
| [A delins is reduced before it is shifted](#a-delins-is-reduced-before-it-is-shifted)                                               |                agrees                |   agrees on g.   |          –           |        –         |
| [Deletions and duplications are written bare](#deletions-and-duplications-are-written-bare)                                         |                agrees                |     differs      |        agrees        |        –         |
| [A repeat resolves to its whole run](#a-repeat-resolves-to-its-whole-run)                                                           |                agrees                |        –         |          –           |      agrees      |
| [Members of a cis allele that meet are written as one change](#members-of-a-cis-allele-that-meet-are-written-as-one-change)         |             **differs**              |        –         |          –           |        –         |
| [The declared CDS end is the stop](#the-declared-cds-end-is-the-stop)                                                               |                  –                   |     differs      |          –           |        –         |
| [A frameshift that starts at the stop is an extension](#a-frameshift-that-starts-at-the-stop-is-an-extension)                       |                agrees                |      agrees      |          –           |      agrees      |
| [An extension needs the stop codon itself to change](#an-extension-needs-the-stop-codon-itself-to-change)                           |                agrees                |        –         |          –           |        –         |
| [A stop formed inside inserted bases is a delins ending in Ter](#a-stop-formed-inside-inserted-bases-is-a-delins-ending-in-ter)     |                agrees                |        –         |          –           |        –         |
| [In-frame changes are written 3'-most](#in-frame-changes-are-written-3-most)                                                        |                agrees                |        –         |          –           | spelling differs |
| [A start-codon change is written specifically](#a-start-codon-change-is-written-specifically)                                       |             **differs**              |     differs      |       differs        |        –         |
| [Edits outside the CDS are statements](#edits-outside-the-cds-are-statements)                                                       |   agrees; silent on no codon left    |      agrees      |       differs        |     differs      |
| [A prediction is bracketed as HGVS writes it](#a-prediction-is-bracketed-as-hgvs-writes-it)                                         |                agrees                |     differs      |          –           |        –         |
| [Frameshift length counts to the first new stop](#frameshift-length-counts-to-the-first-new-stop)                                   |                agrees                |        –         |          –           |        –         |
| [Protein alleles come from the coding change](#protein-alleles-come-from-the-coding-change)                                         |     agrees; silent on derivation     |        –         |          –           |        –         |
| [A stated base the reference lacks is an error](#a-stated-base-the-reference-lacks-is-an-error)                                     |                agrees                | differs intronic |   differs intronic   |        –         |
| [Judged by allele and by the protein left behind](#judged-by-allele-and-by-the-protein-left-behind)                                 |                  –                   |     differs      |          –           |        –         |
| [A description that says nothing matches nothing](#a-description-that-says-nothing-matches-nothing)                                 |     agrees; silent on `p.Arg12?`     |        –         |          –           |     differs      |
| [Judging with no protein sequence is an error](#judging-with-no-protein-sequence-is-an-error)                                       |                  –                   |        –         |          –           |        –         |
| [Versions of one protein accession compare on ours](#versions-of-one-protein-accession-compare-on-ours)                             |                  –                   |        –         |          –           |        –         |
| [A cis allele is the sequence its members make](#a-cis-allele-is-the-sequence-its-members-make)                                     |                  –                   |        –         |          –           |        –         |
| [Overlapping members of a cis allele need one reading](#overlapping-members-of-a-cis-allele-need-one-reading)                       |             **differs**              |        –         |          –           |        –         |
| [An ambiguous description is an error that lists its readings](#an-ambiguous-description-is-an-error-that-lists-its-readings)       |                  –                   |        –         |          –           |        –         |
| [The canonical allele is fully justified](#the-canonical-allele-is-fully-justified)                                                 |                agrees                |        –         |          –           |        –         |
| [Refget accessions are computed over the normalised sequence](#refget-accessions-are-computed-over-the-normalised-sequence)         |                agrees                |        –         |          –           |        –         |
| [Uncertain breakpoints become Range bounds](#uncertain-breakpoints-become-range-bounds)                                             |                agrees                |        –         |          –           |        –         |
| [copyChange is a label](#copychange-is-a-label)                                                                                     |                agrees                |        –         |          –           |        –         |
| [CisPhasedBlock members are sorted before digesting](#cisphasedblock-members-are-sorted-before-digesting)                           |             **differs**              |        –         |          –           |        –         |
| [Reading back gives the normalised variant](#reading-back-gives-the-normalised-variant)                                             |                agrees                |        –         |          –           |        –         |
| [Breakends and fusions are not rendered](#breakends-and-fusions-are-not-rendered)                                                   |                  –                   |        –         |          –           |        –         |
| [The grammar is checked against the biocommons table](#the-grammar-is-checked-against-the-biocommons-table)                         |                  –                   | differs on three |          –           |      agrees      |
| [Forms accepted beyond biocommons](#forms-accepted-beyond-biocommons)                                                               |                agrees                |     differs      |          –           |        –         |
| [Recommended spellings on output](#recommended-spellings-on-output)                                                                 |                agrees                |        –         |       differs        |        –         |
| [A range written backwards is refused when parsed](#a-range-written-backwards-is-refused-when-parsed)                               |                agrees                |     differs      |        agrees        |        –         |
| [Positions number from 1, and an unknown one is not 0](#positions-number-from-1-and-an-unknown-one-is-not-0)                        |                agrees                |        –         |          –           |        –         |
| [A number too large to hold is refused when parsed](#a-number-too-large-to-hold-is-refused-when-parsed)                             |                  –                   |        –         |          –           |        –         |
| [A range past the end of a sequence returns the bases that exist](#a-range-past-the-end-of-a-sequence-returns-the-bases-that-exist) |                  –                   |        –         |          –           |        –         |
| [Interval methods are half-open and 0-based](#interval-methods-are-half-open-and-0-based)                                           |     agrees; silent on intervals      |        –         |          –           |        –         |
| [A mapper keeps its cache](#a-mapper-keeps-its-cache)                                                                               |                  –                   |        –         |          –           |        –         |
| [Refget is its own seam](#refget-is-its-own-seam)                                                                                   |      agrees; silent on the seam      |        –         |          –           |        –         |

**Spec** is the HGVS nomenclature, or VRS or refget where the choice is about them; a dash means it is silent, and a
verdict silent on part of the choice names that part, given in full in the entry. For a tool, a dash means it was not
checked or the choice does not apply to it. The four departures from a specification are in bold.

## Projection between sequences

### A projection states the bases of the target

**Where record and genome differ, the projected edit is written against the target's own bases, as the minimal edit that
produces the same alternate.**

A transcript record and its genome can differ at a base: RefSeq transcripts are curated against submitted mRNAs. When a
variant is projected, the edit is resolved on the source to the bases it removes and the bases it puts there. Where the
target holds different bases over the projected range, weaver writes the edit that turns the target's bases into that
alternate, collapsing to `=` when they already coincide. A duplication or inversion carries the source's bases, because
that is what the molecule becomes. A deletion is positional. Where record and genome agree, which is nearly everywhere,
the edit is carried across unchanged.

**Example** (MUC2 has a record C over a genomic G; SHANK3 a record C over a genomic T):

```text
NM_002457.5:c.12468C>A       →  NC_000011.10:g.1099802G>A        the genome's G, not the record's C
NM_002457.5:c.12468dup       →  NC_000011.10:g.1099802delinsCC   a dup carries the record's C, twice
NM_002457.5:c.12468del       →  NC_000011.10:g.1099802del        a deletion is positional

NM_001372044.2:c.1568=       →  NC_000022.11:g.50697558T>C       no change on the record is a change on the genome
NC_000022.11:g.50697558T>C   →  NM_001372044.2:c.1568=           and the other way round
NM_001372044.2:c.1568C>T     →  NC_000022.11:g.50697558=
```

- **Agrees:** biocommons (`replace_reference`); VariantValidator on substitution, identity, duplication and deletion at
  the MUC2 base.
- **Spec (agrees):** [general](https://hgvs-nomenclature.org/stable/recommendations/general/): "descriptions on
  RNA/protein level should describe the changes observed on that level"; a description is of the sequence it is on.
- **Tests:** `projection_reference_test::a_projection_states_the_genomes_bases_not_the_records`,
  `::a_projection_states_the_transcripts_bases_not_the_genomes`,
  `::the_re_read_is_in_the_targets_orientation_on_the_minus_strand`;
  `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases` (MUC2 and
  SHANK3 with NCBI's alignments); `a_variant_agrees_with_its_projection` (property, where record and genome agree).

### A change across a splice junction has no genomic form

**An `r.` change spanning an exon boundary describes the spliced RNA, and weaver refuses to project it to the genome.**

The `c.` spelling of the same change still projects, to a deletion that includes the intron, and its protein consequence
is still predicted. weaver refuses the `r.` projection because a warning is easy to miss and the two molecules are not
the same.

**Example** (exon 1 of `NM_R.1` ends at r.45):

```text
NM_R.1:r.44_47del   → genome    UnsupportedOperation: "spans a splice junction: the spliced RNA has no single genomic equivalent"
NM_R.1:r.44_47del   → protein   the same frameshift as NM_R.1:c.44_47del
NM_R.1:c.44_47del   → genome    projects, over the intron
NM_R.1:r.45del      → genome    projects: one exon
```

- **Differs:** VariantValidator projects the `c.` reading with the warning "spans at least one intron".
- **Spec (agrees):** [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/) and
  [RNA splicing](https://hgvs-nomenclature.org/stable/recommendations/RNA/splicing/), where exon skipping is written as
  an `r.` deletion across the junction (`r.(3277_3432del)`), an RNA event.
- **Tests:** `rna_test::r_projects_to_the_genome_within_one_exon_only`, `::r_predicts_the_protein_like_c`;
  `rna_across_a_junction_has_a_protein_but_no_genomic_form` (property).

### No normalisation before projecting

**weaver projects the variant where it is written; shifting it is the caller's decision.**

Normalisation is available separately, as `normalize_variant`.

**Example** (SHANK3, record C over a genomic T):

```text
NM_001372044.2:c.1568dup  →  NC_000022.11:g.50697558delinsCC      weaver, at the position given
                          →  c.1569dup first, then the genome      VariantValidator, reporting the automap
```

- **Differs:** VariantValidator 3'-shifts on the transcript first and reports the automapping.
- **Spec (silent):** none; the [3' rule](https://hgvs-nomenclature.org/stable/recommendations/general/) says how to
  write a variant, not that a tool must rewrite what it is given.
- **Tests:** `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
  (positions come back as given).

### Gapped exons are projected locally

**Where an exon aligns to the genome with insertions and deletions, weaver projects through the exon's CIGAR and writes
a local edit.**

**Example** (SHANK3, whose exon alignment is heavily gapped):

```text
NM_001372044.2:c.1568dup  →  NC_000022.11:g.50697558delinsCC                     weaver
                          →  NC_000022.11:g.50695049_50697558delins… (2.5 kb)    VariantValidator
```

- **Differs:** VariantValidator falls back to a whole-region delins on heavily gapped exons.
- **Spec (silent):** none.
- **Tests:** `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases`
  (the SHANK3 dup case); `src/transcript_mapper.rs::test_g_to_n_cigar`.

### A transcript base the genome lacks projects as an insertion

**A transcript base with no genome counterpart (an `I` in its exon's CIGAR) has no genomic position of its own. The
genome can describe only the whole run of such bases, so an edit touching a run is written over all of it: the genome
gets the run, as changed, inserted between its flanking bases. The same holds in the other direction: a genome base the
transcript lacks (a `D`) has no transcript index, and a genomic edit touching the run is written on the transcript over
the whole run. At the edge of an exon's alignment, where a soft-clipped end is supplied as a leading or trailing `I`, or
a `D` opens or closes the alignment, there is no flanking pair and projecting it is a `ValidationError`.**

Such a base has no genomic substitution equivalent. Placing it on a neighbouring genome base, as the alignment-mapper
weaver descends from does, silently replaces the stated base with the neighbour's. Inserting only the changed base loses
where in the run the change was and, read on any other transcript, is a frameshift the record never had. Carrying the
whole run makes the projected genomic variant, applied to the genome, reproduce the record's exon with the change: the
two descriptions are one allele.

**Example** (`TX_GAP3.1` is one exon on `NC_000099.1:g.11_40` whose CIGAR is `15=3I15=`: the record reads `C TCG A` over
n.15_19 where the genome reads `CA` over g.25_26; `TX_CLIP.1` has its first five bases clipped, CIGAR `5I20=`;
`TX_DEL.1` has CIGAR `15=1D14=`, lacking the genome's A at g.26 between its c.15 and c.16):

```text
TX_GAP3.1:n.17C>A      →  NC_000099.1:g.25_26insTAG    the run, as changed, where the run is
TX_GAP3.1:n.17del      →  NC_000099.1:g.25_26insTG
TX_GAP3.1:n.16_18del   →  NC_000099.1:g.25=            the genome already lacks the run
TX_GAP3.1:n.15_17del   →  NC_000099.1:g.25C>G          a flank and part of the run
TX_GAP3.1:n.15_16insA  →  NC_000099.1:g.25_26insATCG   an insertion beside the run carries it
TX_GAP3.1:n.15C>G      →  NC_000099.1:g.25C>G          the flanking bases do not touch the run
TX_CLIP.1:n.1G>C       →  ValidationError              no flanking pair; nothing places it
NC_000099.1:g.26A>G    →  TX_DEL.1:c.15_16insG         the other direction: the run, as changed, on the record
```

- **Differs:** biocommons hgvs places the base on the flanking genome base (the one before as a start, the one after as
  an end) and flags the interval uncertain; its projected substitution then names the neighbour's base.
- **Differs (inside a run):** VariantValidator inserts only the changed base. IGFBP2 `NM_000597.2` has nine bases,
  c.45_53, that GRCh38 lacks; `c.47C>T` and `c.51G>A` both become `NC_000002.12:g.216633567_216633568insT` and `…insA`,
  and `c.49del` becomes `…insTG`, two of the eight bases that remain. HTT `NM_002111.8:c.54G>T`, in its six-base run
  c.52_57, is likewise `NC_000004.12:g.3074876_3074877insT`.
- **Agrees (beside a run, and the whole run):** VariantValidator writes `NM_000597.2:c.45_53del` as
  `g.216633564_216633569=` and carries the run for an edit beside it: `c.53_54insA` is
  `g.216633576_216633577insAGCCGCTGCT`. It also carries the run for the flanking bases themselves (`c.44C>T` is
  `g.216633567delinsTGCCGCTGCT`), which weaver does not: a flanking base is aligned and does not touch the run.
- **Spec (silent):** none; the nomenclature does not address transcript–genome alignment gaps.
- **Tests:** `gapped_exon_test` (every row above, in both directions, through n. and c., on both strands);
  `src/transcript_mapper.rs::a_base_the_genome_lacks_occupies_an_empty_range_between_its_neighbours`,
  `::test_g_to_n_cigar`.

### A repeat is projected as its whole run

**A repeat is widened to the full run of its unit on the source before projecting, and the target counts its own run:
copies the target reads beyond the source's run are added to the count, and a transcript run that crosses an exon
junction is refused in either direction.**

A repeat's count is the total number of units in the run read 3' from its first base, on the sequence it is written on.
That run is not always the same on the other sequence. On the minus strand the genome reads the run from the
transcript's 3' end towards its 5' end, so where the bases before the written start continue the run, the genome counts
them too; on either strand an intron can continue a run that ends at an exon's edge. Carrying `[3]` across unchanged
would then describe a different change. So weaver widens the repeat to its run on the source (stating the unit, so an
unstated one keeps the bases its range named), and on the target adds the copies its own run holds beyond that range,
which the change leaves alone. At an intronic position the transcript has no bases of its own, so on the minus strand
weaver reads the transcript's run on the genome in the transcript's sense, downwards, and then the genome's upwards in
the same way. So whether a run continues past an exon's edge depends on how the repeat is written: one whose ends are
both exonic counts the record's run, one with an intronic or flanking end counts the genome's, where its bases are.
Where the target's bases over the run are not the source's, the edit is the minimal one, as
[for any projection](#a-projection-states-the-bases-of-the-target), and a count known only as a range is refused because
it has no single alternate there. A transcript run that crosses an exon junction is refused, whichever side the repeat
is written on: copies added to or removed from it could go on either side of the intron, so a `c.` count over it has no
single place on the genome, and a genomic repeat whose run the next exon continues would be written as that `c.` count.
Counting only the copies within the starting exon would answer instead, at the price of a `c.` count that is not the
transcript's total, which the specification defines it to be.

**Example** (`NC_REP.1` is `TTTTT GCCATT GCCATT GCCATT AAAAA…`; `NM_REP_MINUS.1` reads it on the minus strand, where the
run is `AATGGC` from c.78. `NM_REV.1` is `ATGAAACAGCAGCAGTTCTAAGG` on the minus strand of `NC_000099.1`, c.4_6 `AAA` at
g.28_30; `NM_EDGE.1`'s exon 1 ends `GTGT` and its intron opens `GT`; `NM_INT.1` is on the minus strand with c.6+4_7-5
`AAAA` at g.23_26):

```text
NM_REP_MINUS.1:c.78AATGGC[4]  →  NC_REP.1:g.6_23GCCATT[4]           the whole run, in plus-strand orientation
NC_REP.1:g.6GCCATT[4]         →  NM_REP_MINUS.1:c.78_95AATGGC[4]
NM_REV.1:c.5A[3]              →  NC_000099.1:g.28_30T[4]            two A become three; the genome's run has c.4 too
NC_000099.1:g.29T[3]          →  NM_REV.1:c.4_6A[4]                 and the other way round
NM_REV.1:c.5A[(3_5)]          →  NC_000099.1:g.28_30T[(4_6)]
NM_EDGE.1:c.6GT[3]            →  NC_000097.1:g.16_21GT[4]           the intron continues the run
NM_INT.1:c.6+5A[4]            →  NC_000094.1:g.23_26T[5]            three A become four, read in the intron
NC_000094.1:g.25T[4]          →  NM_INT.1:c.6+4_7-5A[6]
NM_EDGE.1:c.7TG[3]            →  UnsupportedOperation               the run crosses the junction after c.9
NC_000095.1:g.15A[3]          →  UnsupportedOperation               NM_JX.1's run of A carries on into exon 2
```

- **Others:** not checked.
- **Spec (agrees):** [repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/): the count
  is the total number of units in the run. Silent on a run that crosses an exon junction.
- **Tests:** `transcript_coordinates_test::repeat_on_the_minus_strand_projects_to_its_whole_run`;
  `repeat_projection_test` (every row above, and the minimal edit over a base the genome does not share);
  `a_coding_repeat_does_to_the_genome_what_it_does_to_the_transcript`,
  `a_genomic_repeat_projected_to_the_transcript_is_the_same_change`,
  `an_intronic_repeat_does_to_the_genome_what_it_does_to_the_pre_mrna` (properties).

### A position outside every exon is an error

**Transcript positions resolve through the exon structure only; a position past the last exon is a `ValidationError`,
not an extrapolated coordinate.**

A silently extrapolated coordinate is wrong without saying so.

**Example** (`NM_PLUS10.1` is one exon of 100 bases; c.\*60 is its last base):

```text
NM_PLUS10.1:c.*60T>A   →  projects
NM_PLUS10.1:c.*61A>G   →  ValidationError, from c_to_g, validate and to_spdi_unambiguous alike
```

- **Differs:** biocommons extrapolates in some paths.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "it is not allowed to
  describe variants in nucleotides beyond the boundaries of a reference sequence".
- **Tests:** `transcript_coordinates_test::a_position_in_no_exon_is_an_error_not_an_extrapolation`;
  `transcript_coordinates_round_trip` (property) covers every in-exon position.

### An open CDS end numbers no position

**A transcript whose CDS runs off one end of its record has no start or stop codon there. A position numbered from that
end — `c.N` and `c.-N` from an open 5' end, `c.*N` from an open 3' end — is an error in either direction of projection,
as is a protein consequence that depends on the open end.**

The provider says so with `cds_start_open` and `cds_end_open` on the transcript model. Without them a publisher's open
CDS can only be given as whole, which numbers from a codon the record does not carry, or withheld, which reads every
`c.` name as non-coding. With a 5' end open nothing in `p.` is predicted, since the frame is a fact about the missing
start. With a 3' end open an in-frame change wholly inside the whole codons the record carries, short of its last, is
predicted; a frameshift, or a change reaching that last codon, reads coding sequence the record does not have.

**Example** (`NM_OPEN3.1` carries its CDS's first ten codons and no stop; `NM_OPEN5.1` carries its last codons and no
start):

```text
NM_OPEN3.1:c.5C>T    →  NP_OPEN3.1:p.(Ala2Val)   in frame, inside what the record carries
NM_OPEN3.1:c.5del    →  ValidationError           a frameshift reads past the record's coding bases
NM_OPEN3.1:c.*1G>T   →  ValidationError           no stop codon to number from
NM_OPEN5.1:c.1A>T    →  ValidationError           no start codon to number from
NM_OPEN5.1:c.*1G>T   →  NC_…:g.33G>T              the stop codon is there
NM_OPEN5.1:n.3A>T    →  NC_…:g.3A>T               n. numbering needs no codon
```

- **Others:** not checked.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): `c.` positions are
  numbered from the A of the start codon and `c.*` from the stop; neither is defined without the codon.
- **Tests:** `open_cds_test` (every row above, and the `g.` to `c.` direction); `tests/test_open_cds.py`.

### Intronic positions are carried as given

**An intronic edit has no transcript base, so it is kept as written: projected to the genome it is carried across, and
normalisation leaves it where it is. The bases it states are read on the genome, on the transcript's strand
([below](#a-stated-base-the-reference-lacks-is-an-error)).**

**Example** (`NM_SPLICED.1`, whose intron is genome\[50..60) on `NC_D.1`, where g.53 is A):

```text
NC_D.1:g.53A>G          →  NM_SPLICED.1:c.40+3A>G      carried as given
NM_SPLICED.1:c.40+3A>G  →  NC_D.1:g.53A>G              and back

normalize NM_PLUS10.1:c.30_30+1insA   →  NM_PLUS10.1:c.30_30+1insA    unchanged, not an error
normalize NM_MINUS10.1:c.1+5del       →  NM_MINUS10.1:c.1+5del
```

- **Others:** not checked.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): RNA and coding reference
  sequences "do not contain intron sequences and can therefore not be used to describe variants affecting these
  sequences".
- **Tests:** `projection_reference_test::an_intronic_position_has_no_transcript_base_to_re_read`;
  `transcript_coordinates_test::normalize_leaves_intronic_coding_variants_as_written`;
  `intron_offsets_agree_from_either_exon` (property).

### Mitochondrial and RNA descriptions share the DNA machinery

**`m.` is `g.` on the mitochondrial reference; `r.` is `c.` or `n.` in RNA letters.**

`m.` and `g.` differ only in the letter they are written with: normalisation, SPDI, validation and equivalence are
shared. Every `r.` operation is a conversion to the `c.` (coding transcript) or `n.` (non-coding) spelling and back; the
statements `r.0`, `r.spl`, `r.?` and `r.(=)` map to `p.0?`, `p.?`, `p.?` and `p.(=)`.

`r.(?)` is a statement of its own, not `r.?` predicted: its parentheses belong to the statement, where those of `r.(=)`
mark a prediction. The uncertain page reads it as "RNA was not analysed but … no changes on the RNA level are expected
other than those simply derived from the change on the DNA level", and the deletion page pairs it with that change's own
consequence, `c.1521_1523del r.(?) p.(Phe508del)`. It says the RNA follows a DNA change it does not carry, so weaver
gives it no protein consequence (`r_to_p` is an `UnsupportedOperation`, where `r.?` gives `p.?`), and in equivalence it
is the same as itself and different from `r.?`, alone or among the members of a cis allele. No statement about the RNA
is predicted, so a dict that marks `r.?` predicted is refused rather than read as `r.(?)`.

**Example:**

```text
normalize NC_TEST.1:m.1012_1013insT   →  NC_TEST.1:m.1012dup      as g. would, written back as m.

NM_R.1:r.-3_5delinsuu   ⇄  NM_R.1:c.-3_5delinsTT
NR_R.1:r.10c>g          ⇄  NR_R.1:n.10C>G
NM_R.1:r.(10c>g)        →  NM_R.1:c.10C>G                        the prediction flag is dropped on the way
NM_R.1:r.spl            →  NP_R.1:p.?
NM_R.1:r.0              →  NP_R.1:p.0?
NM_W.1:r.(?)            →  p.: an error; the protein follows the DNA change, not this statement
NM_W.1:r.(?)  vs  NM_W.1:r.?                                    Different
```

- **Differs:** VariantValidator accepts `r.` in exons but rejects `r.*10`, intronic `r.`, `r.spl` and `r.0`, all of
  which HGVS allows.
- **Spec (agrees; silent on what `r.(?)` implies for the protein and how it compares):**
  [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "nucleotide `r.123` relates to `c.123` or
  `n.123`"; [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/) for `r.0`,
  `r.spl`, `r.?`, `r.=`; [uncertain](https://hgvs-nomenclature.org/stable/recommendations/uncertain/) and
  [protein deletion](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/) for `r.(?)`. The readings
  it is silent on are weaver's.
- **Tests:** `transcript_coordinates_test::mitochondrial_variants_share_the_genomic_implementation`;
  `rna_test::r_is_c_in_rna_letters`, `::r_on_a_non_coding_transcript_is_n`,
  `::statements_about_the_transcript_round_trip`, `::r_predicts_the_protein_like_c`,
  `::r_normalises_validates_and_has_alleles_like_c`,
  `::an_r_change_on_a_transcript_the_provider_lacks_is_not_compared_even_with_itself`;
  `rna_is_the_transcript_in_other_letters` (property); `statements_test::r_follows_the_dna_is_its_own_statement`,
  `::r_follows_the_dna_prints_back_as_parsed`, `::r_follows_the_dna_survives_json`,
  `::a_one_member_cis_allele_of_r_follows_the_dna_is_not_the_unknown_effect`,
  `::a_cis_allele_holding_r_follows_the_dna_is_not_one_holding_the_unknown_effect`,
  `::r_follows_the_dna_has_no_protein_consequence_of_its_own`, `::r_follows_the_dna_is_not_the_unknown_effect`,
  `::one_member_cis_alleles_of_r_follows_the_dna_and_the_unknown_effect_differ`, `::r_follows_the_dna_is_itself`,
  `::r_follows_the_dna_is_its_one_member_cis_allele`, `::a_dict_for_r_follows_the_dna_names_its_statement`,
  `::a_dict_for_r_unknown_marked_predicted_is_refused`, `::a_dict_for_a_special_unknown_marked_predicted_is_refused`;
  `tests/test_statements.py::test_a_dict_names_r_follows_the_dna_as_its_own_statement_not_a_prediction`.

### The whole transcript unchanged is a statement, not a change

**`c.=` says the whole coding transcript was analysed and no change was found. weaver reads it as a statement about the
transcript: it predicts the transcript's own protein unchanged, and has no genomic form, no allele and no `r.`
spelling.**

HGVS writes `c.=` for a coding transcript with no variant anywhere in it, and MaveDB's variant grammar accepts it. It
names no position, so nothing that works from positions can take it: the transcript is spliced from many exons, so there
is no one stretch of genome to project it onto, and there are no bases to make an allele, an SPDI or a VRS identifier
from. weaver does not invent any of these (for example, a reference allele over the transcript's span). It answers only
what follows from the statement:

- **Protein.** If nothing changed in the transcript, the protein it encodes is predicted unchanged, so `c.=` gives
  `p.(=)`. This holds even where the record lacks the start codon (a CDS open at its 5' end). The statement covers only
  the protein the provider pairs with the transcript: naming another protein is an error, and so is `c.=` on a
  transcript without a CDS. In equivalence `c.=` is the same as `p.(=)` on that protein, analogous to `p.=`, and not
  equivalent to either on any other protein. A `p.(=)` names no residues, so `c.=` has no protein allele.
- **No `r.` spelling.** `c.=` and `r.=` make different claims: one says the DNA was analysed, the other the RNA. So
  `c_to_r` refuses `c.=` and `r_to_c` refuses `r.=`, as both refuse `r.spl`. For the same reason the two are not
  equivalent. They do not contradict each other, but neither implies the other. `r.=` and `r.(=)` keep their own
  readings, [above](#mitochondrial-and-rna-descriptions-share-the-dna-machinery).
- **Validation.** `c.=` names no bases, so `validate` has nothing to compare. It is true once the transcript is found
  and has a CDS; on a non-coding transcript it is an error. Normalising leaves `c.=` as written, once its transcript is
  found.
- **Alleles in cis.** `c.[=]` alone is `c.=`, the half that the trans form `c.[2376G>C];[=]` names. Beside a change it
  contradicts that change, so `c.[=;10C>G]` is refused when parsed. The same rule refuses `r.[=;10c>g]` and
  `p.[=;Arg2Gly]`. The other statements keep their place: `r.[spl;10c>g]` parses. An identity at a position is not this
  statement but an edit that leaves its bases as they are, so `c.[10=;12C>G]` is accepted and makes what `c.12C>G` makes
  ([cis alleles that overlap](#overlapping-members-of-a-cis-allele-need-one-reading), where the departure from the
  alleles page is recorded).
- **Refused.** A genomic form, a canonical allele, SPDI and VRS are an `UnsupportedOperation` error.

Only `c.` takes the whole-transcript identity. The specification never writes `n.=` and defines no whole-sequence
statement for `g.` or `m.`, and `c.` has no predicted form, so `n.=`, `g.=`, `m.=` and `c.(=)` do not parse. A
dictionary that describes one of them is refused by `Variant.from_dict` too. The identity has its own edit, so a
dictionary for `r.=` names it rather than a `Special` "=" string; one written with the string loads as the identity, as
the other statements written that way load as theirs.

`r_to_p` holds every `r.` statement to the CDS condition `c_to_p` holds `c.=` to: on a transcript without a CDS, `r.=`,
`r.spl` and `r.0` are an error rather than a statement about a protein. It does not check that a named protein is the
transcript's own.

**Example** (`NM_W.1` encodes `NP_W.1`; `NR_W.1` is non-coding):

```text
NM_W.1:c.=   →  NP_W.1:p.(=)
NM_W.1:c.=   →  NP_OTHER.1:  an error: NM_W.1 does not encode it
NR_W.1:c.=   →  p.:          an error: no CDS
NM_W.1:c.=   →  c_to_g, c_to_r, to_spdi, to_vrs: an error
NM_W.1:c.=   vs NM_W.1:c.[=]     Identity
NM_W.1:c.=   vs NP_W.1:p.(=)     Identity
NM_W.1:c.=   vs NP_W.1:p.=       Analogous
NM_W.1:c.=   vs NM_W.1:r.=       Different: DNA analysed is not RNA analysed
NM_W.1:c.=   vs NM_W.1:c.10=     Different: a statement about c.10 is not one about the whole transcript
NM_W.1:c.[=;10C>G]               refused when parsed
NR_W.1:n.=, NC_W.1:g.=, NM_W.1:c.(=)   do not parse
```

- **Differs:** biocommons rejects `c.=`, `c.[=]` and `r.=`, and parses `p.=`. VariantValidator: not checked. The ClinVar
  sample has no whole-sequence statement.
- **Spec (agrees; silent on the protein consequence of `c.=`, its genomic form and allele, `n.=`, and how `c.=` relates
  to `r.=`):** [DNA substitution](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/): "`c.=`
  indicates the entire `NM_004006.2` coding DNA reference sequence was analysed and no change was identified";
  [DNA alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) for `c.[2376G>C];[=]`; the
  specification writes `[=]` only as a whole allele and never beside a change in one. The readings it is silent on are
  weaver's. It writes `r.(=)` for "no changes on the RNA level are expected (RNA not analysed)"
  ([uncertain](https://hgvs-nomenclature.org/stable/recommendations/uncertain/)), which is a different claim from `c.=`.
- **Tests:** `whole_sequence_identity_test::c_identity_prints_back_as_parsed`,
  `::c_identity_is_its_one_member_cis_allele`, `::c_identity_is_the_predicted_protein_unchanged_on_its_own_protein`,
  `::an_r_identity_serialised_as_a_statement_string_still_loads`, `::n_identity_does_not_parse`,
  `::a_predicted_c_identity_does_not_parse`, `::c_identity_beside_a_change_is_refused`,
  `::r_splicing_beside_a_change_parses`, `::c_identity_predicts_its_transcripts_protein_unchanged`,
  `::c_identity_says_nothing_about_a_protein_its_transcript_does_not_encode`,
  `::c_identity_on_a_transcript_without_a_cds_has_no_protein_consequence`,
  `::c_identity_predicts_the_protein_unchanged_on_a_cds_open_at_its_start`, `::c_identity_has_no_r_spelling`,
  `::c_identity_has_no_genomic_form`, `::c_identity_has_no_vrs_allele`, `::c_identity_on_a_coding_transcript_validates`,
  `::c_identity_is_not_r_identity`, `::c_identity_is_analogous_to_its_own_protein_unchanged`,
  `::c_identity_is_not_another_protein_unchanged`, `::c_identity_is_not_one_base_unchanged`,
  `::r_identity_on_a_transcript_without_a_cds_has_no_protein_consequence`, `::a_deserialised_g_identity_is_refused`;
  `tests/test_mavehgvs_compliance.py`, `tests/test_whole_identity.py`.

### An allele in cis projects member by member

**`c.[a;b]` projects to the genome, and `g.[a;b]` to a transcript, one member at a time, by the rules a single variant
follows, and the result is an allele in cis of the projections, still an allele when it has one member. Across strands
the members are written in reverse. If any member fails to project, the allele does, and the error names the member.**

HGVS writes each member of an allele as it would be written alone, numbered against the reference: `c.[2376G>C;3103del]`
contains "two different changes, `c.2376G>C` and `c.3103del`". Nothing one member does moves another's position, so each
can be projected on its own, and the projections together describe the same chromosome. Projection keeps the description
as written, as it does for a single variant ([no normalisation](#no-normalisation-before-projecting)): members that meet
are not merged, and a one-member allele stays one, `c.[7C>T]` becoming `g.[22C>T]`, not `g.22C>T`. Which spelling to
hand on is normalisation's question, not projection's.

The members' order is kept as the writer read it along their strand. HGVS lists members "in genomic order", and on a
minus-strand transcript that order runs the other way along the genome, so the projection reverses the list: an allele
written in order on the transcript is in order on the genome, and projecting it back returns the members as they were,
unless a duplication was rewritten as an insertion (below). The alternative, sorting by position, would write a
recommended order whatever the input, but it would reorder members the writer placed, and the round trip would return
them in another order.

All or nothing follows from what an allele is. A chromosome carrying two changes, one of which has no genomic form, is
not described by the genomic form of the other, so returning the members that did project would describe a different
allele. The error is the failing member's own, an `UnsupportedOperation` for an `r.` member across a splice junction, a
`ValidationError` for a position outside every exon, wrapped with the member and its place in the allele, so a caller
can still match on why. A member that projects to a form the target cannot hold fails the same way: an `r.` change
across a splice junction has [no genomic form](#a-change-across-a-splice-junction-has-no-genomic-form), and neither has
an allele containing one.

One member reads differently alone. A duplication inserts its copy "directly 3' of the original copy", on the strand it
is written on, and on the other strand 3' is the other side. Alone that does not matter: a copy beside the original
makes the same sequence on either side of it. Beside an insertion it does. On `NR_XR.1`, `n.[22_23insA;23_25dup]` puts
the A before the copied `AGG` and the copy after it. Projected one at a time, the duplication is `g.26_28dup`, whose
copy goes between g.28 and g.29, and the insertion is `g.28_29insT`, between the same two bases: two insertions at one
junction, whose order the genome's description cannot say, though the transcript's could. So a duplication whose
projected form would insert where another member inserts is written as the insertion it is, at the junction it uses
where it was written: `g.25_26insCCT`. An insertion of a stated length, `insN[5]`, meets a copy as one of stated bases
does. The insertion written can itself land where another duplication's projected form would insert, so the rewriting
repeats until no duplication meets an insertion; which members are rewritten does not depend on the order they are
written in. The same holds from the genome to a minus-strand transcript. A duplication with nothing inserted at its
junction stays a duplication. Projected back, a rewritten duplication comes back as the insertion, not as a duplication:
the same sequence, another spelling. Where a duplication or an insertion in such an allele has an uncertain position,
whether the copy meets the insertion cannot be told, so the projection across strands is refused rather than guessed.

**Example** (`NM_X.1` is `GGGGG` `ATG AAA CTG GCC TAT CGC TAA` `CCGTATAAGTAAGG` at genomic index 10 of `NC_X.1`, so c.7
is g.22; `NR_XR.1` reads the same forty bases of `NC_X.1` on the minus strand, so n.23 is g.28 and n.29 is g.22, and
`NM_XR.1` is `NR_XR.1` with its CDS from n.4; `NM_S.1` is `NM_X.1` with exon 1 ending at c.5 and a ten-base intron after
it in `NC_S.1`, so c.4 is g.19 and c.6 is g.31):

```text
NM_X.1:c.[7C>T;13T>G]              →  NC_X.1:g.[22C>T;28T>G]
NC_X.1:g.[22C>T;28T>G]             →  NM_X.1:c.[7C>T;13T>G]
NR_XR.1:n.[23A>C;29G>A]            →  NC_X.1:g.[22C>T;28T>G]           the minus strand reverses the members
NC_X.1:g.[22C>T;28T>G]             →  NM_XR.1:c.[20A>C;26G>A]          and back
NM_S.1:c.[4A>T;6A>T]               →  NC_S.1:g.[19A>T;31A>T]           members in two exons
NM_S.1:c.[5+2T>A;7C>T]             →  NC_S.1:g.[22T>A;32C>T]           an intronic member, carried as given
NM_X.1:c.[7C>T]                    →  NC_X.1:g.[22C>T]                 one member, still an allele
NM_X.1:c.[7C>T;13T>G]              →  NM_X.1:r.[7c>u;13u>g]

NR_XR.1:n.[23_25dup;29G>A]         →  NC_X.1:g.[22C>T;26_28dup]        nothing inserted beside the copy
NR_XR.1:n.[22_23insA;23_25dup]     →  NC_X.1:g.[25_26insCCT;28_29insT] the copy written where the transcript puts it
NC_X.1:g.[25_26insA;26_28dup]      →  NM_XR.1:c.[19_20insAGG;22_23insT]
NR_XR.1:n.[26_28dup;22_23insA;23_25dup]
                                   →  NC_X.1:g.[25_26insCCT;28_29insT;22_23insTGG]
                                                                       the second copy meets the first, rewritten

NM_S.1:r.[4a>u;5_6del]             →  UnsupportedOperation, member 2, NM_S.1:r.5_6del: spans a splice junction
NM_X.1:c.[7C>T;*20del]             →  ValidationError, member 2, NM_X.1:c.*20del: outside every exon
```

- **Others:** biocommons 1.5.7 does not parse alleles in cis; VariantValidator and ClinVar not checked.
- **Spec (agrees; silent on a duplication whose side of the original the strands disagree on):**
  [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/): an allele is "a series of variants on
  one chromosome", each member written as alone, listed "in genomic order";
  [duplication](https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/): "a copy of one or more
  nucleotides is inserted directly 3' of the original copy". Writing a duplication as an insertion where the strands
  disagree on its side is weaver's reading; the specification does not address it.
- **Tests:** `cis_projection_test::members_project_to_the_genome_one_by_one`,
  `::members_project_from_the_genome_one_by_one`, `::the_minus_strand_reverses_the_members_on_the_genome`,
  `::the_minus_strand_reverses_the_members_on_the_transcript`, `::members_in_different_exons_project_to_each`,
  `::a_member_in_an_intron_projects_as_given`, `::a_one_member_allele_stays_an_allele`,
  `::a_coding_allele_is_respelled_in_rna_and_back`, `::a_non_coding_allele_is_respelled_in_rna_and_back`,
  `::an_rna_allele_projects_to_the_genome`, `::every_transcript_over_the_members_gets_the_allele`,
  `::a_minus_strand_duplication_alone_stays_a_duplication`,
  `::a_minus_strand_duplication_beside_an_insertion_is_written_as_the_insertion`,
  `::a_genomic_duplication_beside_an_insertion_is_written_as_the_insertion_on_a_minus_strand_transcript`,
  `::a_coding_minus_strand_duplication_beside_an_insertion_is_written_as_the_insertion`,
  `::an_rna_minus_strand_duplication_beside_an_insertion_is_written_as_the_insertion`,
  `::an_insertion_of_a_length_beside_a_duplication_is_met_like_one_of_bases`,
  `::a_rewritten_duplication_can_meet_the_next_whatever_the_order_written`,
  `::a_rewritten_genomic_duplication_can_meet_the_next`,
  `::a_rewritten_duplication_written_in_another_order_meets_the_next`,
  `::an_uncertain_duplication_beside_an_insertion_across_strands_is_refused`,
  `::an_uncertain_transcript_duplication_beside_an_insertion_is_refused`,
  `::a_copy_of_the_genomes_first_base_has_no_junction_before_it`,
  `::a_member_that_fails_to_project_fails_the_allele_and_is_named`, `::a_member_outside_every_exon_fails_the_allele`;
  `a_cis_allele_projected_to_the_genome_and_back_is_itself`,
  `a_duplication_beside_an_insertion_makes_the_same_exon_on_the_genome`,
  `a_cis_allele_in_transcript_order_projects_in_genomic_order` (properties); `tests/test_cis_phased.py`.

## Normalisation

### Shift 3', cyclically over repeats

**A deletion, duplication or insertion is shifted as far 3' as the reference keeps matching it cyclically, on the strand
of the sequence the variant is written on.** On a transcript the shift stays within one exon
([below](#a-shift-stops-at-an-exonexon-junction)); a delins is first reduced to the edit it is
([below](#a-delins-is-reduced-before-it-is-shifted)).

An insertion that repeats the bases immediately before it becomes a duplication, and an insertion that slides has its
bases rotated by how far it moved. An insertion does not slide past the last base of the sequence, after which HGVS has
no position to name; one written or read there, from SPDI or VRS, is a delins of that base, or the duplication it is.

**Example** (on `ACGT` repeated, and on the ten bases `TTCAGCAGTT`):

```text
NC_TEST.1:g.1012_1013insT      →  NC_TEST.1:g.1012dup            repeats the T before it
NC_TEST.1:g.1012_1013insACGT   →  NC_TEST.1:g.1997_2000dup       slides to the end of the repeat
NC_TEST.1:g.1012_1013insAA     →  NC_TEST.1:g.1013_1014insAA     slides one base past an A, still an insertion
NM_PLUS0.1:n.12_13insT         →  NM_PLUS0.1:n.12dup             n. and g. alike

on TTCAGCAGTT:  g.3_5del       →  g.6_8del                       the first CAG shifts onto the second
on GGGGGGGGG:   g.2_3insGA     →  g.3_4insAG                     slid one base, so rotated by one
on TTCAGCAGTT:  g.9_10insTA    →  g.9_10insTA                    not g.10_11insAT, after the last base
on TTCAGCAGTT:  g.10_11insA    →  g.10delinsTA                   and NC_…:9:T:TA reads as the same
```

- **Agrees:** biocommons, VariantValidator and ClinVar.
- **Spec (agrees):** [general, 3' rule](https://hgvs-nomenclature.org/stable/recommendations/general/);
  [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/): "tandem duplications are described
  as a duplication, not an insertion".
- **Tests:** `one_spelling_test::an_insertion_does_not_slide_past_the_last_base_of_the_sequence`,
  `::an_insertion_written_after_the_last_base_is_a_delins_of_the_last_base`,
  `::an_insertion_written_after_the_last_base_that_repeats_it_is_a_duplication`,
  `::an_spdi_insertion_after_the_last_base_reads_as_the_normalised_delins`,
  `::an_spdi_insertion_after_the_last_residue_reads_as_a_delins_of_the_last_residue`;
  `src/normalize.rs::deletion_shifts_3_prime_and_states_no_bases`,
  `::insertion_that_repeats_preceding_bases_becomes_duplication`, `::a_slid_insertion_is_rotated`,
  `::shift_5_mirrors_shift_3` (the last two rows above, read from the 0-based ranges those tests assert);
  `transcript_coordinates_test::normalize_converts_genomic_and_noncoding_insertions_to_duplications`;
  `normalising_preserves_the_edited_sequence` and `ambiguous_range_is_sound_and_complete` (properties).

### A shift stops at an exon/exon junction

**On a transcript, a deletion, duplication or insertion shifts 3' only as far as the last base of the exon it is written
in. No end of it moves into another exon.**

A `c.`, `n.` or `r.` variant is normalised on the transcript's own sequence, which is spliced: the last base of one exon
and the first base of the next are neighbours there, but an intron separates them on the genome. Where the bases either
side of a junction are the same, the transcript cannot tell `c.4del` from `c.6del`, yet on the genome they delete
different bases, one in each exon. Shifting across the junction would rewrite the variant as another change.

The rule reaches every kind of edit that slides:

- An insertion does not slide onto the junction. Between the two exons, `c.5_6ins` names no single place on the genome:
  the end of one exon or the start of the next. An insertion written there is left as written.
- An insertion that repeats the bases before it becomes a duplication only when those bases lie in one exon. A
  duplication of bases from two exons would duplicate the intron with them.
- A deletion or duplication that already spans a junction may shift while each end stays in its exon. Each spelling
  covers the whole intron on the genome, and two spellings one base apart differ only in their first and last bases,
  which are the same base, so on the genome they are one change.
- The transcript's last base is the same kind of edge: an insertion does not slide past it.

A shift still crosses the start or the end of the CDS. HGVS gives no exception there, and within one exon the genome
change is the same either side.

An exon aligned to the genome with gaps ([above](#gapped-exons-are-projected-locally)) is one stretch to the shift,
which can therefore still cross a gap in its alignment, as it can cross a base where the record and the genome differ.

On the ClinVar sample, the descriptions this decides are `NM_000069.3:c.5042_5048delTTTCCAG` and
`NM_001354663.2:c.2449delA`. Each ends on the last base of an exon, the base across the junction repeats the deletion's
first base, and each normalises to itself, with the SPDI ClinVar gives; shifted into the next exon, the first would
delete the intron with it and the second would land 25 kb from ClinVar's position.

**Example** (exon 1 of `NM_J.1` ends at c.5; c.4, c.5 and c.6 are all A, c.3 and c.7 both G; c.5 is g.20 and c.6 is
g.31):

```text
NM_J.1:c.4del         →  NM_J.1:c.5del          g.20del, as c.4del is; c.6del would be g.31del
NM_J.1:c.4dup         →  NM_J.1:c.5dup
NM_J.1:c.3_4insA      →  NM_J.1:c.5dup          slides to the junction and duplicates the base before it
NM_J.1:c.4_5insAC     →  NM_J.1:c.4_5insAC      not c.5_6insCA, on the junction
NM_J.1:c.5_6insAA     →  NM_J.1:c.5_6insAA      written on the junction, so left as written
NM_J.1:c.6_7insAAA    →  NM_J.1:c.6_7insAAA     not c.4_6dup, across two exons
NM_J.1:c.3_6del       →  NM_J.1:c.4_7del        spans the junction; each end stays in its exon
NM_J.1:r.4del         →  NM_J.1:r.5del          r. and n. alike
NM_JM.1:c.4del        →  NM_JM.1:c.5del         the minus strand alike
```

- **Agrees:** biocommons, run with its default of not crossing boundaries, on the deletion, the duplication, the
  insertion that becomes a duplication, `r.`, and both strands; ClinVar, on the two rows above.
- **Differs:** biocommons slides `c.4_5insAC` onto the junction as `c.5_6insCA`; refuses `c.5_6insAA` and `c.3_6del` as
  spanning an exon boundary; and also stops a shift at the CDS start and end. VariantValidator not checked.
- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "the 3' rule is not
  applied when there is a deletion/duplication around exon/exon junctions with identical nucleotides flanking the
  junction, where shifting the variant 3' would place it in the next exon". The
  [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/) page's example is
  `LRG_199t1:c.3921del`, "not described as c.3922del since this would shift the position of the variant to the next
  exon", and the [duplication](https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/) page's is
  `NC_000023.11(NM_004006.2):c.3921dup`, not `c.3922dup`. The specification names deletions and duplications only;
  holding insertions to the same line is weaver's reading of its reason. It gives no exception at the CDS bounds.
- **Tests:** `exon_junction_test::a_deletion_shifts_to_the_last_base_of_its_exon_and_no_further`,
  `::a_deletion_normalised_short_of_the_junction_projects_where_it_was_written`,
  `::a_duplication_shifts_to_the_last_base_of_its_exon_and_no_further`,
  `::an_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon`,
  `::an_insertion_does_not_slide_onto_the_junction`, `::an_insertion_written_on_the_junction_is_left_as_written`,
  `::an_insertion_is_not_rewritten_as_a_duplication_across_the_junction`,
  `::a_deletion_spanning_the_junction_shifts_while_each_end_stays_in_its_exon`,
  `::a_deletion_spanning_the_junction_shifts_no_further_than_its_start_can_stay_in_its_exon`,
  `::a_duplication_spanning_the_junction_shifts_while_each_end_stays_in_its_exon`,
  `::a_deletion_beside_the_junction_does_not_shift_to_span_it`,
  `::an_insertion_does_not_slide_past_the_last_base_of_the_transcript`,
  `::a_shift_crosses_the_cds_start_within_one_exon`, `::a_shift_crosses_the_cds_end_within_one_exon`,
  `::a_non_coding_deletion_stops_at_the_junction`,
  `::a_non_coding_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon`,
  `::an_rna_deletion_stops_at_the_junction`,
  `::an_rna_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon`,
  `::on_the_minus_strand_a_deletion_stops_at_the_junction`,
  `::on_the_minus_strand_a_deletion_normalised_short_of_the_junction_projects_where_it_was_written`,
  `::on_the_minus_strand_a_duplication_stops_at_the_junction`,
  `::on_the_minus_strand_an_insertion_that_slides_to_the_junction_duplicates_the_last_base_of_its_exon`,
  `::on_the_minus_strand_an_insertion_does_not_slide_onto_the_junction`;
  `normalising_on_the_transcript_keeps_the_genomic_allele` and `normalising_beside_a_junction_keeps_the_genomic_allele`
  (properties).

### A delins is reduced before it is shifted

**`normalize` drops the bases a delins shares with the reference at either end, and writes what remains as the edit it
is: a substitution, deletion, insertion, duplication, inversion or delins. That edit is then normalised as usual,
shifted 3' within its exon. A delins that changes a base at each end cannot slide, and stays where it is.**

One change has more than one spelling as a delins: `delCinsCA` and `insA` put the same A after the same C, and a delins
written over a wider range than the change says the same thing again. The 3' rule is stated for every description, and a
delins that is really an insertion has a 3'-most position like any insertion. Reducing first gives each change one
spelling. Deleted bases a delins states are checked against the reference first
([below](#a-stated-base-the-reference-lacks-is-an-error)), and are not written back
([below](#deletions-and-duplications-are-written-bare)).

On a transcript, the reduction respects the exons as a shift does ([above](#a-shift-stops-at-an-exonexon-junction)):

- A delins whose range spans a junction removes the intron on the genome. It keeps at least one base in each exon, so
  what it reduces to removes the intron too, and it is never rewritten as an inversion, which would turn the intron
  round instead.
- A delins inside one exon that reduces to an insertion at the edge of that exon would put the insertion on the
  junction, where it names no single place on the genome. It is written as a delins of the exon's base beside it
  instead, unless the insertion can sit one base further 5' inside the exon, where an insertion slid to the junction
  stops. biocommons writes both `c.4_5delinsAAC` and `c.6delinsCA` below as `c.5_6insC`, one spelling for two changes on
  either side of the intron.
- The last base of a genome is the same kind of edge: HGVS cannot name a position after it.

A position known only to lie in a range, `g.(3_4)`, names no one placement to reduce or shift from, so a variant with
one is left as written.

The sequence a description makes is unchanged: the property tests hold normalisation to the edited sequence, to the
genomic allele across junctions, and to one output for a change and a wider delins spelling it.

A delins is never slid as a delins. Sliding one with the rule for a deletion, while the base after the range repeats the
first deleted base, moves the inserted bases past that base and changes the sequence. Of the 429 delins in the
100,000-variant ClinVar sample, the 40 that rule would slide all change a base at each end and normalise to themselves,
and every delins normalises to the SPDI ClinVar gives.

**Example** (on `TTCAGCAGTT`, and on `NM_J.1`, where exon 1 ends at c.5, c.3 is G, c.4 to c.6 are A and c.7 is G; c.5 is
g.20 and c.6 is g.31):

```text
g.3delinsCGG             →  g.3_4insGG               C becomes CGG: an insertion after the C
g.3delinsCAA             →  g.4_5insAA               the insertion, slid past the A at g.4
g.6delinsCAG             →  g.7_8dup                 slid to the end of the CAG, which it repeats
g.2_5delinsT             →  g.6_8del                 TCAG becomes T: the CAG deleted, the 3'-most copy
g.3_5delinsCTG           →  g.4A>T
g.3_4delinsTG            →  g.3_4inv                 CA becomes its reverse complement
g.3_5delinsTT            →  g.3_5delinsTT            a base changed at each end: it stays
g.3_5delinsCAG           →  g.3_5=

NM_J.1:c.3delinsGA       →  NM_J.1:c.5dup            as c.3_4insA does
NM_J.1:c.4_5delinsAAC    →  NM_J.1:c.5delinsAC       g.20_21insC; c.5_6insC would be on the junction
NM_J.1:c.6delinsCA       →  NM_J.1:c.6delinsCA       g.30_31insC, after the intron
NM_J.1:c.4_7delinsAG     →  NM_J.1:c.5_6del          g.20_31del either way; each end stays in its exon
NM_J.1:c.5_6delinsAGA    →  NM_J.1:c.5_6delinsAGA    removes the intron, which c.5_6insG would not
```

- **Agrees:** biocommons on every genomic row, `c.3delinsGA` included.
- **Differs:** biocommons writes `c.4_5delinsAAC` and `c.6delinsCA` as `c.5_6insC`, and refuses `c.4_7delinsAG` and
  `c.5_6delinsAGA` as spanning an exon boundary. VariantValidator not checked.
- **Spec (agrees):** [deletion-insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/): a delins is
  "a sequence change where, compared to a reference sequence, one or more nucleotides are replaced by one or more other
  nucleotides and which is not a substitution or inversion"; "by definition, when one nucleotide is replaced by one
  other nucleotide, the change is a substitution"; and "for all descriptions, the most 3' position possible of the
  reference sequence is arbitrarily assigned to have been changed (3'rule)". Keeping a delins on the junction's edge
  rather than writing an insertion there is weaver's reading of the exon/exon exception
  ([above](#a-shift-stops-at-an-exonexon-junction)).
- **Tests:** `one_spelling_test::a_delins_that_keeps_its_first_base_is_an_insertion`,
  `::a_delins_that_is_an_insertion_is_shifted_3_prime`,
  `::a_delins_that_is_an_insertion_repeating_the_bases_before_it_is_a_duplication`,
  `::a_delins_that_keeps_its_last_base_is_an_insertion`, `::a_delins_that_is_a_deletion_is_shifted_3_prime`,
  `::a_delins_that_changes_one_base_is_a_substitution`, `::a_delins_that_is_the_reverse_complement_is_an_inversion`,
  `::a_delins_that_changes_a_base_at_each_end_stays_where_it_is`, `::a_delins_is_written_without_its_deleted_bases`,
  `::a_delins_that_changes_nothing_is_an_identity`,
  `::a_delins_reduced_to_an_insertion_after_the_last_base_keeps_the_last_base`,
  `::a_change_with_an_uncertain_position_is_left_as_written`, `::a_change_with_uncertain_ends_is_left_as_written`;
  `exon_junction_test::a_delins_that_is_an_insertion_shifts_to_the_last_base_of_its_exon_and_no_further`,
  `::a_delins_reduced_to_an_insertion_at_the_end_of_its_exon_keeps_the_exons_last_base`,
  `::a_delins_reduced_to_an_insertion_at_the_start_of_its_exon_keeps_the_exons_first_base`,
  `::a_delins_spanning_the_junction_is_trimmed_while_each_end_stays_in_its_exon`,
  `::a_delins_of_the_two_bases_either_side_of_the_junction_stays_where_it_is`,
  `::on_the_minus_strand_a_delins_reduced_to_an_insertion_at_the_end_of_its_exon_keeps_the_exons_last_base`,
  `::on_the_minus_strand_a_delins_that_is_an_insertion_stops_at_the_junction`,
  `::a_delins_with_an_uncertain_position_is_left_as_written`; `a_change_and_a_wider_delins_spelling_it_normalise_alike`,
  `normalising_preserves_the_edited_sequence`, `normalising_a_transcript_variant_twice_is_normalising_it_once` and
  `normalising_beside_a_junction_keeps_the_genomic_allele` (properties).

### Deletions and duplications are written bare

**`normalize` writes `c.306del`, not `c.306delC`, whether or not the edit moved. The same holds for every reference base
an edit can state: a duplication's, an inversion's, the deleted part of a delins and an identity's all come out bare. A
substitution keeps its base, which is part of how it is written. Parsing, projecting and the other operations keep what
was written.**

A deletion written with its bases and one written without are one change, and `normalize` gives one change one spelling.
Written with their bases or without, `g.6_8delCAG` and `g.6_8del` would otherwise each be normalised output, two
spellings of one change. The reference holds the bases, so writing them adds nothing, and they are not filled in either.
For the same reason a range of one base is written as the base: `g.6_6dup` is `g.6dup`.

The bases are dropped only once they have been checked: a stated base the reference lacks is an error in `normalize` as
everywhere else ([below](#a-stated-base-the-reference-lacks-is-an-error)), so `g.3_5delAAA` over `CAG` is refused, not
rewritten as the deletion of `CAG`.

**Example** (on `TTCAGCAGTT`):

```text
NM_TEST:c.4_5del               →  NM_TEST:c.5_6del               bare, shifted
g.3_5delCAG                    →  g.6_8del                       shifted
g.6_8delCAG                    →  g.6_8del                       not moved, and bare all the same
g.6_8del3                      →  g.6_8del                       a count is not written either
g.6_8dupCAG                    →  g.6_8dup
g.3_4invCA                     →  g.3_4inv
g.3C=                          →  g.3=
g.3_5CAG=                      →  g.3_5=
g.3_5delCAGinsTT               →  g.3_5delinsTT
g.3_5delCAGinsN[2]             →  g.3_5delinsN[2]
g.3_5delAAA                    →  ReferenceMismatch              the AAA it states is not there
g.3C>T                         →  g.3C>T                         a substitution keeps its base
```

- **Agrees:** biocommons, as it prints by default (`max_ref_length` 0; it fills the bases into the value it returns), on
  every genomic row but two; VariantValidator on deletions and duplications, not checked on the rest.
- **Differs:** biocommons, run with `validate=False`, normalises `g.3_5delAAA` to `g.3_5del` rather than refusing it,
  and does not parse `delinsN[2]`.
- **Spec (agrees):** [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/): "the recommendation
  is not to describe the variant as NC_000023.11:g.33344591delA, i.e., describe the deleted nucleotide sequence. This
  description is longer, it contains redundant information, and chances to make an error increases";
  [duplication](https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/), the same for `c.20dupT`;
  [substitution](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/): "nucleotides that have been
  tested and found not changed are described as c.123=, g.4567_4569=". Inversions and delins are written bare in every
  example the specification gives.
- **Tests:** `src/normalize.rs::deletion_shifts_3_prime_and_states_no_bases`;
  `one_spelling_test::a_deletion_that_does_not_move_is_written_without_its_bases`,
  `::a_deletion_written_with_its_length_is_written_without_it`,
  `::a_duplication_that_does_not_move_is_written_without_its_bases`, `::a_range_of_one_base_is_written_as_the_base`,
  `::an_inversion_is_written_without_its_bases`, `::an_identity_of_one_base_is_written_without_it`,
  `::an_identity_of_a_range_is_written_without_its_bases`, `::a_delins_is_written_without_its_deleted_bases`,
  `::a_delins_of_a_number_of_bases_is_written_without_its_deleted_bases`, `::a_substitution_keeps_its_stated_base`,
  `::a_deletion_stating_bases_the_reference_lacks_is_refused`,
  `::a_delins_whose_deleted_bases_the_reference_lacks_is_refused`,
  `::a_delins_of_a_number_of_bases_whose_deleted_bases_the_reference_lacks_is_refused`,
  `::an_inversion_stating_bases_the_reference_lacks_is_refused`,
  `::an_identity_stating_a_base_the_reference_lacks_is_refused`, `::one_base_stated_over_two_is_refused`;
  `transcript_coordinates_test::mitochondrial_variants_share_the_genomic_implementation`;
  `stating_the_bases_an_edit_replaces_does_not_change_its_normalised_spelling` (properties);
  `tests/test_api.py::test_normalization`.

### A repeat resolves to its whole run

**`c.10AC[3]` means every existing copy of the unit from c.10 becomes three copies, so its allele and protein reading
cover the run.**

Before this, 0 of 394 ClinVar repeat rows matched on SPDI; after, 386.

**Example** (`NC_REP.1` is `TTTTT GCCATT GCCATT GCCATT AAAAA…`; all four are one allele):

```text
NM_REP_MINUS.1:c.78AATGGC[4]
NC_REP.1:g.6GCCATT[4]           →  the same canonical allele: the run of three becomes four
NC_REP.1:g.18_23dup
NC_REP.1:g.23_24insGCCATT
```

- **Agrees:** ClinVar, on those rows.
- **Spec (agrees):** [repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/DNA/repeated/): the
  bracketed number "represents the total count of repeat units".
- **Tests:** `transcript_coordinates_test::repeat_on_the_minus_strand_projects_to_its_whole_run`;
  `protein_allele_test::substitution_insertion_delins_dup_and_repeat_resolve_on_the_protein`;
  `equivalence_cases_test::test_analogous_repeat_equivalence`, `::test_multi_unit_repeat_equivalence`.

### Members of a cis allele that meet are written as one change

**Normalising a cis allele merges the members that meet, and on `c.` the members in one codon, and writes each group as
the one change it makes, or, where the group makes the reference, as the identity of the bases its members name; it
repeats until no two of the results meet. Statements are set aside first: identity members, and every set of one or more
members that together make the reference, no smaller set of them doing so, whose span holds no base another member names
and which shares no member with another such set, two that share one making the allele ambiguous; each is written as an
identity and grouped with nothing. A member that meets no other is normalised alone, as a single variant is. Identities
that overlap or touch in one exon are written as one. Nothing is merged across an exon junction, and no shift crosses
one or enters a positional identity.**

`normalize` produces the preferred spelling for interchange: what weaver hands to other tools, and what string matching
uses. Parsing, validating, projecting and comparing keep the form as written, `to_vrs` normalises because it too hands
the allele to other tools, and weaver [does not normalise before projecting](#no-normalisation-before-projecting).
Merging is for the readers of that spelling. Consumers that match strings converge only on preferred forms. Tools that
read members one at a time mispredict, which is the specification's reason for the delins. Members that touch are read
differently by readers with other conventions for applying them, and a merged change is read one way. And `to_vrs`
renders the normalised allele, so the HGVS and VRS output share a shape. A caller who wants the structure as written, a
combinatorial screen recording that a position was wild type in a given combination, does not call `normalize`.

Normalising the members one at a time can break the allele. A member's 3' shift stops at the first base that differs
from what it moves over, and a neighbour may change that base, or sit inside the run the member moves along. On
`TTCAGCAGTT`, `g.[3_5del;7A>T]` is a valid allele, but its deletion shifted alone becomes `g.6_8del`, which contains
g.7, and `g.[6_8del;7A>T]` is [refused as an overlap](#overlapping-members-of-a-cis-allele-need-one-reading). On
`NM_X.1`, the deletion of an A in `c.[4del;7C>A]` stops at c.7, a C, which the other member makes an A: shifted one at a
time the allele is `c.[6del;7C>A]`, two members that together make what `c.7del` makes.

Grouping once, on the members as written, is not enough either, because what a group becomes can reach a member it did
not reach before. On `TCTCT`, the first two deletions of `g.[1del;2del;5del]` touch and together delete `TC`, which
justifies over the whole run and is written 3'-most as `g.4_5del`; the third, grouped with nothing, stays `g.5del`. The
result, `g.[4_5del;5del]`, is an overlap. On `CTCTT`, `g.[1del;2del;4del]` grouped once is `g.[3_4del;5del]`, which
makes the right sequence, but whose two members touch, so normalising it again gives `g.3_5del`.

So weaver normalises by what the members do together, and keeps going until nothing changes:

1. The [overlap rule](#overlapping-members-of-a-cis-allele-need-one-reading) is checked on the members as written.
1. Identical members other than insertions and duplications are one edit written twice, and are kept once.
1. Statements are set aside: every identity member, and every set of one or more members that together make the
   reference, no proper subset of it doing so, and whose span, from the first base its members name to the last, holds
   no base another member names. A single member can be such a set, `c.9_12inv` where c.9_12 is its own reverse
   complement. Such a set is written as the identity of the bases its members name. Two such sets that share a member
   make the allele [ambiguous](#an-ambiguous-description-is-an-error-that-lists-its-readings): which pair cancels
   decides where the statement sits, and each choice is a reading. What is set aside is a barrier, as any identity is
   (below).
1. A change that agrees with an identity on some of its bases is cut where the identity stands, and the identity kept:
   `c.[7=;7_8delinsCG]` is written `c.[7=;8T>G]`. The overlap rule accepts such an allele; only `normalize` writes it
   this way.
1. Each member has an extent: the bases it is written over, together with its fully justified region, the region its
   [canonical allele](#the-canonical-allele-is-fully-justified) spans, everywhere it could equally be written. An
   insertion is written over no base: its extent is the empty interval at its junction, which touches the bases on
   either side.
1. Two members meet when their extents overlap or touch, with no unchanged base between them. On `c.`, two members also
   meet when their regions share a codon of the CDS. Two members in different exons never meet, nor do two groups whose
   members lie in different exons.
1. Members that meet are merged into a group. A group is applied to the reference and what it makes is written as one
   change: the bases it shares with the reference at both ends are trimmed, and the remainder is a delins unless it
   reduces to a deletion, insertion, duplication, substitution or inversion, by the rule a variant
   [read back](#reading-back-gives-the-normalised-variant) from VRS or SPDI follows, shifted 3' as a single variant is.
   The group's extent is the justified region of that change. A group with no remainder makes the reference, and is
   written as the identity of the bases its members name (below).
1. Groups are judged again by the same test, and groups that meet are merged. Each round merges every set of groups
   joined by meeting, directly or through another, so the result does not depend on the order they are looked at. This
   repeats until no two groups meet; it ends, because each round merges.
1. A member in a group of its own is normalised as a single variant is: a deletion, duplication or insertion is
   [shifted 3'](#shift-3-cyclically-over-repeats), a substitution is left as written, and a delins is
   [reduced before it is shifted](#a-delins-is-reduced-before-it-is-shifted).
1. Identities, written or left by a group that makes the reference, that overlap or touch within one exon are written as
   one.

The result makes the sequence the allele makes, passes the overlap rule, and normalising it again changes nothing.

Touching counts because of where a shift stops: once the neighbour has changed the base that stopped it, the member's
change could run on, as the deletion in `c.[4del;7C>A]` runs on into c.7. Regions that only touch share no base, so the
overlap test alone would leave that allele as two members. The specification says the same of consecutive bases: they
"are described as deletion/insertion (delins) variants".

An extent holds the bases a member is written over because a substitution or delins left as written still occupies them.
On `ACACCCCGC`, `g.[4del;7_9delinsCGG]` deletes a C from the run g.4_7 and changes g.9: the delins, trimmed to what it
alters, is `g.9C>G`, which does not touch the run. Judged by that, the two would stay apart, the deletion would shift to
`g.7del`, and `g.[7del;7_9delinsCGG]` is an overlap. Judged by the delins as written, g.7_9, they meet and are written
`g.7_9delinsGG`.

The rejected alternative groups once and skips a shift that would create an overlap. It avoids the overlap, but leaves a
change short of its 3' position because of where its neighbour sits: on `TCTCT`, `g.[1del;2del;5del]` would become
`g.[1_2del;5del]`, the deletion of `TC` written at the 5' end of the run it could be written anywhere in, where the same
deletion alone is written at the 3' end.

The codon is for the tools that read the description. weaver's own protein consequence is read from the codons all the
members leave ([the cis allele entry](#a-cis-allele-is-the-sequence-its-members-make)), so `c.[4A>T;6A>T]` is
`p.(Lys2Tyr)` however it is written. A tool that predicts member by member reads `c.4A>T` as `p.(Lys2Ter)` and `c.6A>T`
as `p.(Lys2Asn)`, the conflict the delins page writes its exception against; written `c.4_6delinsTAT`, the description
gives that tool the right answer too. Members one base apart in two codons, `c.[6A>T;8T>G]` in codons 2 and 3, stay
apart, as the same page says. Grouping runs across both conditions: in `c.[4A>T;6A>T;7C>G]` the first two share a codon
and the third touches the second, so all three are one change.

The codon is a property of the coding sequence, so it applies only to `c.` members inside the CDS. A `g.` or `n.`
member, or one in a UTR or an intron, has none. A `g.` allele and its `c.` projection can therefore normalise to
different shapes, `g.[19A>T;21A>T]` staying two members while `c.[4A>T;6A>T]` becomes one delins. Equivalence compares
the sequences, so the two are still the same change.

Two details of the codon condition:

- **A member's codons come from its justified region**, not from the bases it is written over. A deletion in a run that
  crosses a codon boundary is otherwise in whichever codon it is written in: in `NM_X.1`'s `CTG GCC`, `c.9del` and
  `c.10del` are one deletion, and with `c.12C>T` in codon 4, read by their written bases, `c.[10del;12C>T]` would become
  `c.10_12delinsCT` while `c.[9del;12C>T]` stayed two members. The justified region, c.9_10 for either spelling, groups
  both.
- **An insertion that cannot move** has an empty region. One inside a codon, between its first and second or second and
  third base, belongs to that codon, and one between two codons belongs to neither. `c.[4A>T;5_6insG]` becomes
  `c.4_5delinsTAG`; `c.[5A>T;6_7insG]` stays two members.

Nothing is merged across an exon junction, whether the members touch there or share a codon the junction splits. A
merged change would span the junction, and a `c.` change across a junction
[projects over the intron](#a-change-across-a-splice-junction-has-no-genomic-form): it deletes the intron, which the
members did not. It would also hide where each change is. On `NM_S.1`, which is `NM_X.1` with exon 1 ending at c.5,
`c.[4A>T;6A>T]` projects to two substitutions, `g.[19A>T;31A>T]`, and `c.4_6delinsTAT` to `g.19_31delinsTAT`, a
different genome; `c.[5A>T;6A>T]`, which touch across the junction, stay two members, while `c.[6A>T;7C>G]`, both in
exon 2, become `c.6_7delinsTG`. The same holds for `r.`, which normalises as its `c.` or `n.` spelling, and to which the
codon rule applies as it does to `c.`.

No shift crosses an exon junction either, as [a single variant's does not](#a-shift-stops-at-an-exonexon-junction), and
a member's justified region stops there too. Shifted across it, `c.4del` in `NM_S.1`'s `AA|A` would become `c.6del`,
which the [genome comparison](#a-cis-allele-is-the-sequence-its-members-make) judges a different change, g.31 deleted
rather than g.19. Inside an allele the crossing shift would move a deletion to the other exon, `c.[4del;7C>A]` becoming
`c.[6del;7C>A]`, which normalises again to `c.7del`, or onto a neighbour, `c.[4del;6A>T]` becoming `c.[6del;6A>T]`, an
overlap. Stopped at the junction, they are `c.[5del;7C>A]` and `c.[5del;6A>T]`. Like an identity, the junction bounds
the duplication test too. An insertion's shift stops before the junction itself, so normalisation never writes an
insertion on one: `c.4_5insAC` stays `c.4_5insAC`, where one place further it would be `c.5_6insCA`, which projects to
either side of the intron. A duplication's shift stops at the exon's last base, as a deletion's does
([the junction entry](#a-shift-stops-at-an-exonexon-junction)), and that is allowed inside an allele too:
`c.[4dup;6A>T]` is written `c.[5dup;6A>T]`. The copy goes on the junction, but directly 3' of the base it copies, so on
the genome it is `g.20dup`, and no reader can put it in the other exon. An insertion written on a junction, `c.5_6insC`,
is the ambiguity [issue #64](https://github.com/populationgenomics/hgvs-weaver/issues/64) records, and is left as
written. An intronic member's extent is computed on the genome, where it has bases.

A positional identity, `c.6=`, is a barrier: it is merged with no change, and no shift enters it, since it says those
bases are as they were. Identities that overlap or touch are written as one, `c.[7_9=;8=]` as `c.7_9=` and
`c.[7_9=;10=]` as `c.7_10=`, since they state the same bases either way; `c.[10=;20=]` stays two members, and so does
`c.[5=;6=]` on `NM_S.1`, whose two bases are in different exons and would, written `c.5_6=`, name the intron between
them. Nothing reaches across an identity either: a justified region, a group and the test that turns an insertion into a
duplication all stay on one side, so `c.[4A>T;5=;6A>T]` stays three members although c.4 and c.6 share a codon.
`c.[4del;6=]` becomes `c.[5del;6=]`, the deletion stopping at c.5, and `c.[4del;5=]` stays as written. A member that
cannot be applied, a copy number or an uncertain position, is left as written too.

The result is not a canonical form. `c.[5A>G;9G>T]` stays two members, and `c.5_9delinsGACTT`, which makes the same
sequence, stays one delins, because a delins that changes a base at each end
[stays where it is](#a-delins-is-reduced-before-it-is-shifted). Equivalence compares sequences, so neither spelling is a
different change; a normalised string is not a key to compare by. The rejected alternative splits the whole difference
from the reference into blocks of changed bases, which would be canonical. It needs an alignment rule to say which bases
an indel changes, and it would change how a single variant normalises, `c.5_9delinsGACTT` becoming `c.[5A>G;9G>T]`, away
from the end-trimming biocommons and VRS do.

`normalize` writes the members in order of position along the sequence the allele is written on, as the alleles page
lists them "in genomic order", whatever order they were written in. When every member falls into one group, the result
is the plain variant, `c.7_8delinsTG`, not the one-member allele `c.[7_8delinsTG]`. Normalising changes the description,
and the [VRS identifier follows the description](#a-cis-allele-is-the-sequence-its-members-make), so the normalised form
can have a different identifier from the allele written.

A group that makes the reference is written as what it states. `c.[4_5insA;5del]` inserts an A and deletes the A after
it; on `NM_X.1`, whose c.4_6 is `AAA`, the members pass the overlap rule and together make the reference. What remains
is a statement that the bases the members are written over, c.4 and c.5, are unchanged, which is how
[equivalence reads it](#a-cis-allele-is-the-sequence-its-members-make), so `normalize` writes `c.4_5=`, and `to_vrs`,
which normalises, renders the reference Allele over those bases. The statement is about the sequence: c.4 and c.5 read
as the reference does, with nothing between them, which is what the members make. It is not a claim that no member
touched them; read as members of one allele, `c.[4_5=;5del]` would be contradictory, but normalisation replaces the
members with the statement rather than writing it beside them.

A member names the positions written in it: a range every base in it, an insertion its two flanks, a duplication the
bases it copies, and a repeat the [whole run](#a-repeat-resolves-to-its-whole-run) it resolves to. There is no widening
over a run the members sit in: `c.4_6=` would claim c.6, which neither member names. So the statement depends on where
the members were written. On the same run, `c.[5_6insA;6del]` is `c.5_6=`, a different statement from `c.4_5=`, where
the equivalent changes written at different places in a run are one change; that follows from comparing statements by
what they name. Whether members cancel depends on the reference, not the text: on `TTTTTTTTTT`, `g.[4_5insA;5del]` makes
`TTTTATTTTT`, which is `g.5T>A`.

A single variant can make the reference as well, and is written as the same statement. On `NM_X.1`, c.9_12 is `GGCC`,
its own reverse complement, so `c.9_12inv` is `c.9_12=`; `c.7_9delinsCTG` is `c.7_9=`; and `c.4A[3]`, three copies of
the A that c.4_6 already holds, names the run and is `c.4_6=`.

Cancelling members beside a change keep their statement as an identity member, since the statement is part of the allele
([the cis allele entry](#a-cis-allele-is-the-sequence-its-members-make)): `c.[4_5insA;5del;13T>G]` becomes
`c.[4_5=;13T>G]`, and `c.[4_5insA;5del;6A>T]` becomes `c.[4_5=;6A>T]`, although its three members would otherwise meet
in one group whose change, trimmed, is `c.6A>T` alone. That is why statements are found before anything is grouped: a
group's sequence does not record which of its members cancelled. A statement written beside a change is written as one
with any identity it meets, `c.[4_5insA;5del;6=]` becoming `c.4_6=`.

A set of cancelling members is a statement only when its span holds no base another member names. That is what lets the
rest of the allele be written on its own: the members set aside change nothing where any other member acts, so what
remains makes the allele's sequence. Two cases fall outside. A set around another change states nothing, and its members
stay with the change: in `c.[4del;5A>T;6_7insA]` the deletion and the insertion cancel, but c.5, which the substitution
changes, lies between them, and the allele is `c.4A>T`. Its effect is `c.4A>T` whichever way it is read, and nothing it
says is lost. Two sets that share a member are ambiguous: in `c.[3_4insA;5del;6_7insA]` the deletion cancels either
insertion, and which one is not written. Paired with the first, it states c.3 to c.5 and leaves the second insertion,
`c.[3_5=;6dup]`; paired with the second, it states c.5 to c.7 and leaves the first, which stops at the statement,
`c.[4dup;5_7=]`. Both make the allele's sequence, and they state different bases, so the allele is an
`AmbiguousDescriptionError` with the two as its readings rather than `c.6dup`, which would drop a statement the allele
makes either way. Whether a set qualifies depends only on the members, not on the order they are written or looked at,
and the search is confined to members that meet, since members whose extents do not meet change different places and
cannot cancel each other.

The rejected alternatives were to drop a group that makes the reference, which loses the one thing its members still
say, and to refuse an allele that makes the reference, which would also refuse `c.[10=;20=]` and `c.[7_9=;8=]`,
statements with one reading. `validate` passes such an allele when its stated bases are right.

`to_vrs` normalises a cis allele first, so members whose justified regions meet become one Allele. On `GAAAAAG`,
`g.[2del;4A>T]` makes `GATAAG` and `g.[6del;4A>T]` makes `GAATAG`; their deletions both justify over g.2_6, so rendered
member by member, as ga4gh.vrs does, they give one block for two different alleles. Normalised first, they are
`g.3_4delinsT` and `g.4_5delinsT`, two different Alleles. weaver's rendering differs from a block of members translated
one at a time wherever members meet, and that is where such a block can lose this information.

`from_vrs` reads a `CisPhasedBlock` back the same way: each member is trimmed to its minimal edit, every placement of
that edit within the member's justified location is considered, and the placements the overlap rule accepts are
normalised as one allele by this rule. Where they all make one sequence, that is the answer. A block's members are fully
justified one at a time (ga4gh.vrs 2.3.1's translator justifies each Allele, and its block normalisation only sorts them
by digest), and justified locations can overlap where the edits do not: the deletion of `g.[3_5del;7A>T]` on
`TTCAGCAGTT` is justified over g.3_8, which holds g.7. VRS 1.3's `Haplotype` was "a set of non-overlapping Allele
members"; 2.0.1's `CisPhasedBlock` drops the condition. Read back one at a time, the members are either
`g.[6_8del;7A>T]`, an overlap weaver refuses, or overlapping justified spans, which are not a description; grouping is
the only reading that returns an answer. Of the deletion's placements g.3_5, g.4_6, g.5_7 and g.6_8, the first two pass
beside `g.7A>T` and both make `TTCTGTT`, so the block reads back as `g.4_7delinsT`. A block whose members' minimal edits
can be placed in ways that make different sequences, as a ga4gh.vrs block for either `GAAAAAG` allele above can, is
[ambiguous](#an-ambiguous-description-is-an-error-that-lists-its-readings), and reading it back gives both.

**Example** (`NM_X.1` and `NC_X.1` as in [the cis allele entry](#a-cis-allele-is-the-sequence-its-members-make); c.4_12
reads `AAACTGGCC`, codons 2 to 4, and c.4 is g.19; `NM_S.1` is `NM_X.1` with exon 1 ending at c.5, so c.4_9 reads
`AA|ACTG` across the junction):

```text
                                  makes
on TTCAGCAGTT:
X:g.[3_5del;7A>T]                 TTCTGTT      →  X:g.4_7delinsT       one at a time: g.[6_8del;7A>T], an overlap
X:g.[4A>T;6_8del]                 TTCTGTT      →  X:g.4_7delinsT
X:g.4_7delinsT                    TTCTGTT      →  X:g.4_7delinsT
on TCTCT:
X:g.[1del;2del;5del]              TC           →  X:g.3_5del           grouped once: g.[4_5del;5del], an overlap;
                                                                       skipping that shift: g.[1_2del;5del]
on CTCTT:
X:g.[1del;2del;4del]              CT           →  X:g.3_5del           grouped once: g.[3_4del;5del], whose members touch
on AAATAAAA:
X:g.[1_3del;4del;6A>T]            ATAA         →  X:g.2_6delinsT       grouped once: g.[4_7del;6A>T], an overlap
on ACACCCCGC:
X:g.[4del;7_9delinsCGG]           ACACCCGG     →  X:g.7_9delinsGG      by the trimmed delins: g.[7del;7_9delinsCGG]

                                  c.4_12 becomes
NM_X.1:c.[7C>T;8T>G]              AAATGGGCC    →  NM_X.1:c.7_8delinsTG     the members touch
NM_X.1:c.[4del;7C>A]              AAATGGCC     →  NM_X.1:c.7del            one at a time: c.[6del;7C>A]
NM_X.1:c.[7C>A;7_8insC]           AAAACTGGCC   →  NM_X.1:c.6dup            one at a time: c.[7C>A;7dup]
NM_X.1:c.[4A>T;6A>T]              TATCTGGCC    →  NM_X.1:c.4_6delinsTAT    one codon
NM_X.1:c.[4A>T;6A>T;7C>G]         TATGTGGCC    →  NM_X.1:c.4_7delinsTATG   one codon, then touching
NM_X.1:c.[6A>T;8T>G]              AATCGGGCC    →  NM_X.1:c.[6A>T;8T>G]     two codons, c.7 unchanged between
NM_X.1:c.[4A>T;5_6insG]           TAGACTGGCC   →  NM_X.1:c.4_5delinsTAG    the insertion is inside codon 2
NM_X.1:c.[5A>T;6_7insG]           ATAGCTGGCC   →  NM_X.1:c.[5A>T;6_7insG]  the insertion is between codons 2 and 3
NC_X.1:g.[19A>T;21A>T]            TATCTGGCC    →  NC_X.1:g.[19A>T;21A>T]   no codon on g.
NM_X.1:c.[5A>G;9G>T]              AGACTTGCC    →  NM_X.1:c.[5A>G;9G>T]     c.6_8 unchanged between
NM_X.1:c.5_9delinsGACTT           AGACTTGCC    →  NM_X.1:c.5_9delinsGACTT  the same sequence, another string
NM_X.1:c.[10del;12C>T]            AAACTGCT     →  NM_X.1:c.10_12delinsCT   the deletion in codon 4, as written
NM_X.1:c.[9del;12C>T]             AAACTGCT     →  NM_X.1:c.10_12delinsCT   written in codon 3; its region, c.9_10, reaches codon 4

                                  c.4_9 becomes
NM_S.1:c.[5A>T;6A>T]              ATTCTG       →  NM_S.1:c.[5A>T;6A>T]     they touch across the junction
NM_S.1:c.[4A>T;6A>T]              TATCTG       →  NM_S.1:c.[4A>T;6A>T]     codon 2 is split by the intron
NM_S.1:c.[6A>T;7C>G]              AATGTG       →  NM_S.1:c.6_7delinsTG     both in exon 2
NM_S.1:c.[4del;7C>A]              AAATG        →  NM_S.1:c.[5del;7C>A]     the deletion stops at the junction
NM_S.1:c.[4del;6A>T]              ATCTG        →  NM_S.1:c.[5del;6A>T]     the same
NM_S.1:c.4_5insAC                 AACAACTG     →  NM_S.1:c.4_5insAC        not onto the junction, c.5_6insCA
NM_S.1:c.[4dup;6A>T]              AAA|TCTG     →  NM_S.1:c.[5dup;6A>T]     the copy against c.5: g.20dup
NM_S.1:c.[5=;6=]                  AA|ACTG      →  NM_S.1:c.[5=;6=]         identities in two exons stay apart
NM_X.1:c.[4del;6=]                AACTG        →  NM_X.1:c.[5del;6=]       no shift enters c.6
on GAAAAAG:
X:g.[2del;4A>T]                   GATAAG       →  X:g.3_4delinsT           both deletions justify over g.2_6
X:g.[6del;4A>T]                   GAATAG       →  X:g.4_5delinsT

                                  c.*1_*5 becomes
NM_X.1:c.[*1C>A;*3G>T]            ACTTA        →  NM_X.1:c.[*1C>A;*3G>T]   no codon outside the CDS

                                  c.4_12 becomes; c.13 is outside it
NM_X.1:c.[4_5insA;5del]           AAACTGGCC    →  NM_X.1:c.4_5=            the bases the members name
NM_X.1:c.[4_5insA;5del;13T>G]     AAACTGGCC    →  NM_X.1:c.[4_5=;13T>G]    the cancelling group keeps its statement
NM_X.1:c.[4_5insA;5del;6=]        AAACTGGCC    →  NM_X.1:c.4_6=            the two statements touch
NM_X.1:c.[4_5insA;5del;6A>T]      AATCTGGCC    →  NM_X.1:c.[4_5=;6A>T]     the statement is found before grouping
NM_X.1:c.[4del;5A>T;6_7insA]      TAACTGGCC    →  NM_X.1:c.4A>T            the cancelling pair surrounds c.5: no statement
NM_X.1:c.[3_4insA;5del;6_7insA]   AAAACTGGCC   →  ambiguous: c.[3_5=;6dup] | c.[4dup;5_7=]
NM_X.1:c.[9_12inv;13T>G]          AAACTGGCC    →  NM_X.1:c.[9_12=;13T>G]   one member is the statement
NM_X.1:c.[10=;20=]                AAACTGGCC    →  NM_X.1:c.[10=;20=]       two statements, apart
NM_X.1:c.[7_9=;8=]                AAACTGGCC    →  NM_X.1:c.7_9=
NM_X.1:c.[7_9=;10=]               AAACTGGCC    →  NM_X.1:c.7_10=           statements that touch are one
NM_X.1:c.[5_6insA;6del]           AAACTGGCC    →  NM_X.1:c.5_6=            written elsewhere in the run, another statement
NM_X.1:c.9_12inv                  AAACTGGCC    →  NM_X.1:c.9_12=           GGCC is its own reverse complement
NM_X.1:c.7_9delinsCTG             AAACTGGCC    →  NM_X.1:c.7_9=
NM_X.1:c.4A[3]                    AAACTGGCC    →  NM_X.1:c.4_6=            the run it resolves to
on TTTTTTTTTT:
X:g.[4_5insA;5del]                TTTTATTTTT   →  X:g.5T>A                 here the members change a base
```

- **Differs:** a block assembled from members each translated by ga4gh.vrs 2.3.1, whose translator fully justifies an
  Allele by default and whose `normalize.py::_normalize_cis_phased_block` only sorts members by digest, makes
  `g.[2del;4A>T]` and `g.[6del;4A>T]` on `GAAAAAG` one block; weaver merges meeting members first and renders two
  Alleles.
- **Others:** biocommons 1.5.7 does not parse alleles in cis; mavehgvs, VariantValidator and ClinVar not checked.
- **Spec (differs on one point):**
  [deletion-insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/): "changes involving two or more
  consecutive nucleotides are described as deletion/insertion (delins) variants"; "two variants separated by one or more
  nucleotides should be described individually and not as a "delins""; and the exception, "two variants separated by one
  nucleotide, together affecting one amino acid, should be described as a "delins"", which "prevents tools predicting
  the consequences of a variant to make conflicting and incorrect predictions". Grouping by extent until nothing
  changes, the codon rule's reach to members further apart in one codon, and keeping exons apart are weaver's reading.
  The departure: merging can write a delins where the page would describe changes separately, by touch through a
  justified region and by the codon rule, which also merges changes that do not "together affect one amino acid", as a
  frameshift does not (`c.[9del;12C>T]` becomes `c.10_12delinsCT`). `g.[4del;7_9delinsCGG]` on `ACACCCCGC` becomes
  `g.7_9delinsGG`, though `g.[7del;9C>G]` makes the same sequence with an unchanged base between. Merging the members
  that meet is what keeps normalisation sequence-preserving and idempotent, so weaver merges.
- **Tests:** `cis_phased_test::members_normalise_and_validate_one_by_one`, whose two alleles have members that meet no
  other, so its outputs hold under this rule, and `::a_single_variant_that_makes_the_reference_is_its_statement`;
  `a_single_variant_that_makes_the_reference_normalises_to_its_statement` (property). Planned:
  `::members_in_one_run_are_written_as_one_change`, `::touching_members_are_one_delins`,
  `::a_group_can_reduce_to_a_deletion`, `::a_group_can_reduce_to_a_duplication`,
  `::members_with_an_unchanged_base_between_stay_apart`, `::members_in_one_codon_are_one_delins`,
  `::members_one_base_apart_in_two_codons_stay_apart`, `::grouping_runs_across_touch_and_codon`,
  `::a_genomic_allele_has_no_codon`, `::members_outside_the_cds_have_no_codon`,
  `::an_insertion_inside_a_codon_joins_it`, `::a_deletion_across_two_codons_groups_by_its_region`,
  `::members_in_a_codon_split_by_an_intron_stay_apart`, `::an_allele_merged_whole_is_a_plain_variant`,
  `::a_cancelling_allele_is_the_identity_of_the_bases_its_members_name`,
  `::a_cancelling_group_beside_a_change_leaves_an_identity`,
  `::cancelling_members_that_would_meet_a_change_are_still_a_statement`,
  `::cancelling_members_around_a_change_are_part_of_it`, `::cancelling_sets_that_share_a_member_are_ambiguous`,
  `::a_single_member_that_makes_the_reference_beside_a_change_is_a_statement`,
  `::a_change_is_cut_at_an_agreeing_identity`, `::identities_that_meet_are_written_as_one`,
  `::identities_apart_stay_apart`, `::identities_in_two_exons_stay_apart`,
  `::a_duplication_in_an_allele_stops_at_its_exons_last_base`; `normalising_a_cis_allele_keeps_its_sequence`,
  `normalising_a_cis_allele_twice_is_normalising_it_once`, `normalising_a_cis_allele_keeps_its_stated_bases` and
  `a_cis_alleles_statements_do_not_depend_on_member_order` (properties).

## Protein consequences

Protein consequences are read from codons, not from a diff of two protein strings. The rules that follow all bear on
where the stop is. The unit-test examples below use the coding sequence `ATG AAA CTG GCC TAT CGC TAA` (Met Lys Leu Ala
Tyr Arg Ter) followed by a 3'UTR, and their `c.` spellings are read from the 0-based ranges the tests assert.

### The declared CDS end is the stop

**The reference protein ends where the transcript record declares the CDS ends, when that codon really is a stop; only a
loosely marked CDS end falls back to the first stop in the translation.**

A selenocysteine `TGA` inside the CDS is therefore not a stop.

**Example** (SEPN1 has selenocysteine codons at 127 and 462; a toy CDS `ATG AAA TGA CTG TAA` reads Met Lys Sec Leu Ter):

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

- **Agrees:** biocommons and ClinVar.
- **Spec (agrees):** [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/): "the variant
  extends the amino acid sequence at the C-terminal end and is therefore by definition an extension".
- **Tests:** `src/protein.rs::a_frameshift_in_the_stop_codon_is_an_extension`,
  `::stop_loss_is_an_extension_to_the_next_stop`; `real_transcripts_test::biocommons_real_transcript_cases` (EXT01 to
  EXT06).

### An extension needs the stop codon itself to change

**An in-frame change just before the stop that leaves the stop intact is an insertion or deletion, even though the
stop's position is the first residue that differs.**

**Example** (CDS `ATG CTG TTT GTA TTG TGT CGT CTT TAA`, Met Leu Phe Val Leu Cys Arg Leu Ter):

```text
c.21_22insTTGTCT   →  p.Leu8_Ter9insSerLeu    the stop is intact two codons on
toy c.16_18del     →  p.Arg6del               the stop shifts left but is the same stop
```

- **Others:** not checked.
- **Spec (agrees):** [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/), the same
  definition: the stop codon is what changes.
- **Tests:** `src/protein.rs::an_insertion_before_the_stop_that_repeats_the_last_residue_is_not_an_extension`,
  `::an_in_frame_deletion_reaching_the_stop_keeps_it_original`.

### A stop formed inside inserted bases is a delins ending in Ter

**Translation ends before the new frame reads a single reference base, so the change is a delins ending in Ter, not a
frameshift.**

**Example** (the toy CDS; and a two-codon CDS `AAA GGG`):

```text
toy c.4_12delinsCTGTAA     →  p.Lys2_Ala4delinsLeuTer    in frame
toy c.4_12delinsCTGTAAGG   →  p.Lys2_Ala4delinsLeuTer    9 bases out, 8 in: the frame breaks, but after the stop
NM_1.1:c.1_2delinsTA       →  NP_1.1:p.(Lys1Ter)         not p.Lys1fsTer1
```

- **Others:** not checked.
- **Spec (agrees):** [frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/): "variants
  which introduce an immediate translation termination (stop) codon are described as nonsense variant", not a
  frameshift.
- **Tests:** `src/protein.rs::a_stop_inside_the_inserted_bases_is_a_delins_not_a_frameshift`;
  `equivalence_cases_test::test_immediate_stop_normalization`.

### In-frame changes are written 3'-most

**The shared tail is trimmed, so an in-frame deletion or insertion is written at its 3'-most equivalent residues.**

**Example** (CDS `ATG AAA AAA CTG TAA`, Met Lys Lys Leu Ter):

```text
c.4_6del   →  p.Lys3del    the first Lys codon is deleted; the second is named
c.4_6dup   →  p.Lys3dup
```

- **Differs:** ClinVar, in spelling: it writes `Xxx_Yyyins…` at the 3' end of a run where weaver writes `dup`, and
  one-letter repeat forms such as `p.490PRS[1]`; 170 rows in 100,000 are the same protein in another spelling, and
  equivalence judges them Analogous.
- **Spec (agrees):** [protein deletion](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/): "the
  most C-terminal position possible of the reference sequence is arbitrarily assigned to have been changed".
- **Tests:** `src/protein.rs::in_frame_deletion_and_duplication_take_the_3_prime_position`;
  `protein_insertions_are_written_3_prime_most` and `in_frame_indels_describe_the_translated_protein` (properties).

### A start-codon change is written specifically

**`c.1A>G` is `p.(Met1Val)`, not `p.(Met1?)`, because the specific prediction carries more information.**

This is the one choice on this page where weaver's default output departs from an explicit recommendation, and it is a
contentious one. The nomenclature's reasoning is that translation may initiate elsewhere, so the consequence cannot be
predicted; weaver's is that `p.(Met1Val)` says which codon changed and how, is written as a prediction, and can be
reduced to `p.(Met1?)` but not recovered from it. `VariantTransformSettings` with `StartCodonConvention::HgvsQuestion`
rewrites it on request, keeping whether it was a prediction: `p.(Met1Val)` becomes `p.(Met1?)`, a stated `p.Met1Val`
becomes `p.Met1?`.

`p.Met1?` read from text is the same statement weaver writes: an unknown change from `Met1` on, not a substitution of
`Met1` by an unknown residue, which would be `p.Met1Xaa`, the commitment the nomenclature declines to make. The
bracketed `p.(Met1?)` the nomenclature writes is the same statement as a prediction, and is what `c_to_p` writes
([below](#a-prediction-is-bracketed-as-hgvs-writes-it)). So what `c_to_p` and the transform write reads back as the
value they built. The statement names one residue, certain: `p.Met1_Arg3?` and `p.(Met1)?` do not parse.

**Example:**

```text
NM_0001.3:c.1A>T          →  NP_0001.1:p.(Met1Leu)
NM_007199.2:c.1A>G        →  NP_009130.2:p.(Met1Val)     weaver
                          →  NP_009130.2:p.Met1?         biocommons and VariantValidator

transform NP_000051.2:p.(Met1Val)  →  NP_000051.2:p.(Met1?)    with start_codon = HgvsQuestion
transform NP_000051.2:p.Met1Val    →  NP_000051.2:p.Met1?

parse NP_W.1:p.(Met1?)    →  the unknown change from Met1 that c_to_p builds for NM_W.1:c.-3_4del
parse NP_W.1:p.Met1?      →  the same, stated
```

- **Differs:** biocommons and VariantValidator write `p.Met1?`, the nomenclature's `p.(Met1?)` without its parentheses.
  A recorded difference. biocommons reads the `?` of `p.Met1?` as a substitution to an unknown residue (its grammar
  table); weaver reads it as the statement.
- **Spec (differs):**
  [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/): `p.(Met1?)` when
  "the consequence, on the protein level, of a variant affecting the translation initiation codon can not be predicted".
- **Tests:** `mapping_test::test_mapper_c_to_p_start_codon_subst`;
  `real_transcripts_test::biocommons_real_transcript_cases` (INITMET01, recorded);
  `src/transform.rs::test_transform_met1_to_question`, `::test_transform_stated_met1_to_stated_question`;
  `statements_test::a_predicted_p_met1_unknown_parses_to_what_the_mapper_builds`,
  `::the_start_codon_transforms_unknown_change_from_met1_reparses_to_itself`,
  `::the_start_codon_transform_of_a_one_letter_description_reparses_to_itself`,
  `::a_predicted_unknown_change_from_met1_prints_back_as_parsed`.

### Edits outside the CDS are statements

**A deletion of the whole CDS is `p.0?`; an edit that starts in the 5'UTR and reaches into the CDS is `p.(Met1?)`; an
edit entirely upstream is `p.?`; an edit entirely in the 3'UTR is `p.(=)`; an edit that leaves no codon where the
protein first changes is `p.?`.**

weaver does not commit to a consequence upstream because the initiation site is not predictable: a 5'UTR change can
create an upstream start. Downstream of the stop codon nothing can change the protein, so the prediction is no change,
and no residue past the end of the protein is named. An edit that runs to the end of the transcript can leave fewer than
three bases after the last intact codon: translation runs off the edited transcript with no stop and no residue to name
at the first change, so weaver says nothing rather than read the missing codon as a stop.

These statements are about a change. An identity, `=`, changes nothing wherever it lies, so outside the codons it
predicts `p.(=)` too: one in the 5'UTR, in an intron, or covering every codon (`c.-5_*15=`, every base of the
transcript) leaves the protein as it is, not disrupted from `Met1` and not unknown. An identity inside the codons is
still written silently at the residues it covers, `p.(Arg2=)`, and so is one that reaches into them from a UTR without
covering all of them: `c.-2_3=` covers codon 1 and is `p.(Met1=)`, a statement about one residue, where `p.(=)` would
claim the whole protein (planned, with
[the comparison of silent statements](#judged-by-allele-and-by-the-protein-left-behind)). A conversion is not counted as
an identity here: intronic or upstream it stays `p.?`. weaver does not model a conversion's bases (it parses to an edit
with none, printed `=`), so inside the codons one is still read as silent; that is an open defect, not a choice.

**Example** (`NM_X.1` has a five-base 5'UTR and a fourteen-base 3'UTR):

```text
NM_X.1:c.-3G>A             →  NP_X.1:p.?         entirely upstream
NM_X.1:c.-3_2del           →  NP_X.1:p.(Met1?)   reaches into the start codon
NM_X.1:c.-5_*14del         →  NP_X.1:p.0?        the whole CDS
NM_X.1:c.*3G>A             →  NP_X.1:p.(=)       entirely downstream
NM_X.1:c.21A>G             →  NP_X.1:p.(Ter7=)   the stop codon itself, silently
NM_X.1:c.21_*1dup          →  NP_X.1:p.(Ter7=)   across the stop, named up to it
NM_X.1:c.10_*14del         →  NP_X.1:p.?         residue 4 to the end of the transcript
NM_W.1:c.-5_*15=           →  NP_W.1:p.(=)       an identity of the whole transcript
NM_W.1:c.45+1=             →  NP_W.1:p.(=)       an intronic identity

NM_000249.3:c.-7_*46del    →  NP_000240.1:p.0?
NM_022051.2:c.-1_1insGCC   →  NP_071334.1:p.(Met1?)
```

- **Agrees:** biocommons, on `p.0?` and on the start codon, which it writes `p.Met1?`, unbracketed (its own test table).
- **Agrees:** VariantValidator, on `p.(=)` for the 3'UTR (older versions wrote `p.?`); ClinVar writes a 3'UTR change as
  silent at the stop, `p.Ter1648=`, which leaves the same protein; once silent statements are compared by the residues
  they name (planned), that spelling is Different from weaver's `p.(=)`, and the ClinVar gate rows that carry it will
  move.
- **Differs:** ClinVar writes `Met1fs`, `Met1_Glu2insGly…`, committing to a consequence. For a substitution entirely in
  the 5'UTR VariantValidator writes `p.(=)` where weaver writes `p.?`; both are defensible.
- **Spec (agrees; silent on an edit that leaves no codon where the protein changes):**
  [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/): `p.0?` "when you
  predict that no protein is produced"; `p.Met1?` as above; `p.Met1ext-5` shows a 5'UTR change can activate an upstream
  initiation site.
- *Not checked:* biocommons, VariantValidator and ClinVar on an edit that leaves no codon where the protein changes.
- **Tests:** `protein_allele_from_coding_test::edits_around_the_cds_start_are_statements`,
  `::a_change_entirely_in_the_three_prime_utr_leaves_the_protein_unchanged`,
  `::a_silent_change_across_the_stop_names_no_residue_past_it`,
  `::an_edit_that_leaves_no_codon_where_the_protein_changes_says_nothing`,
  `::statements_have_no_allele_and_a_wrong_protein_is_an_error`,
  `::a_coding_variant_agrees_with_every_spelling_of_its_consequence`;
  `open_cds_test::a_3_open_cds_deleting_residues_that_repeat_its_last_is_a_deletion`,
  `::a_cds_with_no_stop_shifted_to_the_end_of_its_record_says_nothing`;
  `real_transcripts_test::biocommons_real_transcript_cases` (WHOLEGENE01/02, INITMET02/03);
  `statements_test::an_identity_across_the_whole_transcript_predicts_the_protein_unchanged`,
  `::an_identity_from_the_5_utr_into_the_cds_predicts_the_protein_unchanged` (to change: `p.(Met1=)` for `c.-2_3=`,
  planned, in the PR that implements the silent statements, which reports each ClinVar gate row it moves),
  `::an_identity_in_the_5_utr_predicts_the_protein_unchanged`, `::an_intronic_identity_predicts_the_protein_unchanged`,
  `::an_identity_across_the_whole_transcript_is_not_no_protein`,
  `::an_identity_inside_the_codons_is_written_silently_at_its_residues`,
  `::an_intronic_conversion_is_an_unknown_effect`, `::an_upstream_conversion_is_an_unknown_effect`.

### A prediction is bracketed as HGVS writes it

**A protein consequence predicted from a coding or RNA change is written in parentheses, `p.(Arg2Gly)`, as are the
predicted statements `p.(=)` and `p.(Met1?)`; `p.?` and `p.0?` are written bare.**

The nomenclature marks a consequence that was predicted rather than observed by putting it in parentheses:
"`p.(Trp24Cys)` amino acid `Trp24` is predicted to change to a `Cys` (no experimental proof, e.g., based on DNA level
data)". Every protein description weaver writes from a `c.` or `r.` change is such a prediction, since weaver has seen
no protein. For a change the rule needs no exception. The statements are where the nomenclature's own spellings vary,
and weaver follows the spelling it gives each one rather than bracketing all of them:

- `p.(=)`: "no changes on the protein level are expected (RNA not analysed)", bracketed.
- `p.(Met1?)`: the start codon disrupted, its consequence unknown, bracketed on the substitution page.
- `p.?`: written bare even where it is the predicted consequence ("the predicted consequences of a deletion of exons 1
  to 2 … `p.?`"). What it says is that no prediction can be made, so there is nothing to bracket.
- `p.0?`: "can be used when you predict that no protein is produced"; the `?` is the hedge, and `p.(0?)` appears
  nowhere. It is what weaver writes for a deleted CDS and for `r.0` and `r.0?`: no RNA predicts no protein. The bare
  `p.0` is an observation, as in the alleles page's `p.[Ser86Arg];[0]`, linked there to "the absence of RNA (`r.0`)".
- `p.(?)`, which weaver parses but never writes: the alleles page writes it for "a variant expected but not yet
  identified", in prose beside the allele `p.[(Ser68Arg)];[?]`, so it is read as `p.?` predicted.

A statement the nomenclature has no example for would take the parentheses, as a change does; there is none today. A
parsed description keeps the parentheses it was written with, and the start-codon transform keeps whether its input was
a prediction. A prediction's parentheses enclose the whole description; those of an uncertain residue or range,
`p.(Gln18)[(70_80)]`, belong to the position, so a description that opens with one is not a prediction, and
`p.((Gln18)[(70_80)])` is. In equivalence parentheses are spelling between two `p.` descriptions: `p.(Met1?)` and
`p.Met1?` are Identity, as `p.(?)` and `p.?` are. Against a coding change the prediction written exactly is Identity and
the same description in another spelling Analogous, so `c.-3_4del` is Identity with `p.(Met1?)` and Analogous with
`p.Met1?`, as `c.10C>G` is with `p.(Arg4Gly)` and `p.Arg4Gly`; Identity is therefore not transitive across the two. The
parentheses of `r.(?)` are not a prediction but part of the statement
([above](#mitochondrial-and-rna-descriptions-share-the-dna-machinery)).

**Example:**

```text
NM_W.1:c.10C>G        →  NP_W.1:p.(Arg4Gly)
NM_W.1:c.-3_4del      →  NP_W.1:p.(Met1?)       biocommons writes p.Met1?
NM_W.1:c.*5T>A        →  NP_W.1:p.(=)
NM_W.1:c.-3C>A        →  NP_W.1:p.?
NM_W.1:c.-3_*2del     →  NP_W.1:p.0?
NM_W.1:r.0            →  NP_W.1:p.0?
```

- **Differs:** biocommons writes the start codon `p.Met1?` (its test table, INITMET02/03). VariantValidator: not
  checked. The ClinVar sample has no start-codon statement.
- **Spec (agrees):** [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/)
  for `p.(Trp24Cys)`, `p.(Met1?)` and `p.0?`;
  [uncertain](https://hgvs-nomenclature.org/stable/recommendations/uncertain/) for `p.(=)`, `p.?` and `p.0?`;
  [protein deletion](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/) for a predicted `p.?`;
  [protein alleles](https://hgvs-nomenclature.org/stable/recommendations/protein/alleles/) for `p.0` observed beside
  `r.0` and for `p.(?)`.
- **Tests:** `statements_test::a_predicted_change_is_written_bracketed`,
  `::a_start_loss_is_predicted_bracketed_as_hgvs_writes_it`,
  `::a_change_upstream_of_the_cds_is_an_unknown_effect_written_bare`,
  `::predicted_p_unknown_is_the_unknown_effect_in_another_spelling`,
  `::predicted_p_unknown_does_not_match_a_coding_change`,
  `::a_3_utr_change_predicts_the_protein_unchanged_and_reparses_to_it`,
  `::a_start_loss_is_the_parsed_predicted_unknown_change_from_met1`,
  `::a_start_loss_is_analogous_to_the_unknown_change_from_met1_stated`;
  `protein_allele_from_coding_test::edits_around_the_cds_start_are_statements` (`p.?`, `p.(Met1?)`, `p.0?`);
  `rna_test::r_predicts_the_protein_like_c` (`r.0` to `p.0?`); `real_transcripts_test::biocommons_real_transcript_cases`
  (INITMET02/03); `src/transform.rs::test_transform_met1_to_question`,
  `::test_transform_stated_met1_to_stated_question`;
  `printed_reparse_test::a_change_at_an_uncertain_residue_is_not_a_prediction`,
  `::a_change_over_an_uncertain_residue_range_is_not_a_prediction`,
  `::a_predicted_change_at_an_uncertain_residue_keeps_both_brackets`.

### Frameshift length counts to the first new stop

**`fsTer N` counts residues from the first changed one to the first stop the new frame reaches; an immediate stop is a
plain nonsense substitution.**

**Example** (the toy CDS):

```text
toy c.17_18del   →  p.Arg6LeufsTer4    Leu Ala Val Ter: four residues counting the stop
toy c.4A>T       →  p.Lys2Ter          AAA → TAA, immediate
toy c.4A>C       →  p.Lys2Gln
toy c.6A>G       →  p.Lys2=            AAA → AAG
```

- **Others:** not checked.
- **Spec (agrees):** [frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/): the stop
  position is counted "starting from the first changed amino acid as codon 1"; "the shortest frameshift variant possible
  contains `fsTer2`".
- **Tests:** `src/protein.rs::a_frameshift_in_the_stop_codon_is_an_extension`,
  `::frameshift_reports_first_changed_residue_and_distance_to_stop`, `::substitutions_nonsense_and_silence`;
  `snv_effect_follows_the_codon_table` (property).

### Protein alleles come from the coding change

**From a `p.` description, only edits that name a sequence have an allele. From a coding variant, every consequence has
one, because the edited transcript's translation is known.**

A frameshift's allele is the residues from the first change to the new stop. The translated CDS must be the protein the
provider serves; a difference is an annotation error and is reported as one. A silent change is the reference allele
over the codons the edit touched; a deleted CDS is the deletion of the whole protein.

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

- **Spec (agrees; silent on deriving a protein allele from a coding change):** VRS 2.0.1 defines an `Allele` on a
  protein `SequenceReference`; HGVS has no allele.
- **Tests:** `protein_allele_test::consequences_have_no_allele`;
  `protein_allele_from_coding_test::consequences_without_a_p_sequence_still_have_an_allele`,
  `::definite_changes_give_the_same_allele_by_either_route`,
  `::statements_have_no_allele_and_a_wrong_protein_is_an_error`; `the_protein_allele_is_the_edited_translation`
  (property).

## Validation and stated bases

### A stated base the reference lacks is an error

**The bases a description states — a substitution's, those written after `del`, `dup` or `inv`, the deleted part of a
`delins`, an identity's `c.5A=`, a repeat's unit — and the residues a `p.` description names are compared with the
reference they cite by every operation that reads the description against a sequence. Where they differ the answer is a
`ReferenceMismatchError` naming the stated and the reference's bases; `validate` is the same comparison answered as true
or false.**

A stated base is a claim about the reference. When it is not there, the description was written against another version,
on the other strand, at another position, or with a typo, and no reading of it is safe: using the reference's base
answers a question the author did not ask, and using the stated one describes a sequence that does not exist. A
projection that quietly takes the reference's base is the worst of these: `c.316+5A>G` over a genomic G becomes `g.…=`,
which reads as "no change" rather than "wrong reference".

The comparison is one function on the edit, and resolving an edit to its bases, which projection, normalisation, the
canonical allele, SPDI, VRS and the protein consequence all do, runs it first. The paths that carry an edit without
resolving it call it too: an edit projected as given, an edit normalisation does not move, a `delins` of unknown bases,
a protein consequence that is a statement (`p.(=)` for a 3'UTR change, `p.?` for an intronic one). Residues have their
own comparison, through which every operation on a `p.` description locates it, so `p_to_c` needs the protein's sequence
as well as the transcript's; an SPDI's deletion has its own too. Two operations read no reference and check nothing:
converting between `c.` and `r.` letters, and judging two descriptions equal because their text is the same.

What a base is compared with is the sequence the position numbers:

- an exonic `c.`, `n.` or `r.` base, the transcript record: a `c.` position numbers the record, and where record and
  genome differ the record's base is the one to state. Projected, the edit is then written against the genome's own
  bases ([A projection states the bases of the target](#a-projection-states-the-bases-of-the-target)); that the two
  sequences differ is a property of the reference pair, not of the variant;
- an intronic one, the genome, read on the transcript's strand: the record has no base there. A range with one end in an
  intron is read on the genome throughout, the only sequence holding every base it names;
- a `g.` or `m.` base, the genome; a `p.` residue, at either end of the interval or stated by the edit, the protein,
  where `Ter` one past a sequence served without its stop is the stop.

A count (`del3`, `delins` with no bases) states nothing, and an edit that states nothing is not read for this.

**Example** (`NM_S.1` is spliced from `NC_S.1`, whose g.15 is c.9+1, a G; c.4 is A; the record has C at c.12 where the
genome has T; `NM_R.1` is the same record on the minus strand of `NC_R.1`; `NP_S.1` is `MKLAYR`):

```text
validate NM_S.1:c.9+1G>A        →  true
validate NM_S.1:c.9+1C>A        →  false                         intronic, read on the genome
c_to_g   NM_S.1:c.9+1C>A        →  ReferenceMismatchError       expected C, found G on NC_S.1
c_to_g   NM_R.1:c.9+1G>A        →  NC_R.1:g.36C>T                the transcript's strand
c_to_g   NM_S.1:c.12T>A         →  ReferenceMismatchError       the genome's T is not the record's C
c_to_g   NM_S.1:c.12C>A         →  NC_S.1:g.27T>A                the record's C, written for the genome
c_to_g   NM_S.1:c.4_6CTG[5]     →  ReferenceMismatchError       the unit is not there
c_to_p   NM_S.1:c.*3A>T         →  ReferenceMismatchError       though the prediction reads no codon
g_to_c   NC_S.1:g.15C>A         →  ReferenceMismatchError
normalize, to_spdi, to_vrs, canonical_allele, equivalent  NC_S.1:g.6C>T  →  ReferenceMismatchError
to_vrs   NC_S.1:g.6_7delCCinsN[3] →  ReferenceMismatchError
from_spdi NC_S.1:5:C:T          →  ReferenceMismatchError       an SPDI deletion states bases too
p_to_c, to_vrs, normalize  NP_S.1:p.Arg2Leu  →  ReferenceMismatchError   residue 2 is Lys
```

- **Agrees:** biocommons (1.5.7) and VariantValidator refuse an exonic mismatch, `NM_000059.4:c.316A>G`, "Variant
  reference (A) does not agree with reference sequence (G)", and read MUC2's `NM_002457.5:c.12468G>A` against the
  record's C; biocommons checks the residues a `p.` interval names.
- **Differs (intronic):** biocommons does not check an intronic base: its validator refuses `NM_000059.4:c.316+5C>T` and
  `c.316+5G>T` alike ("Cannot validate sequence of an intronic variant"), and its `c_to_g` rewrites both to
  `NC_000013.11:g.32319330G>T` and `c.316+5A>G` to `g.32319330=`. VariantValidator reads the base on the genome as
  weaver does, but warns and corrects rather than refusing:
  `IntronBaseMismatchWarning … Variant reference (C) does not agree with reference sequence (G)`, then answers for
  `c.316+5G>T`.
- *Not checked:* ClinVar, whose descriptions do not test the choice: none of the 100,000 sampled states a base its
  reference lacks.
- **Spec (agrees):** [general](https://hgvs-nomenclature.org/stable/recommendations/general/): a description is of the
  reference sequence it cites; VRS 2.0.1 states an `Allele` over the sequence's own bases, which a mismatched
  description cannot reach.
- **Tests:** `stated_bases_test::validate_refuses_an_intronic_base_the_genome_does_not_carry`,
  `::validate_refuses_the_genomes_base_where_the_record_differs`,
  `::validate_reads_a_range_across_an_exon_boundary_on_the_genome`,
  `::c_to_g_refuses_an_intronic_base_the_genome_does_not_carry`,
  `::c_to_g_carries_an_intronic_base_on_the_transcripts_strand`,
  `::c_to_g_refuses_an_exonic_base_the_record_does_not_carry`,
  `::c_to_g_writes_the_genomes_base_where_the_record_differs`,
  `::c_to_g_refuses_a_repeat_unit_the_record_does_not_carry`,
  `::validate_refuses_an_intronic_repeat_unit_the_genome_does_not_carry`,
  `::c_to_p_refuses_a_three_prime_utr_base_it_predicts_no_change_from`,
  `::g_to_c_refuses_an_intronic_base_the_genome_does_not_carry`,
  `::normalize_refuses_a_substitution_whose_base_is_not_there`, `::spdi_refuses_a_base_the_genome_does_not_carry`,
  `::to_vrs_refuses_a_base_the_genome_does_not_carry`, `::canonical_allele_refuses_a_base_the_genome_does_not_carry`,
  `::to_vrs_refuses_deleted_bases_of_a_length_insertion_the_genome_does_not_carry`,
  `::from_spdi_refuses_a_deletion_the_sequence_does_not_carry`, `::p_to_c_refuses_a_residue_the_protein_does_not_carry`,
  `::normalize_refuses_a_residue_the_protein_does_not_carry`, and the file's other tests, one per operation;
  `transcript_coordinates_test::validate_checks_stated_reference_through_transcript_coordinates`;
  `protein_allele_test::validation_checks_the_named_and_stated_residues`;
  `projection_reference_real_test::real_records_that_differ_from_the_genome_project_to_the_targets_bases` (the
  `validate` cases); `tests/test_validate_stated_bases.py::test_a_base_the_genome_lacks_is_an_error_naming_the_genome`.

## Equivalence

### Judged by allele and by the protein left behind

**Two nucleotide descriptions are the same change when they make the same sequence on the genome and state the same
bases unchanged; for two single variants that change the sequence, that is when their canonical alleles are equal. A
coding variant and a protein description agree when they leave the same protein.**

A copy number written alone has no allele, and is compared by its projection to the genome, as it is inside a cis
allele: `NM_X.1:c.10_12copy3` and `NC_X.1:g.25_27copy3` are the same, and `copy3` and `copy4` are Different.

`Identity` is the same text after spelling normalisation, or a prediction written exactly as given; the parentheses of a
prediction and the brackets of a one-member allele are spelling. Cis alleles are compared by
[the sequence their members make](#a-cis-allele-is-the-sequence-its-members-make).
[How it decides](equivalence_logic.md) has the full rules.

A protein description that leaves the protein unchanged is a statement in the same way, and is the same only as a
statement about the same residues. `p.(Leu3=)` and `p.(Ala4=)` leave the same protein but say it of different residues,
so they are Different, as `c.10=` and `c.20=` are; `p.(=)`, about the whole protein, is Different from each. A coding
change is compared as the statement it predicts: `c.9G>A`, silent in codon 3, predicts `p.(Leu3=)`, while `c.*3G>A`, in
the 3'UTR, predicts `p.(=)` ([edits outside the CDS](#edits-outside-the-cds-are-statements)), so the first is the same
as `p.(Leu3=)` and the second as `p.(=)`. A coding allele that changes the sequence and states bases predicts both,
`c.[13T>G;10_12=]` giving `p.[(Ala4=;Tyr5Asp)]`, and is compared as that. The rows marked *planned* are those where the
comparison by protein left behind answers otherwise today.

**Example:**

```text
NM_X.1:c.10_12copy3          vs  NC_X.1:g.25_27copy3                       Analogous: the same stretch on the genome (planned)
NP_X.1:p.(Leu3=)             vs  NP_X.1:p.(Ala4=)                          Different: Leu3 and Ala4 (planned)
NP_X.1:p.(Leu3=)             vs  NP_X.1:p.(=)                              Different: one residue, the whole protein (planned)
NM_X.1:c.9G>A                vs  NP_X.1:p.(Ala4=)                          Different: c.9G>A predicts p.(Leu3=) (planned)
NM_X.1:c.9G>A                vs  NP_X.1:p.(Leu3=)                          Identity: the prediction as weaver writes it
NM_X.1:c.*3G>A               vs  NP_X.1:p.(=)                              Identity
NM_X.1:c.[13T>G;10_12=]      vs  NP_X.1:p.(Tyr5Asp)                        Different: it predicts p.[(Ala4=;Tyr5Asp)] (planned)
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
  `equivalence_cases_test::test_clinvar_regression_tyr165ter`, `::test_analogous_fs_wildcard_unification`;
  `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`;
  `equivalence_test::a_lone_copy_number_compares_by_projection`; `a_variant_agrees_with_its_projection`,
  `equivalence_of_single_variants_is_transitive`, `a_single_statement_never_equals_a_change` and
  `single_statements_are_the_same_when_they_name_the_same_bases` (properties). Planned:
  `equivalence_test::silent_statements_at_different_residues_differ`,
  `::a_silent_residue_is_not_the_whole_protein_unchanged`, `::a_coding_change_predicted_silent_is_its_statement`. These
  rules move rows of the ClinVar gate, ClinVar's `p.Ter=` spellings among them; the PR that implements them reports each
  moved row and checks it against the truth columns.

### A description that says nothing matches nothing

**`p.?` matches no other description; `p.Met1?` says only where the change starts and matches only a description whose
change starts there.**

It does not match ClinVar's `Met1fs`, a commitment weaver did not make. Two ClinVar rows judge Different for this
reason, on purpose. `p.Met1?` and `p.Met1Xaa`, an unknown residue in place of Met1, are Analogous, not the same text.

A `?` after any other residue, `p.Arg12?` or `p.(Arg12?)`, is a recorded extension beyond the specification, which
writes the form only at the start codon. weaver reads it as `p.Met1?` is read, a change from that residue on of which
nothing more is said, and it parses, prints, reloads from a dict and compares as that: Analogous to a change starting at
`Arg12`, from a `p.` description or predicted from a coding change, and Different from one starting anywhere else. The
alternatives were to refuse it, which loses a description others write, or to read it as a substitution to an unknown
residue, `p.Arg12Xaa`, the stronger claim that only that one residue changed; the weaker reading claims less than the
text might mean and never more.

**Example:**

```text
NM_X.1:c.-3_2del     vs  NP_X.1:p.Leu3fs    not equivalent: the deletion is p.(Met1?), which promises nothing about Leu3
NM_X.1:c.-5_*14del   vs  NP_X.1:p.0?        equivalent
NM_X.1:c.-5_*14del   vs  NP_X.1:p.0         equivalent
NP_W.1:p.Arg12?      vs  NP_W.1:p.Arg12Gly  Analogous: the change starts at Arg12
NM_W.1:c.34C>G       vs  NP_W.1:p.Arg12?    Analogous: predicted p.(Arg12Gly)
NP_W.1:p.Arg12?      vs  NP_W.1:p.Arg13Gly  Different
```

- **Differs:** ClinVar commits.
- **Spec (agrees; silent on a `?` after any residue but `Met1`):**
  [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/) defines `p.?` and
  `p.(Met1?)` as statements of what is not known. The reading of `p.Arg12?` is weaver's.
- **Tests:** `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`;
  `statements_test::a_start_loss_is_the_parsed_predicted_unknown_change_from_met1`,
  `::a_start_loss_is_analogous_to_the_unknown_change_from_met1_stated`,
  `::an_unknown_change_from_a_residue_is_an_unknown_effect`, `::an_unknown_change_over_a_range_does_not_parse`,
  `::an_unknown_change_from_an_uncertain_residue_does_not_parse`,
  `::an_unknown_change_from_met1_is_not_the_same_text_as_an_unknown_residue_there`,
  `::an_unknown_change_from_a_residue_prints_back_as_parsed`, `::an_unknown_change_from_a_residue_survives_json`,
  `::an_unknown_change_from_a_residue_matches_a_change_starting_there`,
  `::an_unknown_change_from_a_residue_does_not_match_a_change_starting_after_it`,
  `::an_unknown_change_from_a_residue_matches_a_coding_change_predicted_there`,
  `::a_predicted_unknown_change_from_a_residue_is_the_stated_one`;
  `printed_reparse_test::an_unknown_change_from_a_residue_past_met1`.

### Judging with no protein sequence is an error

**Two protein descriptions are compared on the protein they leave, which needs the sequence; without it the judgement is
an error, not a guess.**

Earlier versions reconciled the residues two descriptions happened to name.

**Example:**

```text
NP_0001.1:p.Ala201_Val202insGlyProGlyAla  vs  NP_0001.1:p.Gly198_Ala201dup
    → Err, not Unknown: residues 199 and 200 are named by neither, and may well make these one change
```

- **Spec (silent):** none.
- **Tests:** `equivalence_test::equivalence_needs_the_sequence`.

### Versions of one protein accession compare on ours

**ClinVar often names an older `NP_` version; both descriptions are read against the protein the provider serves. A
different accession is not the same protein.**

**Example:**

```text
NP_X.1:p.(Lys2Gln)   vs  NP_X.0:p.Lys2Gln   equivalent
NP_X.1:p.Lys2Gln     vs  NP_Y.1:p.Lys2Gln   not equivalent
```

- **Spec (silent):** none.
- **Tests:** `protein_allele_from_coding_test::a_coding_variant_agrees_with_every_spelling_of_its_consequence`.

### A cis allele is the sequence its members make

**`c.[a;b]` names the sequence that applying every member to the reference makes. Two descriptions are the same change
when they make the same sequence, whatever order the members are written in, however the change is split between them,
and whether or not it is written as an allele at all.**

HGVS defines an allele as "a series of variants on one chromosome" and writes each member as it would be written alone:
`c.[2376G>C;3103del]` contains "two different changes, `c.2376G>C` and `c.3103del`". Each member is therefore numbered
against the reference, not against a sequence an earlier member has already edited, and applying the members in any
order gives one sequence. "Variants should be listed in genomic order", from the alleles page's answer on a simplified
haplotype format, is a convention for writing them, not part of what they mean.

The specification goes further and treats different splits of one change as spellings of one allele. Its delins page
recommends `c.235_237delinsTAT` over `c.[235A>T;237G>T]`, because a tool reading the members one at a time would predict
"two different substitutions at one position", `p.[Lys79*;Lys79Asn]`, where the codon in fact becomes Tyr. It gives
`c.[850_869del;874_881del;887_897del;901_902insG]` as "an alternative description" of `c.850_901delinsTTCCTCGATGCCTG`.
And asked whether a BRCA1 change from `..CATGACA..` to `..CATATAACA..` may be written `c.[2077G>A;2077_2078insTA]`, it
answers that "the correct description of this variant is `NM_007294.3:c.2077delinsATA`": a correction of the spelling,
for the same sequence. An earlier answer added that the allele form was acceptable; the page has withdrawn that, so
reading `c.[2077G>A;2077_2078insTA]` as an allele with one reading is weaver's choice
([below](#overlapping-members-of-a-cis-allele-need-one-reading)). Matching members one by one cannot see any of this. It
finds no member of `c.[7C>T;8T>G]` in `c.7_8delinsTG`, and in the repeat below it pairs the two deletions, which happen
to be one allele, and fails on the substitutions, which are not, although the two descriptions make the same sequence.

So weaver applies both descriptions and compares what they make. The comparison need not cover the whole reference:
outside the stretch from the first base either description touches to the last, both are the reference, so the two
sequences are equal exactly when they are equal inside it. Nothing about the reference widens the stretch: a deletion
written at one end of a run and its equal written at the other each touch their own end, so the stretch already holds
both. The comparison does not build the edited sequences, and what it reads is bounded by the edits and the repeats they
sit in, the bases full justification already reads for one indel, so two members a megabase apart cost no more than two
side by side: the bases between them, which both descriptions leave as the reference, are compared by position, not
read.

A nucleotide allele is compared on the genome. Each member is projected to the genome on its own, and the two projected
descriptions are applied there. A transcript is used only when the provider has no genome for it, and then only when
nothing the comparison reads reaches outside the exons or crosses an exon junction or an alignment gap: no member, and
no region a member could equally be written over, since a deletion in a run that spans a junction is in either exon.
Otherwise the comparison is refused. Comparing on the transcript whenever both descriptions are on one would answer a
different question. On `NM_S.1`, which is `NM_X.1` with exon 1 ending at c.5, `c.[4A>T;6A>T]` and `c.4_6delinsTAT` make
one transcript. The first projects to `g.[19A>T;31A>T]`; the second is a `c.` change across a junction, which
[projects over the intron](#a-change-across-a-splice-junction-has-no-genomic-form), to `g.19_31delinsTAT`. Each
description is the same change as its own projection, and the two projections are different genomes; had the transcript
made the two descriptions the same change, equivalence would stop being transitive. A single variant meets the same
case: c.4_6 of `NM_S.1` is `AA|A`, and `c.5del` and `c.6del` make one transcript but delete g.20 and g.31.

Projection keeps the reading the allele has where it is written. A duplication inserts its copy 3' of the original on
its own strand, so on the minus strand the copy lands 5' of the original on the genome. `NR_XR.1:n.[23_25dup;22_23insA]`
is accepted on the transcript, the A going before what is copied, as in `c.[7_9dup;6_7insA]`. Written as `g.26_28dup`,
the duplication would insert at g.28_29, where the A, projected as `g.28_29insT`, also goes, and the pair would be
[ambiguous](#overlapping-members-of-a-cis-allele-need-one-reading). So projection is member by member, except that a
minus-strand duplication beside an insertion at the junction its genomic form would use is written on the genome as the
insertion it is, `g.25_26insCCT`, and the allele projects to `g.[25_26insCCT;28_29insT]`, which makes what the
transcript reading makes and normalises to `g.28_29insCCTT`.

A one-member allele makes what its member makes, so `c.[7C>T]` is `c.7C>T`, and an allele and a single variant are the
same change when they make the same sequence. The brackets of a one-member allele are spelling, as the parentheses of a
prediction are, so the two are Identity, as `c.=` and `c.[=]` are.

An allele says more than the sequence it makes when some of its members change nothing. An identity member, `c.10_12=`,
states that its bases are unchanged, and members that cancel each other state the same of the bases they are written
over. So a nucleotide description has two parts to compare, the sequence it makes and the bases it states, and two
descriptions are the same when both match, projected to the genome. `c.[7C>T;10_12=]` makes what `c.7C>T` makes and also
states c.10 to c.12, so the two are Different. A description that makes the reference is all statement: `c.10=` and
`c.20=` are Different, and `c.10=` and `NC_X.1:g.25=`, its projection, are Identity. Members that cancel name the bases
they are written over, an insertion naming its two flanks, with no widening over a repeat, as for the comparison's
stretch: on `NM_X.1`, whose c.4_6 is `AAA`, `c.[4_5insA;5del]` is `c.4_5=` and not `c.4_6=`, though c.6 is in the same
run, and beside a change they still state it, so `c.[4_5insA;5del;13T>G]` is `c.[4_5=;13T>G]` and not `c.13T>G`. Where
cancelling members are written is part of what they state: on the same run, `c.[5_6insA;6del]` is `c.5_6=`, so it is
Different from `c.[4_5insA;5del]`, although the two make one sequence and each pair, read as a change, would be written
at any place in the run. Which members count as cancelling, and the cases where none do, is set out with
[normalisation](#members-of-a-cis-allele-that-meet-are-written-as-one-change), which writes them as identities. Whether
members cancel is a fact about the reference, not the text: on `TTTTTTTTTT` the same pair makes `g.5T>A`, a change, and
states nothing.

The comparison is the equality of one value, the pair of sequence and stated bases, so equivalence between nucleotide
descriptions is transitive; across nucleotide and protein, and through a protein residue left unnamed, it is not, as the
paragraph on frameshifts below and [the bracketing entry](#a-prediction-is-bracketed-as-hgvs-writes-it) record. The
alternatives were to ignore what is stated, which would make a check at c.10 the same as a check at c.20, and to ignore
it only beside a change, which would make `c.[7C>T;10_12=]`, a screen's record that c.10 to c.12 were assayed in that
combination, the same as `c.7C>T`, which records nothing about them.

Read from the sequence, a protein consequence comes from the codons all the members leave, as the delins page asks:
`c.[4A>T;6A>T]` turns the codon `AAA` into `TAT`. An allele that makes the reference predicts what its statement
predicts. `c.[4_5insA;5del]`, read as `c.4_5=`, is `p.(Lys2=)`, silent at the residue its bases are in, as
[an identity inside the codons](#edits-outside-the-cds-are-statements) is; in a UTR it is `p.(=)`. `c.[10=;20=]`,
statements in two codons, is `p.[(Ala4=;Ter7=)]`, the predicted allele written as the protein alleles page writes one,
the parentheses inside the brackets. An allele that changes the sequence and states bases predicts both: each stated
base in a codon the change leaves as it was adds its residue, silent, to the predicted allele, so `c.[13T>G;10_12=]` is
`p.[(Ala4=;Tyr5Asp)]` and `c.[4_5insA;5del;13T>G]` is `p.[(Lys2=;Tyr5Asp)]`. A stated base in a codon the change alters
adds nothing, since the change names that residue: `c.[4_5insA;5del;6A>T]` is `p.(Lys2Asn)`. Nor does one outside the
codons, or in a codon read in a new frame or past a new stop; `p.(=)` is about the whole protein and does not stand
beside a change.

Some members cannot be applied to a sequence: a statement (the whole-sequence `=`, `?`, `spl`, `0`, `r.(?)`); a position
that is uncertain or imprecise, `g.(?_5)_(10_?)del`; inserted bases that are not stated, `insN[20]`; a copy number; and
a conversion, whose bases weaver does not model, and which is not an identity though it prints as one. A positional
identity, `c.10_12=`, is not among them: it replaces its bases with themselves, so `c.[7C>T;10_12=]` makes what `c.7C>T`
makes, and states c.10 to c.12 besides. An allele with members that cannot be applied is compared in two parts: the
members that can be applied by the sequence they make, as above, and the members that cannot as an unordered collection
of their projections to the genome, counted with repeats. Two alleles are the same change when both parts match. An
applicable part that makes the reference is compared as the statement it makes, so `c.[4_5insA;5del;10_12copy3]` is
`c.[4_5=;10_12copy3]` and not `c.[13=;10_12copy3]`. Order still does not matter, so `g.[22C>T;30_35copy3]` is
`g.[30_35copy3;22C>T]`, and a split of the applicable part still does not, so `c.[7C>T;8T>G;10_12copy3]` is
`c.[7_8delinsTG;10_12copy3]`; a member that cannot be applied matches only its own projection, so `30_35copy3` is not
`30_35copy4`. The projection, not the text as written, is compared, so `c.[7C>T;10_12copy3]` and `g.[22C>T;25_27copy3]`,
the same members on one genome, are one allele though `10_12copy3` and `25_27copy3` are different text. Such a member is
still held to the [overlap rule](#overlapping-members-of-a-cis-allele-need-one-reading) on its outer bounds, so
`g.[22C>T;20_30copy3]` is refused rather than compared in two parts.

A `p.` frameshift or extension can be applied. Equivalence already turns a single one into
[the protein it leaves](equivalence_logic.md#the-protein-a-description-leaves): residues it does not name match any
residue, and with no length the end is open. A member of an allele is read the same way, a frameshift replacing
everything from its residue on, so a member after it overlaps it. Compared as text instead, an allele would be judged by
a different rule from the variant it contains: on `NP_X.1`, `p.Leu3ProfsTer3` leaves `MK`, `P`, one unnamed residue,
then the end, and is Analogous with `p.Leu3delinsProAlaTer`, which leaves `MKPA`; but `p.[Lys2Gly;Leu3ProfsTer3]`,
leaving `MGP` and one residue, would be Different from `p.[Lys2Gly;Leu3delinsProAlaTer]`, leaving `MGPA`, by text alone.
Read as proteins, the two alleles are Analogous as their frameshift and delins are. The match through an unnamed residue
is "can be the same", so it is not transitive: `p.Leu3ProfsTer3` is Analogous with both `p.Leu3delinsProAlaTer` and
`p.Leu3delinsProGlyTer`, which are Different from each other.

The VRS identifier follows the normalised description, not the sequence. `to_vrs`
[normalises a cis allele first](#members-of-a-cis-allele-that-meet-are-written-as-one-change), so members that meet
become one Allele, and an allele whose members all meet is rendered as the Allele of the plain variant: `c.[7C>T;8T>G]`
and `c.7_8delinsTG` share one identifier. Members that stay apart are rendered as a `CisPhasedBlock` of their Alleles,
so `c.[5A>G;9G>T]` and `c.5_9delinsGACTT`, the same change, still have different identifiers; equivalence is where those
are found to name one change. Other producers build a block from members translated one at a time (MaveDB, for one,
translates each member with ga4gh.vrs and assembles the block), so weaver's identifier differs from theirs wherever
members meet: MaveDB's stored blocks for adjacent multi-variants, `c.[7C>T;8T>G]` among them, will not match weaver's
Allele. The block's identifier does not depend on member order
([below](#cisphasedblock-members-are-sorted-before-digesting)).

**Example** (`NM_X.1` is `GGGGG` `ATG AAA CTG GCC TAT CGC TAA` `CCGTATAAGTAAGG`, Met Lys Leu Ala Tyr Arg Ter, at genomic
index 10 of `NC_X.1`, so c.7 is g.22; `NR_XR.1` reads the same forty bases of `NC_X.1` on the minus strand, so n.23 is
g.28 and n.29 is g.22; `NP_X.1` is `MKLAYR`; `NM_S.1` is `NM_X.1` with exon 1 ending at c.5 and a ten-base intron after
it in `NC_S.1`, so c.4 is g.19 and c.6 is g.31):

```text
c.1_21 of NM_X.1                       ATG AAA CTG GCC TAT CGC TAA
NM_X.1:c.[7C>T;8T>G]                   ATG AAA TGG GCC TAT CGC TAA
NM_X.1:c.7_8delinsTG                   ATG AAA TGG GCC TAT CGC TAA    Analogous: the same sequence
NM_X.1:c.[8T>G;7C>T]                   ATG AAA TGG GCC TAT CGC TAA    Analogous: order is how it is written
NM_X.1:c.7_8delinsTA                   ATG AAA TAG GCC TAT CGC TAA    Different from c.[7C>T;8T>G]

NM_X.1:c.[7C>T;13T>G]                  ATG AAA TTG GCC GAT CGC TAA
NC_X.1:g.[22C>T;28T>G]                 the same bases                 Identity: the members project exactly to it
NR_XR.1:n.[23A>C;29G>A]                the same bases                 Analogous: the minus strand reverses the order
NM_X.1:c.7C>T                          ATG AAA TTG GCC TAT CGC TAA    Different from c.[7C>T;13T>G]: a member short
NM_X.1:c.[7C>T]                        ATG AAA TTG GCC TAT CGC TAA    Identity with c.7C>T
NM_X.1:c.[7C>T;10_12=]                 ATG AAA TTG GCC TAT CGC TAA    Different from c.7C>T: it also states c.10_12
NM_X.1:c.10=  vs  NM_X.1:c.20=         the reference for both         Different: statements about c.10 and c.20
NM_X.1:c.10=  vs  NC_X.1:g.25=         the reference for both         Identity: it projects exactly to it
NM_X.1:c.[4_5insA;5del]  vs  c.4_5=    the reference for both         Analogous: the members name c.4 and c.5
NM_X.1:c.[4_5insA;5del]  vs  c.4_6=    the reference for both         Different: c.6 is in the run, but not named
NM_X.1:c.[4_5insA;5del]  vs  c.10=     the reference for both         Different
NM_X.1:c.[4_5insA;5del]  vs  c.[5_6insA;6del]                         Different: c.4_5= and c.5_6=, placed where written
NM_X.1:c.[4_5insA;5del;13T>G]          ATG AAA CTG GCC GAT CGC TAA    Analogous with c.[4_5=;13T>G]
                                                                      Different from c.13T>G: c.4 and c.5 stated
on TTTTTTTTTT:
X:g.[4_5insA;5del]  vs  X:g.5T>A       TTTTATTTTT for both            Analogous: here the members change a base

NM_X.1:c.[4_5insA;5del]                →  NP_X.1:p.(Lys2=)            the statement c.4_5=, silent at Lys2
NM_X.1:c.[10=;20=]                     →  NP_X.1:p.[(Ala4=;Ter7=)]
NM_X.1:c.[13T>G;10_12=]                →  NP_X.1:p.[(Ala4=;Tyr5Asp)]  the change and the statement
NM_X.1:c.[4_5insA;5del;6A>T]           →  NP_X.1:p.(Lys2Asn)          the stated bases are in the changed codon

NM_X.1:c.[4A>T;6A>T]                   ATG TAT CTG GCC TAT CGC TAA    →  NP_X.1:p.(Lys2Tyr)
NM_X.1:c.4_6delinsTAT                  ATG TAT CTG GCC TAT CGC TAA    Analogous with c.[4A>T;6A>T]
                                                                      not p.[(Lys2Ter;Lys2Asn)], what the members predict alone

on TTCAGCAGTT:
X:g.[3_5del;7A>T]                      TTCTGTT
X:g.[4A>T;6_8del]                      TTCTGTT                        Analogous: the same sequence, though 7A>T and 4A>T differ alone
X:g.4_7delinsT                         TTCTGTT                        Analogous

NP_X.1:p.[Lys2Gly;Leu3Val]             MGVAYR
NP_X.1:p.Lys2_Leu3delinsGlyVal         MGVAYR                         Analogous
NP_X.1:p.Lys2Gly                       MGLAYR                         Different
NP_X.1:p.[Lys2Gly;Leu3ProfsTer3]       MGP, one residue, end
NP_X.1:p.[Lys2Gly;Leu3delinsProAlaTer] MGPA, end                      Analogous: the unnamed residue matches A

NC_X.1:g.[22C>T;30_35copy3]  vs  NC_X.1:g.[30_35copy3;22C>T]          Analogous: 22C>T by sequence, the copy number by projection
NC_X.1:g.[22C>T;30_35copy3]  vs  NC_X.1:g.[22C>T;30_35copy4]          Different: three copies, four copies
NM_X.1:c.[7C>T;10_12copy3]   vs  NC_X.1:g.[22C>T;25_27copy3]          Analogous: the copy number projects to that text
NC_X.1:g.[22C>T;20_30copy3]                                           refused: g.22 is inside the counted stretch

c.4_9 of NM_S.1                        AA|ACTG                        exon 1 ends at c.5
NM_S.1:c.[4A>T;6A>T]                   TA|TCTG                        Identity with NC_S.1:g.[19A>T;31A>T]
NM_S.1:c.4_6delinsTAT                  TA|TCTG                        Different from c.[4A>T;6A>T]: on the genome
                                                                      it is g.19_31delinsTAT, the intron deleted
NM_S.1:c.5del  vs  NM_S.1:c.6del       AACTG for both                 Different: g.20del and g.31del

to_vrs_variation NM_X.1:c.[7C>T;8T>G]    →  an Allele, ga4gh:VA.…, that of c.7_8delinsTG: the members touch
to_vrs_variation NM_X.1:c.[5A>G;9G>T]    →  a CisPhasedBlock of two Alleles, ga4gh:CPB.…
to_vrs_variation NM_X.1:c.5_9delinsGACTT →  an Allele, ga4gh:VA.…     the same change, another identifier
```

- **Differs:** mavehgvs 0.8.1 compares the members as written: `c.[7C>T;8T>G]` is not `c.7_8delinsTG`, and members out
  of order are a parse error unless `relaxed_ordering` sorts them.
- **Others:** biocommons 1.5.7 does not parse alleles in cis; VariantValidator and ClinVar not checked.
- **Spec (silent):** HGVS describes, it does not compare. What the comparison rests on:
  [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/): "a series of variants on one
  chromosome", and, in its answer on a simplified haplotype format, listed "in genomic order";
  [deletion-insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) for `c.[235A>T;237G>T]`, the
  "alternative description" of `c.850_901delinsTTCCTCGATGCCTG` and `c.2077delinsATA`. It says which spelling to write,
  not how to compare two; that reading is weaver's.
- **Tests:** `cis_phased_test::cis_alleles_compare_as_sets_of_their_members_alleles`, whose rows on order, the genome
  spelling and a missing or extra member hold under this rule, and whose rows comparing a copy-number member as text are
  to change, and `::statements_about_different_bases_differ`; planned: `::adjacent_members_are_the_delins_they_make`,
  `::members_in_either_order_are_one_allele`, `::members_on_the_minus_strand_are_the_genomic_allele`,
  `::members_in_one_repeat_match_as_the_sequence_they_make`, `::a_cis_allele_is_not_one_of_its_members`,
  `::a_one_member_cis_allele_is_its_member`, `::protein_members_are_the_delins_they_make`,
  `::members_in_one_codon_predict_one_residue`, `::members_without_a_sequence_match_as_text_in_any_order`,
  `::members_on_one_transcript_compare_without_the_genome`, `::the_vrs_identifier_follows_the_description`,
  `::a_cancelling_allele_states_the_bases_its_members_name`, `::cancelling_depends_on_the_reference`,
  `::a_change_is_never_a_statement`, `::an_identity_beside_a_change_is_part_of_the_allele`,
  `::cancelling_members_beside_a_change_are_a_statement`, `::a_change_with_a_statement_predicts_both`,
  `::a_cancelling_allele_predicts_its_statements_protein`; `equivalence_is_transitive_over_statements_and_changes`
  (property).

### Overlapping members of a cis allele need one reading

**Members of a cis allele that reach for the same reference base are read every way the description allows. With no
reading the allele is contradictory, a `ValidationError` naming the two members; with readings that make different
sequences it is [ambiguous](#an-ambiguous-description-is-an-error-that-lists-its-readings), an
`AmbiguousDescriptionError` that lists each one; with one it is accepted.**

Because every member is numbered against the reference, two members that reach for the same reference base describe a
chromosome that either cannot exist or cannot be told from another, and the two cases call for different answers.

**Contradictory: no reading.** `c.[7_9del;8T>G]` deletes c.8 and changes it: on that chromosome there is no c.8 to
become G. `c.[7C>T;7C>G]` puts two bases at one position of one chromosome. `c.[7_9=;8T>G]` says c.8 is both unchanged
and changed, and `c.[7_8=;7_8insA]` says c.7 and c.8 are unchanged, which includes that they are still side by side, and
puts an A between them. `c.[7_9delinsA;8T>G]` replaces three bases with one, which is not any of the three, so there is
no base left for c.8 to be, and `c.[7_9delinsTT;8_9insA]` puts an A between c.8 and c.9, a junction no base of `TT`
stands on either side of. No sequence answers any of these, so each is a `ValidationError`.

**Ambiguous: several readings.** `c.[7_8insA;7_8insG]` puts an A and a G between c.7 and c.8 but not in which order, and
`c.7_8insAG` and `c.7_8insGA` are different sequences. A duplication inserts at the junction 3' of what it copies, which
is how the specification defines it, so `c.[7_9dup;9_10insA]` is two insertions at one junction too.
`c.[7_8insA;7_8insA]` is one insertion written twice, `c.7_8insA`, or two, `c.7_8insAA`, and `c.[7_9dup;9_10insCTG]`
inserts `CTG` at c.9_10 twice over in the same way. `c.[7_9dup;8T>G]` changes a base the duplication copies: the copy is
of the reference, which normalises to `c.7_8insGGC`, or of the edited bases, `c.8delinsGGCG`. An inversion reads the
bases it covers in the same way. `c.[7_9inv;8T>A]` inverts `CTG`, read as written, to `CAG`, the substitution's A then
landing where an A already is, or inverts the changed bases `CAG` back to `CTG`, the reference; `c.[7_9inv;8T>G]` is
`CGG` or `CCG` by the same two readings. The same holds for any member that reads a stretch another member changes: a
conversion, which the specification no longer uses but weaver parses, and an insertion of a copied range (`ins5_16`),
which weaver does not yet parse. An insertion between two bases another member deletes, or replaces keeping the length,
is ambiguous too: `c.[7_9del;7_8insA]` either loses the A with the bases around its junction, `c.7_9del`, or keeps it in
their place, `c.7_9delinsA`, and `c.[7_8delinsTG;7_8insA]` is `c.7_8delinsTG` or `c.7_8delinsTAG` by the same two
readings. Each of these has answers, and weaver does not pick one: it returns them all.

**One reading: accepted.** Nearness is not the problem: `c.[7_8del;9_10del]` is two deletions meeting at a junction,
together a deletion of four bases, and `c.[2077G>A;2077_2078insTA]`, the allele the specification's BRCA1 question asks
about, changes a base and inserts beside it. Members that agree on every base they share have one reading, and the
redundancy is dropped: one edit written twice, `c.[7C>T;7C>T]` or `c.[7C>T;7delinsT]`; overlapping edits that put the
same base at each shared position, `c.[7_8delinsTG;7C>T]`, which is `c.7_8delinsTG`; overlapping deletions, which remove
each shared base once, `c.[7_9del;8_10del]`, which deletes c.7 to c.10; and identities, `c.[7_9=;8=]`, which is
`c.7_9=`. An insertion at the outer flank of a deleted or replaced stretch meets it rather than overlapping it:
`c.[7_9del;6_7insA]` puts the A where c.7_9 were, `c.7_9delinsA`, and `c.[7_8delinsTG;8_9insA]` puts it after the G,
`c.7_8delinsTGA`. Insertions at one junction whose orders make one sequence and that are not identical,
`c.[7_8insA;7_8insAA]`, have one reading too.

So each member is read as the reference interval it replaces and the bases it puts there; an insertion replaces the
empty interval at a junction between two bases, and a duplication is an insertion at the junction 3' of the bases it
copies. A substitution, a deletion, a delins that keeps the length and a positional identity are read base by base: each
reference base gets one base, or none. A delins that changes the length gives its bases no such reading and replaces its
stretch as a whole. A duplication and an inversion read the bases they cover. Identical members other than insertions
and duplications are one edit written twice, and are read once before anything else is checked. The allele is:

- contradictory when two members read base by base disagree at a base, one keeping or changing a base the other changes
  or removes; when a delins that changes the length shares a base with another member, or another member inserts inside
  it, since no base of the delins corresponds to the junction; when a member inserts inside a positional identity's
  stretch; or when an identity names a base an inversion changes; contradictory comes first when both apply;
- ambiguous when a member inserts between two bases another deletes, or replaces base by base keeping the length; when
  two members insert at one junction and their orders make different sequences, when identical insertions could be one
  written twice, or when a member changes, removes or inserts among the bases a duplication or an inversion reads;
- accepted otherwise: `c.[7_9delinsA;7_9delinsA]`, read once, is as accepted as `c.[7C>T;7C>T]`. Where a change agrees
  with an identity, the identity stays and the change is cut where it stands: `c.[7=;7_8delinsCG]` is written
  `c.[7=;8T>G]`, and `c.[8=;7_9delinsTTA]`, agreeing in the middle, `c.[7C>T;8=;9G>A]`.

A repeat replaces its [whole run](#a-repeat-resolves-to-its-whole-run). A positional identity, `c.10_12=`, is an edit
that replaces its bases with themselves: it can be applied, it takes part in this rule, and normalisation never drops a
base it names. Combinatorial variant screens use it to record that a position in a given combination was assayed and is
wild type. That record is part of the allele, beside a change as alone, so `c.[7C>T;10_12=]` is Different from `c.7C>T`
([the cis allele entry](#a-cis-allele-is-the-sequence-its-members-make)). The whole-sequence `c.=` beside a change says
no change was found and names one, a contradiction, and stays refused
([the `c.=` entry](#the-whole-transcript-unchanged-is-a-statement-not-a-change)).

A member that cannot be applied, a copy number or a position that is uncertain or imprecise, still names a stretch, and
is held to this rule on its outer bounds. `g.[22C>T;20_30copy3]` changes a base inside the stretch whose copies it
counts, and `g.[(?_5)_(10_?)del;7A>T]` changes g.7, which lies inside even the bases the deletion certainly removes:
both are contradictory. An unknown outer bound, the `?` of `(?_5)`, reaches the end of the sequence on that side, as a
bound for this check only: nothing is built or read for it.

One case sits on the line. `c.[7_8insA;8T>G]` inserts between c.7 and c.8 and changes c.8. Read against the reference
the insertion is at a junction and c.8 is a base, so they share nothing, and the result has one reading, the same as
`c.7_8delinsCAG`; it is accepted, the mirror image of `c.[2077G>A;2077_2078insTA]`. A reading that counts an insertion's
flanking bases as part of it, as mavehgvs does, refuses both. The specification answers that the 2077 change is written
`c.2077delinsATA` and has withdrawn an earlier note calling the allele form acceptable, so accepting it is weaver's
choice, made because the allele has one reading, and refusing it departs from no recommendation.

The rule reads the members as written. Shifting them one at a time could make an accepted allele overlap, which is why
normalisation writes [members that meet as one change](#members-of-a-cis-allele-that-meet-are-written-as-one-change).

**Example** (`NM_X.1` as [above](#a-cis-allele-is-the-sequence-its-members-make); its c.4_12 reads `AAACTGGCC`; each
reading normalised):

```text
                                c.4_12 becomes
NM_X.1:c.[7_8del;9_10del]       AAACC                          accepted, c.8_11del
NM_X.1:c.[7C>T;7_8insTA]        AAATTATGGCC                    accepted, c.7delinsTTA: the 2077 case
NM_X.1:c.[7_8insA;8T>G]         AAACAGGGCC                     accepted, c.8delinsAG
NM_X.1:c.[7_9dup;6_7insA]       AAAACTGCTGGCC                  accepted, c.6_7insACTG: the A goes before what is copied
NM_X.1:c.[7C>T;7C>T]            AAATTGGCC                      accepted, c.7C>T
NM_X.1:c.[7_8delinsTG;7C>T]     AAATGGGCC                      accepted, c.7_8delinsTG: the two agree at c.7
NM_X.1:c.[7_9del;8_10del]       AAACC                          accepted, c.8_11del: c.8 and c.9 deleted once
NM_X.1:c.[7_9del;6_7insA]       AAAAGCC                        accepted, c.7_9delinsA: the A at the deletion's flank
NM_X.1:c.[7_8delinsTG;8_9insA]  AAATGAGGCC                     accepted, c.7_8delinsTGA: the A after the G
NM_X.1:c.[7_9=;8=]              AAACTGGCC                      accepted, c.7_9=
NM_X.1:c.[7=;7_8delinsCG]       AAACGGGCC                      accepted, c.[7=;8T>G]: the identity stays
NM_X.1:c.[8=;7_9delinsTTA]      AAATTAGCC                      accepted, c.[7C>T;8=;9G>A]: the delins cut at c.8
NM_X.1:c.[7_9delinsA;7_9delinsA] AAAAGCC                       accepted, c.7_9delinsA: one edit written twice
NM_X.1:c.[7_8insA;7_8insAA]     AAACAAATGGCC                   accepted, c.7_8insAAA

NM_X.1:c.[7C>T;10_12=]          AAATTGGCC                      accepted, left as written; Different from c.7C>T

NM_X.1:c.[7C>T;7C>G]                                           contradictory: two bases at c.7
NM_X.1:c.[7_9del;8T>G]                                         contradictory: no c.8 to change
NM_X.1:c.[7_8delinsTG;8del]                                    contradictory: G at c.8, and no c.8
NM_X.1:c.[7_9=;8T>G]                                           contradictory: c.8 unchanged and changed
NM_X.1:c.[7_8=;7_8insA]                                        contradictory: c.7 and c.8 unchanged, and parted
NM_X.1:c.[7_9inv;8=]                                           contradictory: the inversion makes c.8 an A
NM_X.1:c.[7_9delinsA;8T>G]                                     contradictory: the delins changes the length over c.8
NM_X.1:c.[7_9delinsTT;8_9insA]                                 contradictory: an insertion inside such a delins
NC_X.1:g.[22C>T;20_30copy3]                                    contradictory: g.22 is inside the counted stretch

NM_X.1:c.[7_8insA;7_8insG]      AAACAGTGGCC | AAACGATGGCC      ambiguous: c.7_8insAG | c.7_8insGA
NM_X.1:c.[7_8insA;7_8insA]      AAACATGGCC | AAACAATGGCC       ambiguous: c.7_8insA | c.7_8insAA
NM_X.1:c.[7_9dup;8T>G]          AAACGGCTGGCC | AAACGGCGGGCC    ambiguous: c.7_8insGGC (a copy of the reference)
                                                                        | c.8delinsGGCG (a copy of the edited bases)
NM_X.1:c.[7_9dup;9_10insA]      AAACTGACTGGCC | AAACTGCTGAGCC  ambiguous: c.6_9dup | c.9_10insCTGA
NM_X.1:c.[7_9dup;9_10insCTG]    AAACTGCTGGCC                   ambiguous: c.7_9dup, one copy written twice
                                | AAACTGCTGCTGGCC                         | c.9_10insCTGCTG, two copies
NM_X.1:c.[7_9inv;8T>A]          AAACAGGCC | AAACTGGCC          ambiguous: c.8T>A (the reference inverted)
                                                                        | c.7_9= (the changed bases inverted)
NM_X.1:c.[7_9inv;8T>G]          AAACGGGCC | AAACCGGCC          ambiguous: c.8T>G | c.8T>C
NM_X.1:c.[7_9del;7_8insA]       AAAGCC | AAAAGCC               ambiguous: c.7_9del | c.7_9delinsA
NM_X.1:c.[7_8delinsTG;7_8insA]  AAATGGGCC | AAATAGGGCC         ambiguous: c.7_8delinsTG | c.7_8delinsTAG
```

- **Agrees:** mavehgvs 0.8.1 refuses every case found contradictory or ambiguous above that was checked against it
  (`c.[7C>T;7C>G]`, `c.[7_9del;8T>G]`, `c.[7_8delinsTG;8del]`, `c.[7_9delinsA;8T>G]`, `c.[7_8insA;7_8insG]`,
  `c.[7_9dup;8T>G]`, `c.[7_9dup;9_10insA]`, `c.[7_9del;7_8insA]`, `c.[7_8delinsTG;7_8insA]`, with "multi-variant has
  overlapping changes" or "multiple changes at same position") and accepts `c.[7_8del;9_10del]`. The other rows: not
  checked; it does not parse `inv`.
- **Differs:** mavehgvs counts an insertion's flanking positions as part of it and refuses any overlap, so it also
  refuses the accepted `c.[7C>T;7_8insTA]`, `c.[7_8insA;8T>G]`, `c.[7_9dup;6_7insA]`, `c.[7_9del;8_10del]`,
  `c.[7_9del;6_7insA]` and `c.[7_8delinsTG;8_9insA]`, and `c.[2077G>A;2077_2078insTA]`, for which the specification
  recommends `c.2077delinsATA`; it refuses any `=` among the members. `bcftools consensus`, applying records rather than
  an allele, does not refuse: it applies the first of two overlapping records and skips the second with "The site …
  overlaps with another variant, skipping…"
  ([consensus.c](https://github.com/samtools/bcftools/blob/develop/consensus.c)).
- **Others:** biocommons 1.5.7 does not parse alleles in cis; VariantValidator and ClinVar not checked.
- **Spec (differs):** on a positional identity in cis. The
  [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) page, in a note on descriptions of
  variants in trans, calls `LRG_199t1:c.[2376G>C;3103=];[2376=;3103del]` "not correct", and in a separate note on the
  `NM_004006.2` example says "the indication "not changed" is used only when one variant was identified". Both are about
  trans descriptions, and the page's own grammar admits a positional `=` inside an allele; weaver reads them as
  discouraging `=` beside a change and accepts it in cis, for the screens that record a wild-type position. The page
  does not address members that overlap;
  [duplication](https://hgvs-nomenclature.org/stable/recommendations/DNA/duplication/): "a copy of one or more
  nucleotides is inserted directly 3' of the original copy";
  [deletion-insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/delins/) answers the 2077 question with
  `c.2077delinsATA` and says "The previous format "con" is no longer used".
- **Tests:** none yet; planned: `cis_phased_test::members_meeting_at_a_junction_are_accepted`,
  `::an_insertion_beside_a_changed_base_is_accepted`, `::an_insertion_before_a_changed_base_is_accepted`,
  `::an_insertion_before_what_a_duplication_copies_is_accepted`, `::members_that_agree_on_a_base_are_one_change`,
  `::overlapping_deletions_delete_each_base_once`, `::an_insertion_inside_a_same_length_delins_is_ambiguous`,
  `::an_insertion_inside_a_deletion_is_ambiguous`, `::an_insertion_at_the_flank_of_a_deletion_is_accepted`,
  `::an_insertion_at_the_flank_of_a_delins_is_accepted`, `::an_identity_beside_an_agreeing_change_stays`,
  `::a_change_agreeing_with_an_identity_inside_it_is_cut_there`, `::an_edit_written_twice_is_one_edit`,
  `::two_changes_to_one_base_are_contradictory`, `::a_change_inside_a_deletion_is_contradictory`,
  `::a_change_under_an_identity_is_contradictory`, `::an_insertion_inside_an_identity_is_contradictory`,
  `::a_delins_that_changes_the_length_under_another_member_is_contradictory`,
  `::a_member_without_a_sequence_is_held_to_its_outer_bounds`, `::two_insertions_at_one_junction_are_ambiguous`,
  `::identical_insertions_at_one_junction_are_ambiguous`, `::a_change_to_what_a_duplication_copies_is_ambiguous`,
  `::a_change_to_what_an_inversion_reads_is_ambiguous`, `::an_identity_inside_an_inversion_is_contradictory`,
  `::an_insertion_where_a_duplication_inserts_is_ambiguous`, `::the_error_names_both_members`.

### An ambiguous description is an error that lists its readings

**Where a description has several readings, weaver does not pick one. It raises `AmbiguousDescriptionError` in Python,
`HgvsError::Ambiguous` in Rust, distinct from a `ValidationError`, which carries every reading as an unambiguous
description, so a caller can test for ambiguity, catch it and choose.**

weaver gives an answer or an error, never a guess ([`PRODUCT.md`](../PRODUCT.md)). A contradictory description has no
answer, and the error says so. An ambiguous one has several, and choosing among them is the guess the principle rules
out; refusing it with a plain error would throw away what is known, that the change is one of a short list. So the error
is the list. A caller that knows which reading was meant, from the assay or the source, takes it; one that does not
reports the ambiguity rather than a sequence nobody wrote.

The readings come from [cis alleles whose members overlap](#overlapping-members-of-a-cis-allele-need-one-reading): an
insertion inside a deleted or replaced stretch, two insertions at one junction, identical insertions, a change to what a
duplication copies or an inversion reads. They also come from
[cancelling members](#members-of-a-cis-allele-that-meet-are-written-as-one-change) that can be paired more than one way:
the readings then make one sequence, and state different bases. A `c.` range that crosses an exon junction
[projects over the intron](#a-change-across-a-splice-junction-has-no-genomic-form); whether it is instead ambiguous
between including the intron and not is issue #64. A `CisPhasedBlock` whose members' placements make different sequences
is another source, [read back](#members-of-a-cis-allele-that-meet-are-written-as-one-change) by `from_vrs`.

VRS meets the same cases from the other side. A `CisPhasedBlock`'s members are unordered (VRS 2.0.1 marks `members`
`ordered: false`, though its description calls the block "an ordered set"), so `c.[7_8insA;7_8insG]` and
`c.[7_8insG;7_8insA]` would make one block, which silently conflates what HGVS leaves ambiguous. `to_vrs` therefore
refuses an ambiguous allele before a block is built, and `from_vrs`, reading two members that insert at one interbase
point, returns the ambiguity with both orders as its readings, where the orders make different sequences, and with one
copy and two as its readings where the two members are identical. VRS states a duplication's inserted bases against the
reference, so which bases are copied does not arise inside VRS; only its overlap with another member does.

**Example** (`NM_X.1` as [above](#a-cis-allele-is-the-sequence-its-members-make)):

```text
NM_X.1:c.[7_8insA;7_8insG]       →  AmbiguousDescriptionError: NM_X.1:c.7_8insAG | NM_X.1:c.7_8insGA
NM_X.1:c.[7_9dup;8T>G]           →  AmbiguousDescriptionError: NM_X.1:c.7_8insGGC | NM_X.1:c.8delinsGGCG
NM_X.1:c.[3_4insA;5del;6_7insA]  →  AmbiguousDescriptionError: NM_X.1:c.[3_5=;6dup] | NM_X.1:c.[4dup;5_7=]
NM_X.1:c.[7C>T;7C>G]         →  ValidationError: contradictory, no reading
```

- **Others:** mavehgvs 0.8.1 refuses `c.[7_8insA;7_8insG]` and `c.[7_9dup;8T>G]` as overlapping changes and gives no
  readings; other tools not checked.
- **Spec (silent):** HGVS does not say what a tool does with a description it cannot read one way.
- **Tests:** none yet; planned: `cis_phased_test::the_ambiguity_lists_every_reading`,
  `::an_ambiguity_is_not_a_validation_error`, `::an_ambiguity_carries_its_readings`,
  `tests/test_cis_phased.py::test_an_ambiguity_is_an_ambiguous_description_error_with_its_readings`,
  `::to_vrs_refuses_an_ambiguous_allele`;
  `vrs_round_trip_test::a_block_inserting_twice_at_one_point_reads_back_as_ambiguous`,
  `::a_block_of_two_identical_insertions_reads_back_as_ambiguous`.

## Alleles, SPDI and VRS

### The canonical allele is fully justified

**A change is trimmed to what it alters and then widened over the whole region in which it could equally be written, so
every spelling of one change is one allele with one identifier.**

`to_spdi_unambiguous` renders it. The plain `to_spdi` renders the 3'-normalised variant instead, with an insertion
placed at its second flank as SPDI counts.

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
- **Spec (agrees):** VRS 2.0.1 "Normalization" (fully justified alleles); refget/SPDI interbase coordinates.
- **Tests:** `src/allele.rs::deletion_anywhere_in_a_run_is_the_same_allele`,
  `::duplication_and_insertion_of_the_unit_are_the_same_allele`,
  `::literal_insertion_stays_literal_and_substitution_is_trimmed` (read from the 0-based ranges they assert);
  `transcript_coordinates_test::plain_spdi_places_an_insertion_at_its_second_flank`;
  `ambiguous_range_is_sound_and_complete` and `vrs_and_spdi_round_trip` (properties).

### Refget accessions are computed over the normalised sequence

**Letters uppercased, everything else dropped, as the refget specification defines.**

Without a `Refget` lookup the accession is computed from the whole sequence, and `from_vrs` then needs the accession
passed. The digest is checked against the sequence either way.

**Example:**

```text
ACGT             →  SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2     the specification's own example
acgt             →  the same
AC⏎G T⏎          →  the same
chr19 of the NCBI GRCh38 FASTA  →  SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl, the specification's published value
```

- **Spec (agrees):** [refget](https://samtools.github.io/hts-specs/refget.html), "Sequence Normalization" and the
  `ga4gh` identifier.
- **Tests:** `src/vrs.rs::refget_accession_is_over_the_normalised_sequence`,
  `::digest_and_identifier_match_the_spec_example`;
  `tests/test_digest_table.py::test_digest_follows_the_refget_normalisation_rule`;
  `vrs_round_trip_test::the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked`.

### Uncertain breakpoints become Range bounds

**A deletion with uncertain breakpoints renders as an Allele with Range bounds and an empty state, unnormalised; a
duplication with uncertain breakpoints has no Allele, since its bases are unknown, and renders as a CopyNumberChange
gain.**

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
  `(A_B)_(C_D)` "where `B_C` describes the minimal extent and `A_D` the maximal", `?` for unknown; VRS 2.0.1 `Range` and
  `CopyNumberChange`.
- **Tests:** `vrs_range_test::uncertain_breakpoints_become_ranges`,
  `::an_imprecise_deletion_has_an_empty_literal_state_and_json_ranges`,
  `::only_deletions_may_have_uncertain_breakpoints`, `::unknown_bounds_round_trip_as_question_marks`;
  `vrs_copy_change_test::an_imprecise_duplication_is_a_gain_over_range_bounds`,
  `::to_vrs_variation_gives_a_change_for_an_imprecise_dup_only`; `spec_summary_test::test_uncertain_intervals`.

### copyChange is a label

**VRS 2.0.1 defines `copyChange` as a string enum; weaver writes labels, and reads labels, bare EFO CURIEs and the 2.0.0
`MappableConcept` form.**

The specification's published `CopyNumberChange` example digest only reproduces when `copyChange` is digested as the
pre-release string `EFO:0030071`; both digests are pinned in tests so the machinery is anchored to a published value
while output follows the current schema. EFO's high-level loss is `EFO:0020073`, not in the `00300xx` block.

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

**The identifier does not depend on the order the members are written; the `members` array keeps the order of the
normalised allele `to_vrs` renders, by position along the sequence the allele is written on.**

The schema's `minItems: 2` is relaxed to one when a block is read, a departure from VRS, because HGVS allows the
degenerate `c.[145C>T]` and other producers render it as a block of one; `to_vrs` itself renders it as an Allele. The
in-trans form `[..];[..]` is rejected: it describes two molecules.

**Example:**

```text
NM_X.1:c.[7C>T;13T>G]   and  NM_X.1:c.[13T>G;7C>T]    →  the same ga4gh:CPB. id, members in the normalised order
NM_X.1:c.[7C>T;13T>G]   and  NC_X.1:g.[22C>T;28T>G]   →  the same id
NM_X.1:c.[7C>T]                                        →  an Allele: normalised, it is c.7C>T
NM_X.1:c.[7C>T];[13T>G] →  UnsupportedOperation: describes two molecules (alleles in trans); parse each [..] as a cis allele instead
```

- **Spec (differs on minItems):** VRS 2.0.1 computed identifiers, "order arrays of digests and ids by Unicode Character
  Set values"; [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/) for cis `[a;b]` against
  trans `[a];[b]`.
- **Tests:** `src/vrs.rs::cis_phased_block_digest_matches_the_spec_validation_data`,
  `::a_cis_phased_block_carries_its_members_in_the_order_given`, to change to the normalised order;
  `cis_phased_test::the_identifier_is_the_same_whichever_way_the_members_are_written`,
  `::alleles_in_trans_are_two_molecules_and_are_refused`.

### Reading back gives the normalised variant

**`from_vrs` and `from_spdi` give the 3'-normalised variant on the allele's own sequence.**

An insertion trimmed to before the first base becomes a delins of that base, since HGVS has no insertion before base 1;
a replacement that is the reverse complement of what it replaces reads back as `inv`; a `LengthExpression` comes back as
`insN[n]`; a Range of copies is refused, since HGVS has no spelling for it.

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

- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/) (no position before 1);
  [inversion](https://hgvs-nomenclature.org/stable/recommendations/DNA/inversion/);
  [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/) for `insN[n]`.
- **Tests:** `vrs_round_trip_test::nucleotide_variants_round_trip_to_their_normalised_form`,
  `::alleles_from_other_producers_parse`, `::an_insertion_before_the_first_base_reads_back_as_a_delins_of_that_base`;
  `vrs_length_test::length_expressions_read_back_as_the_recommended_spelling`;
  `vrs_copy_number_test::from_vrs_checks_the_count_it_is_given`; `vrs_and_spdi_round_trip` (property).

### Breakends and fusions are not rendered

**`Adjacency`, `Terminus` and `DerivativeMolecule` describe breakends and fusions this grammar has no input for.**

**Example:** none; there is no HGVS input that would produce them.

- **Spec (silent):** VRS 2.0.1.
- **Tests:** none.

## Parsing

### The grammar is checked against the biocommons table

**580 inputs over the 92 rules the two grammars share, agreeing on all but three.**

weaver accepts a terminator inside an amino acid sequence, which biocommons rejects and ClinVar writes. And biocommons
takes the `?` of `p.Met1?` as a substitution; weaver reads it with its own rule, as the statement
[above](#a-start-codon-change-is-written-specifically), so its substitution rule rejects a bare `?`. biocommons rejects
`p.(?)`, which the protein alleles page writes for a variant expected but not yet identified; weaver reads it as the
unknown effect, predicted.

**Example:**

```text
rule aat13_seq, input TerGly        →  accepted by weaver, rejected by biocommons
NP_…:p.Glu26_Glu27insTerGlu         →  parses; ClinVar writes it
rule pro_subst, input ?             →  rejected by weaver, accepted by biocommons; p.Met1? still parses
rule p_posedit, input (?)           →  accepted by weaver, rejected by biocommons
```

- **Differs:** biocommons, on those three inputs.
- **Agrees:** ClinVar.
- **Spec (silent):** the nomenclature pages above; the table is biocommons's reading of them.
- **Tests:** `grammar_test::biocommons_grammar_table` (the difference is listed in `KNOWN_DIFFERENCES`);
  `parser_round_trips_canonical_hgvs` and `parser_never_panics` (properties).

### Forms accepted beyond biocommons

**Alleles in cis for every coordinate system; predicted RNA changes and RNA statements; protein statements; uncertain
genomic breakpoints; insertions of a stated length; copy number; `r.` positions with `*`.**

**Example** (each parses and prints back as written):

```text
NM_004006.2:c.[145C>T;147C>G]      NP_X.1:p.[Lys2Leu;Ala4del]      NM_X.1:r.[7c>u;13u>g]
NM_R.1:r.(10c>g)                   NM_R.1:r.0   NM_R.1:r.?   NM_R.1:r.spl   NM_R.1:r.=
NP_X.1:p.0   NP_X.1:p.0?   NP_X.1:p.Met1?   NP_X.1:p.(Met1?)   NP_X.1:p.(?)   NM_R.1:r.(?)
NC_000009.11:g.(?_108337304)_(108337428_?)del
NC_TEST.1:g.10_11insN[20]   NC_TEST.1:g.10_11insN[(20_30)]   NC_TEST.1:g.10_12delinsN[5]   NC_TEST.1:g.10_11ins(20)
NC_000014.8:g.88401076_88459508copy4
NM_R.1:r.*5u>a
```

- **Differs:** biocommons rejects these.
- **Spec (agrees):** [alleles](https://hgvs-nomenclature.org/stable/recommendations/DNA/alleles/);
  [RNA substitution](https://hgvs-nomenclature.org/stable/recommendations/RNA/substitution/);
  [protein substitution](https://hgvs-nomenclature.org/stable/recommendations/protein/substitution/);
  [protein alleles](https://hgvs-nomenclature.org/stable/recommendations/protein/alleles/) for `p.(?)`;
  [uncertain positions](https://hgvs-nomenclature.org/stable/recommendations/uncertain/);
  [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/).
- **Tests:** `cis_phased_test::cis_alleles_parse_and_print_back_in_every_coordinate_system`;
  `rna_test::predicted_changes_keep_their_parentheses`, `::statements_about_the_transcript_round_trip`;
  `spec_summary_test::test_uncertain_intervals`, `::test_spec_summary_variants`;
  `statements_test::predicted_p_unknown_prints_back_as_parsed`, `::r_follows_the_dna_prints_back_as_parsed`;
  `vrs_length_test::the_older_spelling_gives_the_same_allele`, `::an_uncertain_length_is_a_range`.

### Recommended spellings on output

**`insN[20]` rather than `ins(20)`; `extTer8` rather than `ext*8`; a range of repeat copies as `[(70_80)]` rather than
`(70_80)`; three-letter residues rather than one-letter; bare `del` and `dup`.**

A description read in an older or shorter spelling is stored as the spelling it prints, so the printed text parses back
to the same value. A position, a substituted residue, and the residues of a frameshift or extension are written in
three-letter code: `p.M1?` is stored and printed as `p.Met1?`. Where only a stop can stand, `X` or `*` is `Ter`: after
`fs`, after `ext` with a count (`ext*5` and `extX+5`, not the offset of `extXaa-5`), and at the start of an extension
past the stop (`p.X30Qext*5` is `p.Ter30GlnextTer5`). Elsewhere a one-letter `X` is `Xaa`. The count after `fs` or `ext`
is a number, so `ext+5` is `ext5` and `fsTer007` is `fsTer7`. The residues an insertion, deletion-insertion or repeat
names are kept and printed as written, one letter or several (`p.Arg2delinsKL`, `p.Arg2_Lys3insX`).

**Example:**

```text
NC_TEST.1:g.10_11ins(20)       →  NC_TEST.1:g.10_11insN[20]
NC_TEST.1:g.10_12del3insN[5]   →  NC_TEST.1:g.10_12delinsN[5]
NM_W.1:c.10_12CAG(70_80)       →  NM_W.1:c.10_12CAG[(70_80)]
NP_W.1:p.R97Pfs*23             →  NP_W.1:p.Arg97ProfsTer23
NM_000051.3:c.9170_9171delGA   →  NP_000042.3:p.(Ter3057PheextTer4)     biocommons and VariantValidator write ext*4; compared as equal
```

- **Differs:** VariantValidator and biocommons write `ext*`.
- **Spec (agrees):** [insertion](https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/);
  [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/), where "both notations are
  acceptable"; [deletion](https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/);
  [protein repeated sequences](https://hgvs-nomenclature.org/stable/recommendations/protein/repeated/) for
  `p.(Gln18)[(70_80)]`; [standards](https://hgvs-nomenclature.org/stable/background/standards/) for the three-letter
  code and "use 'Xaa' only", since `X` has been used for a stop.
- **Tests:** `vrs_length_test::length_expressions_read_back_as_the_recommended_spelling`;
  `real_transcripts_test::biocommons_real_transcript_cases` (`ext*` compared equal to `extTer`);
  `printed_reparse_test::a_repeat_of_a_range_of_copies_prints_as_hgvs_writes_it`,
  `::a_repeat_of_a_range_of_copies_parses_as_hgvs_writes_it`, `::a_one_letter_frameshift_is_the_three_letter_one`,
  `::a_frameshift_to_a_stop_written_x`, `::an_extension_from_a_stop_written_x_is_the_one_written_ter`,
  `::an_n_terminal_extension_from_an_unknown_residue`, `::an_extension_count_written_plus_is_the_bare_count`,
  `::an_extension_count_with_leading_zeros_is_the_number`, `::a_frameshift_count_with_leading_zeros_is_the_number`,
  `::an_extension_count_after_x_is_the_distance_to_a_stop`;
  `statements_test::a_dict_for_an_extension_from_a_stop_written_x_loads_at_ter`,
  `::a_dict_for_an_extension_by_a_count_after_x_loads_as_the_distance_to_a_stop`,
  `::a_dict_for_an_n_terminal_extension_from_x_loads_from_an_unknown_residue`.

### A range written backwards is refused when parsed

**A range whose start is after its end is a `ParseError`, not a variant that projects to a plausible-looking genomic
range and panics when its bases are sliced.**

Transcript positions order by region (5'UTR and CDS before the 3'UTR), then base, then intronic offset. A range built in
code that runs backwards is a `ValidationError` where it is resolved. Published variants do contain them, as typos.

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
- **Spec (agrees):** [general recommendations](https://hgvs-nomenclature.org/stable/recommendations/general/), positions
  in a range are given 5' to 3'.
- **Tests:** `inverted_range_test::a_range_written_backwards_is_refused_when_parsed`,
  `::a_range_built_backwards_is_an_error_wherever_it_is_resolved`.

### Positions number from 1, and an unknown one is not 0

**`g.0`, `m.0`, `n.0`, `r.0`, `c.0`, `c.-0` and `c.*0` are a `ParseError`: a sequence numbers its bases from 1, and a
coding sequence runs `c.-1` to `c.1`, with `c.*1` the first base after the stop. A g. or m. position that is not known,
`?`, is held as unknown, never as a number.**

Were `?` held as position 0, a written 0 would read as `?`, and `g.0A>G` would be a substitution somewhere unknown
rather than an error; a transcript position 0 has no base to convert to, and lands on a neighbour's. A dict gives an
unknown g. or m. position as `"?"` and a known one as its number; a dict with position 0 is refused like the text.

A count after `fs` or `ext` starts at 1 for the same reason. `fsTer1` puts the stop at the first changed residue and
`extTer1` one past the residue that replaced the stop, so `fsTer0` and `extTer0` name a stop before either; `ext-0` puts
the new start at the old one. Each is a `ParseError`, and a dict with a count of 0 is refused.

**Example:**

```text
NC_W.1:g.0A>G               →  ParseError
NC_W.1:g.0_5del             →  ParseError
NM_W.1:c.-0A>G              →  ParseError
NM_W.1:c.*0A>G              →  ParseError
NP_W.1:p.Arg97ProfsTer0     →  ParseError
NP_W.1:p.Met1ext-0          →  ParseError
NC_W.1:g.(?_5)_(10_?)del    →  NC_W.1:g.(?_5)_(10_?)del
```

- **Spec (agrees):** [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "numbering is `g.1`,
  `g.2`, `g.3`, ..., etc., from the first to the last nucleotide", the same for `m.` and `n.`, and "there is no
  nucleotide `c.0`". [Frameshift](https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/): the stop is
  counted "starting at the first amino acid changed by the frameshift (codon 1)";
  [extension](https://hgvs-nomenclature.org/stable/recommendations/protein/extension/): `p.Met1ext-5` places a new start
  upstream, numbered like the coding sequence it sits in, which has no 0.
- *Not checked:* biocommons, VariantValidator and ClinVar.
- **Tests:** `position_zero_test::a_genomic_substitution_at_zero`, `::a_genomic_range_from_zero`,
  `::an_uncertain_genomic_bound_at_zero`, `::a_mitochondrial_substitution_at_zero`, `::a_coding_substitution_at_zero`,
  `::a_coding_substitution_at_minus_zero`, `::a_coding_substitution_at_star_zero`, `::a_coding_range_to_zero`,
  `::a_non_coding_substitution_at_zero`, `::an_rna_substitution_at_zero`,
  `::a_dict_for_a_coding_position_of_zero_is_refused`, `::a_frameshift_to_a_stop_zero_residues_on`,
  `::an_extension_to_a_stop_zero_residues_on`, `::an_extension_to_a_start_zero_residues_upstream`,
  `::an_extension_by_a_count_of_zero`, `::a_dict_for_a_frameshift_to_a_stop_zero_residues_on_is_refused`,
  `::an_unknown_genomic_bound_reparses_to_itself`, `::an_unknown_genomic_bound_survives_a_dict`,
  `::a_dict_for_an_unknown_genomic_bound_names_it_with_a_question_mark`,
  `::a_dict_for_a_genomic_position_of_zero_is_refused`.

### A number too large to hold is refused when parsed

**A position, residue number, offset, repeat count, inserted length, copy count or stop distance too large for the
integer weaver holds it in is a `ParseError` naming the number, not a value read as 0 or kept as text.**

Read as 0, an out-of-range count turned a repeat or an insertion into a different variant without a word; kept as text,
a stop distance could not be compared or projected. No real sequence comes near these limits, so a number past them is a
typo or a corrupted record, and the answer is an error. A number before a protein repeat's count is not in the grammar:
it was taken as the count's lower bound, and names neither the repeated residues nor the count. An SPDI or VRS location
is held to the same line: a deletion length or range past what an index can hold, an SPDI or exact VRS Allele position
past the end of its sequence, and a VRS coordinate no g. position can number (an end at interbase 0, one past
`g.2147483647`) are each a `ValidationError`.

**Example:**

```text
NM_W.1:c.99999999999A>G                      →  ParseError    a position is held in 32 bits, signed
NM_W.1:c.10_12CAG[99999999999]               →  ParseError
NC_W.1:g.10_11insN[99999999999999999999999]  →  ParseError
NP_W.1:p.Arg97ProfsTer4294967296             →  ParseError    a stop distance is held in 32 bits
NP_W.1:p.(Gln18)5[(6_7)]                     →  ParseError    was read as p.(Gln18)[(5_6)]
```

- **Spec (silent):** HGVS sets no bound on a number.
- *Not checked:* biocommons, VariantValidator and ClinVar.
- **Tests:** `src/parser.rs::tests::a_genomic_position_past_i32`, `::a_residue_number_past_i32`,
  `::a_transcript_position_past_i32`, `::an_intronic_offset_past_i32`, `::a_repeat_count_past_i32`,
  `::a_repeat_range_end_past_i32`, `::an_inserted_length_past_usize`, `::a_copy_count_past_i32`,
  `::a_frameshift_count_past_u32_is_a_parse_error`, `::an_extension_count_past_u32_is_a_parse_error`,
  `::a_protein_repeat_with_a_number_before_its_count_is_a_parse_error`;
  `vrs_round_trip_test::an_spdi_deletion_length_past_usize`, `::an_spdi_deletion_length_running_past_usize`,
  `::an_spdi_deletion_running_past_usize`, `::an_spdi_insertion_past_the_sequence`, `::an_spdi_insertion_past_i32`,
  `::an_spdi_insertion_past_u32`, `::an_spdi_insertion_at_the_last_index`,
  `::an_spdi_deletion_length_past_the_sequence`, `::a_vrs_deletion_ending_at_interbase_zero_is_refused`,
  `::a_vrs_deletion_whose_end_range_starts_at_interbase_zero_is_refused`, `::a_copy_number_at_the_last_u32`,
  `::a_copy_number_past_i32`, `::a_copy_number_past_u32`, `::a_deletion_range_from_the_last_index`,
  `::a_length_expression_at_the_last_index`, `::a_copy_number_at_the_last_g_position_keeps_it`.

## The data contract

### A range past the end of a sequence returns the bases that exist

**`get_seq` must return what is there rather than error or return nothing.**

The core pages through sequences in fixed blocks, and a short final block is how it learns where the sequence ends.

**Example** (a 40-base sequence):

```text
get_seq(ac, 35, 100)   →  the last five bases
get_seq(ac, 100, 110)  →  ""
```

- **Spec (silent):** none; weaver's own interface.
- **Tests:** `reference_store_agrees_with_direct_slicing` (property);
  `tests/test_refget.py::test_sequences_come_from_the_server_with_ranges_clamped`.

### Interval methods are half-open and 0-based

**HGVS positions are 1-based and inclusive, with `c.` skipping the non-existent position 0; interval methods return
0-based half-open ranges. Exon `reference_end` is inclusive.**

**Example** (`NM_PLUS10.1` and `NM_MINUS10.1` are one exon at genomic 1000..1099 with the CDS at transcript index 10):

```text
NM_PLUS10.1:c.1A>G      →  (1010, 1011)
NM_PLUS10.1:c.*1A>G     →  (1040, 1041)     the base after the stop
NM_MINUS10.1:c.1A>G     →  (1089, 1090)
NM_MINUS10.1:c.1_3del   →  (1087, 1090)     low to high on the genome
```

- **Spec (agrees; silent on 0-based half-open intervals and an inclusive exon end, which are weaver's interface):**
  [numbering](https://hgvs-nomenclature.org/stable/background/numbering/): "there is no nucleotide `c.0`".
- **Tests:** `transcript_coordinates_round_trip` and `intron_offsets_agree_from_either_exon` (properties);
  `transcript_coordinates_test::spdi_interval_honours_cds_end_anchor`, `::spdi_interval_honours_non_zero_cds_start`,
  `::spdi_interval_on_minus_strand`.

### A mapper keeps its cache

**A `VariantMapper` holds sequence blocks and refget accessions for as long as it lives. Build one and reuse it.**

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

**A sequence source need know nothing about digests; the `Refget` lookup is a separate object, and a provider without
one still works for everything but naming the sequence behind a digest.**

**Example:**

```python
table = DigestTable.from_fasta("genome.fa")            # or a RefgetProvider, or any object with the two methods
mapper = weaver.VariantMapper(provider, refget=table)
allele = mapper.to_vrs(weaver.parse("NC_TEST.1:g.7G>C"))
str(mapper.from_vrs(allele))                           # "NC_TEST.1:g.7G>C", the accession looked up from the digest

weaver.VariantMapper(provider).from_vrs(allele)                        # DataProviderError: no lookup
weaver.VariantMapper(provider).from_vrs(allele, "NC_TEST.1")           # works: the caller names the sequence
```

- **Spec (agrees; silent on keeping the digest lookup apart from the sequence source):**
  [refget](https://samtools.github.io/hts-specs/refget.html) defines the digest and the lookup of a sequence by it.
- **Tests:** `vrs_round_trip_test::the_sequence_is_named_by_the_caller_or_looked_up_and_always_checked`;
  `tests/test_refget.py::test_refget_lookups_both_ways`;
  `tests/test_digest_table.py::test_the_table_agrees_with_what_the_mapper_computes_and_names_the_sequence_back`.
