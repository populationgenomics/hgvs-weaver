# Glossary

**Related:** [`PRODUCT.md`](../PRODUCT.md), what weaver is and the principles it keeps; [`choices.md`](choices.md),
every choice weaver makes where another tool could answer otherwise; [`equivalence_logic.md`](equivalence_logic.md), how
two descriptions are judged.

weaver reads HGVS descriptions of changes to a sequence, checks them against the sequences they cite, moves them between
a transcript, its protein and the genome, and decides whether two of them name the same change. HGVS, VRS and the tools
that read them each have their own words, and several everyday words — reference, identity, statement, normalisation —
mean different things in different places. This page fixes one word for each concept weaver's code, docs and errors talk
about. Where another word is in common use for the same thing, it is listed under _Avoid_: not because it is wrong
elsewhere, but so that one concept has one name here.

The terms are grouped by area, and within an area build on each other in order.

## Sequences and records

**Record**: A sequence published under an accession: a genome, a transcript or a protein, as its publisher wrote it.
_Avoid_: entry (for a sequence)

**Reference**: The sequence a description is read on: the genome for `g.` and `m.`, the transcript record for `c.`, `n.`
and `r.`, and the protein record for `p.`. An intronic position is read on the genome, on the transcript's strand.
_Avoid_: reference sequence (except when quoting HGVS), "reference" for the genome specifically

**Genome**: The chromosome or contig record a transcript is placed on, and the reference of a `g.` description. _Avoid_:
reference (when the genome specifically is meant), assembly (for a single record)

**Transcript record**: A transcript's own published sequence. It can differ from the genome where the transcript is
placed, by a base or by a gap. _Avoid_: transcript (when the sequence specifically is meant), mRNA

**Protein record**: A protein's own published sequence.

**Transcript model**: What a provider supplies about a transcript record besides its sequence: its accession and gene,
its coding sequence, and its placement where it has one. _Avoid_: transcript data, annotation

**Placement**: Where a transcript lies on a genome: the genome's accession, the strand, and how each exon aligns to it.
Where an edit sits is its position, not its placement. _Avoid_: alignment (for the whole), mapping

**Placed record**, **unplaced record**: A transcript record with a placement, or one no assembly places, such as a
GenBank mRNA. An unplaced record is read on its own sequence alone, and anything that needs the genome is refused.

**Strand**: Which strand of the genome a transcript is read from, plus or minus.

**Exon**: A stretch of a transcript record and the stretch of the genome it aligns to.

**Exon alignment**: How an exon's bases pair with the genome's, base by base, including where one sequence holds bases
the other lacks.

**Gap**, **unaligned run**: A run of bases inside an exon that one sequence holds and the other lacks. A position on it
has no single counterpart on the other sequence, so it projects as the whole run. _Avoid_: indel (for an alignment gap)

**Soft-clipped end**: An exon end whose bases the alignment leaves unpaired with the genome. A position on it has no
counterpart and is not projected.

**Intron**: The stretch of genome between two exons of a placed transcript, which the transcript record does not hold.

**Exon border**: The exon base an intronic offset counts from: an exon's last base for `+`, its first for `-`. _Avoid_:
exon boundary

**Junction**: The point between two neighbouring bases, where an insertion goes.

**Exon junction**: The junction on a transcript record between the last base of one exon and the first of the next, with
an intron between them on the genome. _Avoid_: splice junction, exon/exon boundary

**Transcript**: An RNA a gene gives, as a transcript model and a transcript record describe it.

**Coding sequence**, **CDS**: The stretch of a transcript record that is translated, from the first base of its first
codon to the last base of its stop codon, or to the end of the record where it runs off. _Avoid_: coding region, ORF

**Start codon**, **initiator**: The first codon of a coding sequence when translation starts there. Its residue is Met
whatever its bases.

**Internal CDS start**: A coding sequence numbered from a codon inside a protein, where translation does not start:
residue 1 is that codon's own residue. MAVE targets are written this way.

**Open CDS end**: An end of a coding sequence that runs off the transcript record: the record does not carry the start
codon (open at the 5' end) or the stop codon (open at the 3' end). No position is numbered from a missing end, and a
protein consequence that needs bases past it is refused. Only a coding sequence is open.

**Reference stop**: The stop codon that ends the reference protein: the codon the transcript model names as the CDS's
last, where that is a stop, and otherwise the first stop the translation reaches. _Avoid_: terminator, declared stop

**Selenocysteine**, **Sec**: A residue encoded by an in-frame TGA before the reference stop, read as Sec, not as a stop.

**5' UTR**, **3' UTR**: The parts of a transcript record before and after its coding sequence.

**Read-through**: Translation that carries on past the reference stop into the 3' UTR, or off the end of the record,
because the stop was changed.

**Reference protein**: The protein a coding sequence encodes, residue by residue, with its stop after the last residue;
or a protein record read the same way. A coding sequence open at its 5' end has none.

## Positions

**Coordinate system**: What a description's positions count along: `g.` a genome, `m.` the mitochondrial genome, `c.` a
coding transcript, `n.` a non-coding transcript, `r.` a transcript's RNA, `p.` a protein. _Avoid_: coordinate type
(except as the API's name for it), transcript space

**Coding** (`c.`): The coordinate system of a transcript with a coding sequence, numbered from the CDS's first base,
with `-` before it and `*` after its last. _Avoid_: cDNA

**Position**: A place on a sequence, in a coordinate system. A position in one coordinate system is never a position in
another without a projection.

**Anchor**: What a transcript position counts from: the transcript's first base, the CDS's first base, or the CDS's last
base.

**Intronic offset**: The signed distance from an exon border into an intron: `c.88+2` is two bases into the intron after
`c.88`.

**Interbase**: The convention that numbers the points between bases, so a range is half-open, `[start, end)`, and an
insertion has a length-zero range. VRS calls it inter-residue. _Avoid_: inter-residue (except when quoting VRS),
half-open (as a name for the convention)

**Uncertain**: Not known exactly: a position given as a range of possible positions, `g.(3_4)`, or a breakpoint given as
`?`. An uncertain position names no single change. _Avoid_: imprecise, inexact, fuzzy

**Footprint**: What a position occupies on the other sequence of a projection: one base; or, at a gap, an empty point
together with the run of bases the other sequence lacks there.

## Descriptions and edits

**Description**: What an HGVS string says: an accession, a coordinate system, a position and an edit, as text or as the
value parsed from it. _Avoid_: variant (in prose; the API's type is `Variant`), expression (except for the HGVS a VRS
Allele carries)

**Change**: What a description does to its reference. One change has many spellings. _Avoid_: event, mutation, variant

**Spelling**: One written form of a change. `c.5del` and `c.4_5delinsA` can be two spellings of one change. _Avoid_:
notation, form

**Edit**: The operation a description names at its position: a substitution, deletion, insertion, duplication, delins,
inversion, conversion, repeat, copy number or identity. _Avoid_: mutation type

**Delins**: An edit that replaces a range with new bases. _Avoid_: indel, deletion-insertion

**Repeat**, **repeat unit**, **run**: A repeat, `CAG[12]`, says the whole tandem run of its unit at that position
becomes the given number of copies.

**Identity**: The `=` edit: the bases at the position are as the reference has them, `c.10_12=`. An identity states the
bases it names. _Avoid_: no change

**Silent change**: A coding change that leaves its protein as it was; its protein consequence is an identity at the
residues it spans, `p.(Ala4=)`. _Avoid_: synonymous change (except when quoting)

**Whole identity**: `c.=`, `r.=`, `p.=`: the whole sequence was analysed and found unchanged.

**Statement**: A description that says what the sequence is, or how it is affected, without changing its bases. Its
opposite is a change. A statement either names bases, as an identity does and cancelling members do, or names none:
`c.=`, `p.(=)`, `p.?`, `p.0?`, `p.(Met1?)`, `r.spl`, `r.0`, `r.(?)`. A statement that names no bases has no sequence to
compare, so within an allele it is compared by its spelling.

**Written bases**: The reference bases, repeat unit or residues a description writes, such as the `C` of `c.306delC`.
They are checked against the reference. _Avoid_: stated bases (for these)

**Written bare**: Written without the reference bases an edit could carry: `c.306del`, not `c.306delC`.

**Stated bases**: The bases a description asserts are unchanged, through an identity or through members of an allele
that cancel.

**Located edit**: An edit positioned on a transcript record for prediction: either in an intron, or over a range of the
record with its new bases, spelled out or known only by their count, and the exon junctions the range spans.

**Prediction**: A consequence weaver infers rather than one that was observed: the `p.` description it derives from a
coding change is a prediction, in parentheses where HGVS brackets one (`p.?` and `p.0?` stay bare).

**Protein consequence**: The `p.` description a coding change predicts, read from the codons of the edited transcript.
_Avoid_: outcome, effect

**Frameshift**: A change after which translation reads another frame, written `fs` with the distance to the first new
stop.

**Extension**: A change that lengthens a protein, written `ext`: past a lost stop at its C-terminus, or before Met1 at
its N-terminus. _Avoid_: stop loss (for the description)

**Amino acid code**: The one-letter or three-letter names of residues. A `p.` description writes one code throughout.

**Stop**, **Ter**: The end of translation, written `Ter` or `*`. _Avoid_: terminator, `X` (for a stop)

**Ambiguous residue**: An `X` that could be read as a stop or as any residue. It is refused unless the caller says how
to read it.

**Trans**: Changes on two different molecules, `[..];[..]`. weaver refuses them as not one description.

## Alleles

**Allele**: What one molecule carries, a single change or several in cis, judged by the sequence it makes and the bases
it states.

**Cis allele**: Several changes on one molecule, written `c.[a;b]`. _Avoid_: allele in cis (in running prose), phased
allele, haplotype

**Member**: One change of a cis allele.

**Member that cannot be applied**: A member no sequence can be built from: a copy number, an uncertain position, a
conversion, an insertion of unwritten bases, a statement about the whole sequence. _Avoid_: unapplied member

**Extent**: The stretch of the reference a member's edit touches; an insertion's extent is the empty point at its
junction.

**Members that meet**: Members whose extents touch or overlap.

**Cancelling members**: Members that together leave the reference as it was. They state the bases they span.

**Reading**: One way of applying a description's members to the reference, giving one sequence. To read a description on
a sequence is to apply it there. _Avoid_: interpretation; "reading" for the reading frame or a translation

**Contradictory description**: A description with no reading. It is refused, saying why.

**Ambiguous description**: A description with several readings that make different sequences or state different bases.
It is refused, listing the readings; weaver never picks.

**Canonical allele**: A change trimmed to what it alters on its genome, or for a `p.` description its protein, then
widened over its region of ambiguity, so that every spelling of one change has one canonical allele. VRS and VOCA call
this fully justified. _Avoid_: fully justified normalisation, normalised allele

**Region of ambiguity**: The span over which an insertion or deletion could equally be written.

**Reference allele**: An allele that leaves the reference unchanged.

**Protein allele**: A cis allele of `p.` changes, `p.[a;b]`. The VRS Allele of a protein change is not a protein allele;
it is that change's VRS Allele.

## Operations

**Parse**: Read an HGVS string into a description, refusing what the grammar does not allow.

**Validate**: Check a description against its reference: its positions exist and its written bases are there. _Avoid_:
verify

**Project**, **projection**: Rewrite a description onto another sequence that numbers the same bases: a transcript onto
the genome it is placed on, or the genome onto a transcript. A coding change's protein consequence is predicted, not
projected. _Avoid_: map, mapping, transform, convert, lift over

**Respell**: Write an `r.` description as the `c.` or `n.` description of the same positions in DNA letters, or back.

**Correction**: What a statement becomes when it is projected onto a sequence that differs from its own there: the edit
that turns the target's bases into the ones the statement asserts.

**Normalise**, **normalisation**: Rewrite a description into HGVS's preferred spelling of its change: shifted 3',
written bare, a delins reduced to the edit it is, an insertion of a copy written as a duplication. _Avoid_: normalize
(in prose; the API spells it so), canonicalise; "normalisation" for the canonical allele or for spelling comparison

**3' rule**, **shift**: Moving an edit to the 3'-most position at which it writes the same change. _Avoid_: roll, slide,
justify

**Reduce**: Drop the bases a delins replaces with themselves at either end, then write what is left as the edit it is.

**Consistent variants**: For a description on a gene symbol or an unversioned accession, the description on each record
the store holds under that name that validation finds it consistent with; a versioned accession gives the description
itself.

**Translation**: Reading a coding sequence codon by codon into residues.

## Equivalence

**Equivalence**: Whether two descriptions say the same thing about a sequence: they make the same sequence and state the
same bases, or leave the same protein.

**Verdict**: An equivalence answer. **Identity**: the two are the same after spelling comparison, or one projects
exactly onto the other, or one is a coding change and the other its prediction. **Analogous**: written differently, but
equivalent. **Different**: both understood, and not equivalent. **Unknown**: weaver does not return it, refusing a pair
it cannot judge with an error; the validation run records such a pair as Unknown. _Avoid_: level, "identity" (lower
case, for the verdict)

**Spelling comparison**: The first equivalence step: comparing two descriptions as text, ignoring what is only spelling
— a prediction's parentheses, the amino acid code, a one-member allele's brackets, the order of members. _Avoid_: basic
normalisation, spelling normalisation

**Protein left**: The protein a `p.` description implies once applied to its reference; two `p.` descriptions on one
protein are compared by it.

**Unknown length**: A protein whose residues are known up to a point and then continue for an unknown number. _Avoid_:
open (for a protein)

**Paralogue**: A transcript of another protein that carries an identical copy of the sequence; it adds nothing to
equivalence.

## Data

**Store**: Anything holding the records, transcript models and catalogue weaver reads: hgvs-weaver-data's built dataset,
a directory of RefSeq files, a refget server, a test's tables. A store is read through a provider. _Avoid_: dataset,
database, data source

**Provider**: The adapter that answers weaver's questions about records and transcript models from a store. _Avoid_:
data source, backend, HDP

**Sequence cache**: The copies of sequence a provider has served that weaver keeps while it works. _Avoid_: store (for
the cache)

**Accession**: The identifier of a record, such as `NM_000051.3`. A versioned accession names one record; an unversioned
one names every version of it. _Avoid_: identifier (for an accession specifically), ID

**Identifier kind**: What an identifier names: a genome, a transcript, a protein, or a gene symbol. _Avoid_: identifier
type, type

**Wrong reference kind**: A description on an accession of the wrong kind for its coordinate system, such as `p.` on a
transcript.

**Catalogue**: What a store records about accessions: the versions held under an unversioned accession, and each
version's status and tags.

**Accession status**: Whether a version is **current**, **superseded** by a later one, or **suppressed**, withdrawn by
its publisher.

**Historical version**: A version that is not current: superseded or suppressed. _Avoid_: retired, old

**Accession tag**: A curated label on a transcript, such as MANE Select or RefSeq Select. A protein carries the tags of
the transcripts that encode it.

**Shard**, **historical shard**: One release's records in a built store; the historical shard holds NCBI's historical
versions beneath the current release's.

## VRS and identifiers

**SPDI**: `accession:position:deletion:insertion`, interbase. weaver writes it 3'-normalised, or over the canonical
allele as the unambiguous form.

**VRS Allele**: A GA4GH VRS value for one change at one location, with its state. Always capitalised, with "VRS", to
keep it apart from an allele.

**State**: What a VRS Allele says the location holds: literal bases, or a length expression over a repeat.

**CisPhasedBlock**: VRS's form of a cis allele: a set of VRS Alleles on one molecule.

**Copy number count**, **copy number change**: VRS's forms for a copy number edit: an absolute count, or a gain or loss.

**Refget accession**: `SQ.` followed by the digest of a sequence, naming it by its content. _Avoid_: digest (for the
accession)

**Computed identifier**: The `ga4gh:` identifier of a VRS value, derived from its content.

## Checking weaver

**Choice**: A decision recorded in [`choices.md`](choices.md) where another tool could reasonably answer otherwise, with
its example, the specification it rests on, its tests, and which tools agree.

**Spec verdict**: A choice's standing against the specification: it **agrees**, the specification is **silent**, or it
**differs**, and why.

**biocommons hgvs**: The Python `hgvs` package weaver descends from, used as a comparison. _Avoid_: ref-hgvs, the
reference implementation, the standard implementation

**Validation run**: The run that compares weaver's answers with ClinVar's and biocommons hgvs's over a sample of
ClinVar's descriptions. _Avoid_: harness

**Gate**: The 100,000-description validation run every change to the core is checked with before it merges: its rows are
compared with a baseline's, and each moved row is checked against the truth columns.

**Baseline**: A validation run's output on the code a change starts from.

**Moved row**: A ClinVar description whose answer differs between a baseline and a validation run.

**Truth columns**: The answers in a validation run that come from outside weaver, ClinVar's protein and SPDI, against
which a moved row is checked.
