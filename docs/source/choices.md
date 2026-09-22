# The Choices Weaver Makes

HGVS leaves room for judgement, and the tools that read it have settled on different answers. This
page lists every deliberate choice weaver makes where another tool could reasonably answer
otherwise, why weaver chose as it did, and what the other tools do. "Agrees" and "differs" are
stated only where the other tool's behaviour was checked, either by running it (VariantValidator's
REST API, the biocommons `hgvs` grammar and regression tables, ClinVar's 100,000-variant sample) or
by reading its code; anything else is marked as not checked.

The tools referred to: **HGVS** is the nomenclature (hgvs-nomenclature.org); **biocommons** is the
Python `hgvs` package weaver was originally a port of; **VariantValidator** is the web service;
**ClinVar** is the descriptions in its variant summary; **VRS** is the GA4GH Variation
Representation Specification 2.0.1.

## Projection between sequences

**A projection states the target's bases.** A transcript record and its genome can differ at a
base (RefSeq transcripts are curated against submitted mRNAs). When a variant is projected, the
edit is resolved on the source to the bases it removes and the bases it puts there; where the
target holds different bases over the projected range, weaver writes the minimal edit that turns
the target's bases into that alternate. MUC2 `NM_002457.5:c.12468C>A`, over a record C and a
genomic G, is `g.1099802G>A`; the record's `c.1568=` on SHANK3, over a genomic T, is
`g.…T>C`; and the genome's `g.…T>C` projected onto that record is `c.1568=`. A duplication or
inversion carries the source's bases, because that is what the molecule becomes: `c.12468dup`
is `g.1099802delinsCC`, not a duplication of the genome's G. A deletion is positional.
*Biocommons agrees* (`replace_reference`). *VariantValidator agrees* on substitution, identity,
duplication and deletion at the MUC2 base. Where record and genome agree, which is nearly
everywhere, the edit is carried across unchanged.

**A change across a splice junction has no genomic form.** `r.44_47del` spanning an exon boundary
describes the spliced RNA. weaver refuses to project it to the genome and says why, while its
`c.` spelling still projects (to a deletion that includes the intron) and its protein is still
predicted. *VariantValidator differs*: it projects the `c.` reading with the warning "spans at
least one intron". weaver refuses because a warning is easy to miss and the two molecules are
not the same.

**No transcript-side normalisation before projecting.** weaver projects the variant where it is
written. *VariantValidator differs*: it 3'-shifts on the transcript first (`c.1568dup` becomes
`c.1569dup`) and reports the automapping. Normalisation is available separately
(`normalize_variant`) and is the caller's decision.

**Gapped exons are projected locally.** Where a transcript exon aligns to the genome with
insertions and deletions, weaver projects through the exon's CIGAR and writes a local edit.
*VariantValidator differs* on heavily gapped exons: for SHANK3's `c.1568dup` it falls back to a
2.5 kb whole-region delins from `g.50695049` to `g.50697558`. weaver's `g.50697558delinsCC` is
the deliberate answer.

**A repeat is projected as its whole run.** `c.10AC[3]` on the minus strand is widened to the full
run of the unit on the transcript before projecting, so the other strand reads it from the right
end. Without this the genomic repeat started mid-run. *Not checked against other tools.*

**A position outside every exon is an error, not an extrapolation.** Transcript positions resolve
through the exon structure only; a `c.` position that falls in no exon and has no intronic offset
is a typed error. *Biocommons* extrapolates in some paths. weaver does not, because a silently
extrapolated coordinate is wrong without saying so.

**Intronic positions are carried as given.** An intronic edit projected onto a transcript has no
transcript base to compare against, so the edit is kept as written; projected onto the genome,
the bases it states are checked against the genome like any other. Normalisation leaves intronic
edits where they are.

**`m.` is `g.` on the mitochondrial reference.** The two differ only in the letter they are
written with; normalisation, SPDI, validation and equivalence are shared. `r.` is `c.` (on a
coding transcript) or `n.` (on a non-coding one) in RNA letters, and every `r.` operation is a
conversion to that spelling and back. *HGVS agrees* on the numbering; *VariantValidator* accepts
`r.` in exons but rejects `r.*10`, intronic `r.`, `r.spl` and `r.0`, all of which HGVS allows and
weaver accepts.

## Normalisation

**3' shift, cyclic over repeats.** A deletion, duplication or insertion is shifted as far 3' as
the reference keeps matching it cyclically, on the strand of the sequence the variant is written
on. An insertion that repeats the bases immediately before it is written as a duplication, and an
insertion that slides has its bases rotated by how far it moved. *HGVS, biocommons, VariantValidator
and ClinVar agree.*

**A delins is not shifted.** Only pure deletions, duplications and insertions slide; a delins with
bases on both sides is left where it was written. Shifting it was wrong on 72 ClinVar rows.
*Biocommons agrees.*

**`del` and `dup` are written bare.** Normalised output writes `c.306del`, not `c.306delC`. Bases
the input stated survive when the edit does not move and are dropped when it does, since they
would then be stale; nothing is filled in. *HGVS and VariantValidator agree; biocommons differs*
(it fills the deleted bases in). There is no option: this is the one spelling.

**A repeat resolves to its whole run.** `c.10AC[3]` means every existing copy of the unit from
`c.10` becomes three copies, so its SPDI and protein reading cover the run. Before this, 0 of 394
ClinVar repeat rows matched on SPDI; after, 386. *ClinVar agrees* on those rows.

## Protein consequences

Protein consequences are read from codons, not from a diff of two protein strings. The rules that
follow all bear on where the stop is.

**The declared CDS end is the stop.** The reference protein ends where the transcript record
declares the CDS ends, when that codon really is a stop; only a loosely marked CDS end falls back
to the first stop in the translation. A selenocysteine `TGA` inside the CDS is therefore not a
stop: SEPN1 `NM_020451.2:c.943G>A` is `p.(Gly315Ser)`. *Biocommons differs* (`p.?`, since it sees
two in-frame stops). *Not checked* against VariantValidator.

**A frameshift whose first changed residue is the stop is an extension.** `Ter902ArgextTer92`,
not `Ter902ArgfsTer93`. *HGVS and biocommons agree*; five ClinVar rows moved from mismatch to
Analogous.

**An extension needs the stop codon itself to change.** An in-frame insertion just before the stop
that leaves the stop intact is an insertion, not an extension, even though the stop's position is
the first residue that differs.

**A stop formed inside inserted bases is a delins ending in Ter**, not a frameshift: translation
ends before the new frame reads a single reference base.

**In-frame changes are written 3'-most.** The shared tail is trimmed, so an in-frame deletion or
insertion is written at its 3'-most equivalent residues. *HGVS agrees.* *ClinVar differs* in
spelling: it writes `Xxx_Yyyins…` at the 3' end of a run where weaver writes `dup`, and one-letter
repeat forms such as `p.490PRS[1]`; 170 rows in 100,000 are the same protein in another spelling,
and equivalence judges them Analogous.

**A start-codon change is written specifically.** `c.1A>G` is `p.(Met1Val)`. *HGVS and
VariantValidator prefer `p.Met1?`*; weaver keeps the specific prediction because it carries more
information, and `VariantTransformSettings(start_codon=HgvsQuestion)` rewrites it on request.

**An edit across the CDS start.** A deletion that removes the whole CDS is `p.0?`; an edit that
starts in the 5'UTR and reaches into the CDS disrupts the start codon and is `p.Met1?`; an edit
entirely upstream is `p.?`. *ClinVar differs*: it writes `Met1fs`, `Met1_Glu2insGly…`, committing to
a consequence; weaver does not commit because the initiation site is not predictable. *Biocommons
agrees* on `p.Met1?` (its own test table). For a substitution entirely in the 5'UTR *VariantValidator
writes `p.(=)`* where weaver writes `p.?`; both are defensible, and weaver's is the more honest
since a 5'UTR change can create an upstream start.

**Frameshift length.** `fsTer N` counts to the first stop the new frame reaches; an immediate stop
is a plain nonsense substitution. *HGVS agrees.*

**Protein alleles.** From a `p.` description, only edits that name a sequence have an allele:
frameshift, extension, `p.?`, `p.0?` and `p.Met1?` are unsupported. From a coding variant, every
consequence has one (`protein_allele`, `protein_vrs`), because the edited transcript's translation
is known: a frameshift's allele is the residues from the first change to the new stop. The
translated CDS must be the protein the provider serves; a difference is an annotation error and
is reported as one. A silent change is the reference allele over the codons the edit touched; a
deleted CDS is the deletion of the whole protein.

## Validation and stated bases

**Stated bases are checked by `validate`, and nowhere else.** A variant's stated reference
(`c.123A>G`'s A) is compared with the sequence by `validate`, and by `c_to_p`, which errors on a
mismatch. Canonical alleles, SPDI and VRS use the sequence's own bases and ignore what was stated,
for nucleotides and for protein residues alike, because the allele describes the sequence, not
the description. *VRS agrees.*

## Equivalence

**Judged by allele and by the protein left behind.** Two nucleotide variants are the same change
exactly when their canonical alleles are equal. A coding variant and a protein description agree
when the protein the variant leaves is the protein the description leaves; two descriptions agree
when they leave the same protein. `p.Tyr165Ter`, `p.Ala164_Tyr165insTer` and the `c.` deletion
behind them are one change. `Identity` is the same text, or a prediction written exactly as given.
[How it decides](equivalence_logic.md) has the full rules. *Biocommons differs*: it compares
normalised text.

**A description that says nothing matches nothing.** `p.?` matches no other description.
`p.Met1?` says only where the change starts and matches a description whose change starts there;
it does not match ClinVar's `Met1fs`, a commitment weaver did not make. Two ClinVar rows judge
Different for this reason, on purpose.

**Judging with no protein sequence is an error, not a guess.** Earlier versions reconciled the
residues two descriptions happened to name; now the sequence is required.

**Versions of one protein accession are compared on ours.** ClinVar often names an older `NP_`
version; the descriptions are read against the protein the provider serves.

**A cis allele compares as a set.** `c.[a;b]` and `c.[b;a]` are the same change; a single-member
`c.[a]` compares as `a`; a cis allele never compares equal to a single variant of two or more
members.

## Alleles, SPDI and VRS

**The canonical allele is fully justified.** A change is trimmed to what it alters and then
widened over the whole region in which it could equally be written, so every spelling of one
change is one allele with one identifier. *VRS and VOCA agree*; this is their normalisation.
`to_spdi_unambiguous` renders it. The plain `to_spdi` renders the 3'-normalised variant instead,
with an insertion placed at its second flank as SPDI counts.

**Refget accessions are computed over the normalised sequence.** Letters uppercased, everything
else dropped, as the refget specification defines; chr19 of the NCBI FASTA, uppercased, hashes to
the specification's published `SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl`. Without a `Refget` lookup the
accession is computed from the whole sequence; `from_vrs` then needs the accession passed. The
digest is checked against the sequence either way.

**Uncertain breakpoints.** A deletion with uncertain breakpoints, `g.(?_100)_(200_?)del`, renders
as an Allele with Range bounds and an empty state, unnormalised; a duplication with uncertain
breakpoints has no Allele (its bases are unknown) and renders as a `CopyNumberChange` gain. An
unknown bound prints as `?`; a parenthesised exact position is exact. *VRS agrees* on both
renderings.

**`copyChange` is a label.** VRS 2.0.1 defines `copyChange` as a string enum (`gain`, `loss`, …);
2.0.0 used EFO codes in a `MappableConcept`. weaver writes labels and reads labels, bare EFO
CURIEs and the 2.0.0 object form. The specification's published `CopyNumberChange` example digest
only reproduces when `copyChange` is digested as the pre-release string `EFO:0030071`; both digests
are pinned in tests so the machinery is anchored to a published value while output follows the
current schema. EFO's high-level loss is `EFO:0020073`, not in the `00300xx` block.

**`CisPhasedBlock` members are sorted before digesting**, as the digest-serialisation rule requires,
so the identifier does not depend on the order written; the `members` array keeps that order. The
schema's `minItems: 2` is relaxed to one because HGVS allows the degenerate `c.[145C>T]`. The
in-trans form `[..];[..]` is rejected: it describes two molecules.

**Reading back.** `from_vrs` and `from_spdi` give the 3'-normalised variant on the allele's own
sequence. An insertion trimmed to before the first base becomes a delins of that base, since HGVS
has no insertion before base 1; a replacement that is the reverse complement of what it replaces
reads back as `inv`; a `LengthExpression` comes back as `insN[n]`; a Range of copies is refused,
since HGVS has no spelling for it.

**What is not rendered.** `Adjacency`, `Terminus` and `DerivativeMolecule` describe breakends and
fusions this grammar has no input for.

## Parsing

**The grammar is checked against biocommons's table**, 580 inputs over the 92 rules the two
grammars share, and agrees on all but one: weaver accepts a terminator inside an amino acid
sequence (`insTerGlu`), which *biocommons rejects* and *ClinVar writes*.

**Forms accepted beyond biocommons.** Alleles in cis for every coordinate system; `r.(123a>g)`
predicted RNA changes and the statements `r.0`, `r.?`, `r.spl`, `r.=`; `p.0`, `p.0?`, `p.Met1?`;
uncertain genomic breakpoints with `?`; insertions of a stated length, `insN[20]`, `insN[(20_30)]`,
`delinsN[12]` and the older `ins(20)`; `copyN`; `r.` positions with `*`.

**Recommended spellings on output.** `insN[20]` rather than `ins(20)`; `extTer8` rather than `ext*8`
(*VariantValidator writes `ext*`*; the two are compared as equal); bare `del` and `dup`.

## The data contract

**A range past the end of a sequence returns the bases that exist.** `get_seq` must return what is
there rather than error or return nothing, because the core pages through sequences in fixed
blocks and a short final block is how it learns where the sequence ends.

**Interval methods are half-open and 0-based**; HGVS positions are 1-based and inclusive, with
`c.` skipping the non-existent position 0. Exon `reference_end` is inclusive.

**A mapper keeps its cache.** A `VariantMapper` holds sequence blocks and refget accessions for as
long as it lives. Build one and reuse it.

**Refget is its own seam.** A sequence source need know nothing about digests; the `Refget` lookup
is a separate object, and a provider without one still works for everything but naming the
sequence behind a digest.
