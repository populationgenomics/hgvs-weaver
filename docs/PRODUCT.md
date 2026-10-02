# hgvs-weaver — product north star

## What it is

weaver is an HGVS engine. It parses a variant description, checks it against the reference it cites, projects it
between transcript, protein and genome, normalises it, renders it as SPDI and GA4GH VRS, and decides whether two
descriptions name the same change. The core is Rust; the Python package wraps it with typed exceptions and a
protocol for data access, so it needs no database.

```text
"NM_000051.3:c.123A>G" ──parse──▶ Variant ──project / normalise / judge──▶ Variant, SPDI, VRS, verdict
                                      ▲
                               DataProvider: transcript models and sequences
                               (hgvs-weaver-data, a GFF3 and FASTA, a refget server, …)
```

## Why these principles

1. **Positions are types, not integers.** A transcript position, a genomic position and a protein position are
   different types, each 0-based internally and 1-based in the HGVS it writes, and a transcript position carries its
   anchor and intronic offset with it. Mixing systems is a compile error in Rust, not an off-by-one at run time. The
   price is a conversion at every boundary; the return is that an entire class of defect cannot be written.
1. **Consequences are read from the sequence, not from a string.** A protein consequence comes from the codons the
   edited transcript encodes, not from a diff of two protein strings; a canonical allele is the change fully justified
   over its region of ambiguity, so every spelling of one change is one allele with one identifier. What two
   descriptions mean is settled by applying both to the reference and comparing what is left.
1. **A choice is written down, with who agrees.** HGVS leaves room for judgement and the tools that read it have
   settled on different answers. Every deliberate choice weaver makes where another tool could reasonably answer
   otherwise is recorded in [`source/choices.md`](source/choices.md) with a worked example, the specification it rests
   on, the tests that assert it, and whether biocommons, VariantValidator and ClinVar agree. A consistent decision is
   one people can file bugs against; an unrecorded one is a surprise.
1. **An answer or an error, never a guess.** A position outside every exon, a stated base the record does not carry, a
   change across a splice junction asked for a genomic form: each is an error naming why, not a nearest answer. The
   data provider is held to the same line, and weaver does not fall back from one transcript version to another or from
   a record to the genome.

## Shape

- **Core** (`hgvs-weaver/`, the `hgvs-weaver` crate): grammar, position types, `TranscriptMapper`, `VariantMapper`,
  normalisation, protein consequences, canonical alleles, equivalence. Everything the engine decides lives here and is
  tested here.
- **Binding** (the root crate, `weaver._weaver`): pyo3 over the core, with the Python-visible surface declared in
  `weaver/_weaver.pyi`. The Python package adds the `DataProvider`, `Refget` and `TranscriptSearch` protocols and the
  commands.
- **Commands**: `weaver-validate`, the comparison of weaver's answers against ClinVar's and the biocommons `hgvs`
  package's over a store; `weaver-gate`, which builds that store from NCBI's files and runs the comparison against a
  baseline. The 100,000-variant run is the regression gate for any change to the core.
- **Languages**: Rust for the engine, Python for the binding, the protocols and the commands.

## Relationship to hgvs-weaver-data

Reference data is another package's concern.
[hgvs-weaver-data](https://github.com/populationgenomics/weaver-data-provider) builds a store from NCBI's published
alignments and serves it as weaver's `DataProvider`; it depends on weaver, and weaver's runtime never imports it. The
validation commands do, through an extra and a dependency group, which is the one place the arrow is allowed to point
the other way. If changes to the protocol between the two keep needing to land in both repositories at once, they
belong in one.

## Non-goals

- Not a data source. weaver carries no transcript models and no sequences; everything it knows about a transcript comes
  through the `DataProvider` it is given.
- Not a resolution policy. Which transcript version to try when a cited one is absent, or how to read a legacy
  numbering, is decided by the application using weaver.
- Not a liftover. A projection goes between a transcript and the sequence it is placed on; relating two assemblies, or
  two placements of one transcript, is outside the engine.
- Not a clinical interpreter. weaver says what a description means and whether two agree; what that means for a
  patient is not its question.

## Scope

HGVS `g.`, `m.`, `c.`, `n.`, `r.` and `p.` descriptions; SPDI and GA4GH VRS 2.0 alleles; projection through RefSeq and
Ensembl transcript models on GRCh38 and GRCh37, as a data provider supplies them.
