# hgvs-weaver Code Review

> Prepared 2026-02-24. Covers all Rust source files in `hgvs-weaver/src/`.
> Updated 2026-02-25 to reflect fixes applied in subsequent commits.

---

## Open findings

| # | Severity | File | Description |
|---|----------|------|-------------|
| 12 | Minor | `structs.rs` | `SimplePosition.end` for uncertain positions not yet parsed/formatted |

---

## Detailed findings

### 5. `DataProvider::c_to_g` interface contract mismatch (Major) — resolved

`DataProvider::c_to_g` has been removed. `TranscriptMapper` now owns c./n. position
resolution (`position_to_n`, `position_to_g`, `interval_to_n`, `interval_to_g`), and
`BaseOffsetInterval::spdi_interval`, `VariantEquivalence`, `AltSeqBuilder` and
`VariantMapper::get_c_indices` all route through it. Regression coverage lives in
`hgvs-weaver/tests/transcript_coordinates_test.rs`.

---

### 12. `SimplePosition.end` not yet parsed/formatted (Minor)

**File:** `hgvs-weaver/src/structs.rs`

`SimplePosition.end` represents the uncertainty window for positions like `(1_3)_(7_10)`, where the exact breakpoint is unknown. The field is intentional and correctly models this HGVS notation.

**Remaining work:** The parser and formatter do not yet round-trip uncertain genomic positions. Until they do, `end` will remain `None` in practice.

---

## Type system opportunities

### A. Replace `strand: i32` with a `Strand` enum

`Exon::alt_strand` and `Transcript::strand()` return `i32` where only `1` and `-1` are meaningful. Every consumer checks `== 1` or `== -1` with no exhaustiveness guarantee. A `Strand { Plus, Minus }` enum would make invalid strand values unrepresentable and enable exhaustive matching.

### C. Anchor consistency enforcement

`BaseOffsetPosition` can be constructed with any `Anchor` regardless of the enclosing variant type (e.g., a `NVariant` with `Anchor::CdsStart`). Consider a builder or `From` impl that enforces the correct anchor per coordinate system.

---

## Test coverage gaps

The following areas have limited or no test coverage:

| File | What is still untested |
|------|------------------------|
| `transcript_mapper.rs` | `n_to_c`, `c_to_n` |
| `altseq.rs` | `AltSeqBuilder::build_altseq` for most edit types |
| `altseq_to_hgvsp.rs` | Most protein consequence cases (frameshift, delins, synonymous) |

`transcript_mapper.rs` now has unit tests for `g_to_n` (exonic, intronic, CIGAR, and minus-strand ordering). `mapper.rs` is covered by `mapping_test.rs` and `test_shift.rs`, including multi-base insertion shifting.
