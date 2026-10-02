# Writing tests

What a test should assert, what it may depend on, how it reads, and what its data may contain — independent of language.
The mechanics of the test runner live with the language: pytest in [`python.md`](python.md#pytest). Much of it is
distilled from Google's [Testing on the Toilet](https://testing.googleblog.com/) series; the rule headings track the
post titles closely enough to search.

## A test earns its place by the regression it catches

Before a test is kept, name the regression — a change in behaviour some caller could observe — that would make it fail.
If no such change exists, the test is decoration; if the only change that fails it moves code around without changing
what a caller sees, the test is pinning structure ([below](#test-at-the-surface-a-caller-uses)); if the regression it
catches is one an existing test already catches, it is duplicate coverage. The common ways a test that looks like
coverage is not:

- Logic lands with no test that reaches it.
- A test runs the new code path but its assertion would pass with the old behaviour too.
- A test asserts something so weak — the result is not `None`, the call did not raise when the raise is not the
  regression — that the regression it nominally guards passes through it.

Two procedures turn that question into evidence. Writing the test before the change and watching it fail proves it can
fail. Refactoring tests inverts it: break the production code first, confirm the failures are the assertions you expect,
and refactor the tests while they are red — one that starts passing mid-refactor was broken by the refactor.

Coverage shows what is untested. It does not show that covered code is tested: a division is covered by any divisor, a
short-circuit condition by its first operand alone. Read uncovered lines for gaps; never write a test to move the
number, or to catch a change no caller could observe — arithmetic in a log line, a tuning constant. If coverage is gated
at all, gate the lines a change touched, not a project-wide target.

## Test at the surface a caller uses

The unit under test is the surface a caller uses, not the functions behind it. Here there are three: the Rust crate's
public API — `parse` and `format`, `TranscriptMapper`, `VariantMapper` over a `DataProvider`; the Python package, the
same operations through the extension and a provider implementing the protocol; and the commands, `weaver-validate` and
`weaver-gate`. The common failure is to decompose the mapper into helpers and write a test per helper: many tests and no
evidence a variant projects correctly, for two reasons.

The first is that a helper's test has no ground truth but the author's reading of the helper. The expected value is
derived by re-running the helper's logic in the author's head, so the test and the code are one belief written twice. At
the surface, the expectation is what the user would see, which the author can derive from the task rather than from the
code.

The second is that tests of helpers pin the decomposition. A refactor that leaves every caller-visible behaviour alone
should leave every assertion alone — setup may change, assertions should not. A refactor that rewrites half the suite
shows the suite was asserting structure. A helper with a single caller, every path of which is reachable through that
caller, needs no tests of its own; test it through the caller. A behaviour shown to work only through the helper carries
false confidence, because the same path may be wired wrongly at the surface.

A helper earns its own tests when the surface cannot reach it cheaply — a parser with many cases, a numeric routine
whose branches need a hundred inputs to cover — or when more than one caller shares it. Even then, test only the inputs
a caller can produce. The tells of a suite that has crossed the line: a private function tested that a surface test
already reaches or that has one caller, tests that mock sibling helpers, and a test module in which nothing invokes the
surface.

A validation rule is the other place tests multiply — one test per rule. Before writing either, ask whether the invalid
value needs to be expressible at all ([`general.md`](general.md#make-the-invalid-state-unrepresentable)).

## Assert invariants, not the current shape

A test that pins today's values — an error's exact wording, a whole formatted transcript model — is a *change detector*:
it fails on every legitimate change and catches no defect. Name the property that has to hold for any valid state. Where
production code already enforces that property, test the enforcement — hand it a violating input, assert it raises —
rather than restating its current output: an assertion that production code makes unreachable can never fail. Swapping
in the invariant replaces the assertion, not the test.

```python
# Good — the property: projecting to the genome and back lands on the variant parsed
def test_a_coding_variant_round_trips_through_the_genome(provider: weaver.DataProvider) -> None:
    ...
    assert mapper.g_to_c(mapper.c_to_g(variant), variant.ac).format() == variant.format()

# Good — the enforcement can regress; the property it guards cannot
def test_a_reference_base_the_record_does_not_carry_is_refused(provider: weaver.DataProvider) -> None:
    ...
    with pytest.raises(weaver.TranscriptMismatchError, match='expected A, found G'):
        mapper.c_to_g(weaver.parse('NM_000001.1:c.5A>T'))

# Bad — change detector: every change to the error's wording edits this string
def test_mismatch_message(provider: weaver.DataProvider) -> None:
    assert str(error) == 'TranscriptMismatchError: expected A, found G at transcript indices 4..5'
```

- **Narrow assertions.** Assert the fields the behaviour is about, not equality of the whole object, which silently
  tests every other field too: a field added to a transcript model fails a test about CDS bounds. At most one
  whole-object equality per type, for the common case.
- **A pinned value is right where a policy freezes it.** A choice recorded in [`choices.md`](../source/choices.md) — the
  3' rule, a deletion written bare — is a decision, so the test that pins its output is not a change detector: a moved
  value is the defect, and the page names the test that would catch it. The test is whether a valid future state can
  change the value without a decision; if it can, assert the property instead.
- **A floating-point number that came through arithmetic** is compared with a tolerance (`pytest.approx`), never for
  equality: rounding error fails correct code. A number parsed straight from a fixture compares exactly.
- **Results, not interaction sequences.** A test that asserts the exact calls made on mocked collaborators is the same
  trap: it restates the implementation and passes when the behaviour is wrong
  ([below](#test-doubles-in-order-of-fidelity) for the exception).

## One behaviour per test, named for it

A test covers one behaviour, not one method. A method that answers a lookup and also caches the answer has two
behaviours, and each gets a test; a test named after the method is the warning sign. The other tell is a second call to
the system under test after an assertion — a second scenario, whose failure the first now hides. The exception is when
the sequence *is* the behaviour: a cache's eviction order, a variant normalised and then projected.

The name states the scenario and the expected outcome, so the name alone says what broke and the list of names
enumerates what the unit does. A name that has to say two outcomes is two tests.

A failure is actionable when the name and the message let someone start investigating without adding a print and
rerunning. Use the assertion that carries the detail: in pytest, `assert a == b` gets the rewritten diff and
`assert result.ok` gets "assert False"; `pytest.raises` takes a `match`, and a bare one accepts any error of that type.
Where a test needs a precondition — a file written, a store opened — assert it first with a message naming what was
unavailable, so a failed precondition is not misreported as the behaviour failing.

## A test reads without its reader leaving it

A test has no tests of its own. A reader can only verify it by inspection, and inspection needs everything the assertion
depends on in view. Where production code is DRY, test code leans the other way: repeat a line rather than hide it
behind shared state, a helper that performs the action, or a loop over assertions — a `parametrize` row is a literal in
and a literal out that fails on its own; a loop inside one test is neither.

- **No logic in tests.** Inputs and expected outputs are literals; a test that computes its expectation can share the
  production bug it exists to catch. A fixture-shipped constant asserted on by name is a literal; a computed expectation
  is not. Where logic is unavoidable it moves to a test helper, and a non-trivial test helper gets tests of its own.
- **What the assertion depends on is set or named in the test.** Fixtures and builders carry irrelevant construction;
  the value an assertion depends on is passed in the test body — `transcript(cds_start=CDS_START)` then assert on
  `CDS_START` — never left to a fixture two hundred lines away or a builder's default, even when the default is right.
- **A builder, not a helper per combination.** A zero-argument builder returns an object with the logically required
  fields set; each test overrides the fields it is about. A positional helper that grows an argument per test, or a
  helper per combination, makes every irrelevant field a change to every test. In Rust: a baseline struct and
  `..Default::default()` or a builder; in Python: a baseline dataclass or dict plus `dataclasses.replace` or a keyword
  factory.
- **Values that could expose the bug.** Never a type's default — a missing int reads as zero, so `0` cannot tell stored
  from absent. A distinct value for every argument, so a swapped or reused argument is caught. Then the empty case, the
  boundaries — an exon's first and last base, a CDS end, a minus-strand transcript, an intronic offset of one — and the
  input that takes the complex path: a gapped alignment, a repeat that crosses an exon boundary.
- **Parametrise over inputs to one behaviour, never over behaviours.** A table that grows a dimension per behaviour
  hides which case failed and can carry a bug in its own indexing that the code under test shares.
- **Brittleness is a missing assertion.** A test that depends on an irrelevant detail — set iteration order, an error's
  exact text — wants an assertion that names the salient property: unordered equality, a predicate, a match on the part
  of the message that matters.

## Deterministic by construction

A test that passes and fails on the same code is a defect worse than a reliable failure: its failures are discounted and
its passes doubted. Root-cause it. There is no automatic rerun here, and a flaky test is quarantined only with an issue
filed against it. Hermeticity ([below](#tests-run-with-no-cloud-and-no-network)) removes the largest source; what
remains is in the test's own hands:

- **Sleep is not synchronisation.** Wait on the event or condition itself, with a timeout. The timeout bounds a failing
  run and never slows a passing one; a fixed delay does the reverse, and grows until it fails again. Most code under
  test is single-threaded — test that separately, and reserve thread choreography for the concurrency itself: a mapper
  shared between worker processes.
- **Time, randomness and the environment are inputs.** A random sample is drawn from a seeded generator, so a failure
  reproduces — the property tests under `hgvs-weaver/tests/properties` print the seed that failed; the same for an
  environment variable or a hostname.
- **Each test owns its state.** Tests run in any order, and a parallel runner must stay possible, so a test assumes
  nothing about what another left behind and leaves nothing another could find: a unique path under `tmp_path`, a fresh
  mapper, no module global mutated and unrestored.

## Test doubles, in order of fidelity

A double's fidelity is how closely its behaviour matches production's. Prefer, in order: the real implementation; a
*fake*, a working lightweight implementation held to the real one's contract; and last a *stub* or *mock*, whose
behaviour is scripted in the test. Stop at the highest fidelity that keeps the test hermetic — in one process, or against a gated
local process ([below](#tests-run-with-no-cloud-and-no-network)). A hand-written stand-in scripted in the test body is a
stub or mock for this purpose, whatever it is called: its behaviour is the author's belief about the dependency.

- **A synthetic provider over a mocked one.** A `DataProvider` over a few hand-written transcripts on a synthetic
  chromosome is a fake: it answers every method from its own tables, so a test drives the real mapper end to end. A
  mock that scripts `get_transcript` for one call answers only the question the author thought of.
- **Don't mock types you don't own.** A mock of pysam, the biocommons `hgvs` package or the extension module encodes the
  author's reading of that library; an upgrade that changes its behaviour leaves the mock answering as before and ships
  the bug green. Use the real library against real inputs, or a fake its owner ships; failing both, wrap it in a class
  you own, mock the wrapper, and test the wrapper against the real library so the slow tests are confined to one suite.
- **Verify commands, never queries.** A test that verifies a getter was called proves nothing about what was done with
  its result and duplicates the assertion on the outcome. Assert the state or the return value. Verifying an interaction
  is right only when the interaction *is* the behaviour: a transcript model fetched once while the mapper caches it.
- **Verify only the arguments the behaviour depends on.** Pinning every argument of a verified call makes every test
  that names the call a casualty of any change to any argument, and hides which argument the test is about; match the
  rest.
- **Signs of a mock too many.** More than one or two collaborators mocked, a mock scripting more than a method or two,
  or a reader who has to step through the code under test to follow the test.

## Tests run with no cloud and no network

CI has no credentials and no reference data, so a test that needs a bucket, a live service or the GRCh38 files does not
run there — which, for a required check, reads as a pass. Tests are hermetic by construction: a synthetic provider for
the mapper, recorded answers for the biocommons comparison, and the ClinVar gate kept out of CI as the regression check
it is ([`../../CLAUDE.md`](../../CLAUDE.md) says when it runs).

Where a test does need real reference data, the gate is visible: it skips, naming the files it wanted, on a machine
without them. It never fails for their absence, and never passes vacuously by asserting nothing when they are absent.

## A fixture is an independent oracle

A fixture the code under test produced can only show that the code round-trips its own output; a defect shared between
the writing and reading halves is invisible to it. The same holds for the property a test asserts, not only the data it
asserts on ([above](#test-at-the-surface-a-caller-uses)). Build fixtures from something the code does not control: a
projection test checks weaver's answer against coordinates worked out from the synthetic transcript's own definition,
not against what an earlier weaver printed; the grammar is checked rule by rule against the biocommons `hgvs` grammar
table; the equivalence cases and the ClinVar sample are other tools' answers. Where another tool was run to settle a
case, [`choices.md`](../source/choices.md) records which, so the oracle stays named.

## Test data is public

This repository is public, and everything in it ships. Public reference data — an assembly's sequence, a published
annotation — is fine to excerpt; nothing that could plausibly have come from a patient or a case is. No variant tables,
phenotype text or sample identifiers, anonymised or not. Synthesise, and make the synthetic origin unmistakable: a
chromosome named `NC_000099.1`, genes named `FWD` and `REV`.

## Where tests live

The Rust crate's unit tests sit beside the code under `hgvs-weaver/src`, its integration tests under
`hgvs-weaver/tests`, and its property tests under `hgvs-weaver/tests/properties`. The Python tests live in `tests/` and
run against the extension built by `uv run --no-sync maturin develop`. A behaviour of the core is tested in Rust; the
Python suite tests the binding, the protocol and the commands, not the core again.
