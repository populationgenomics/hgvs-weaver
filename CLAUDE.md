# hgvs-weaver development notes

## Product

See [`docs/PRODUCT.md`](docs/PRODUCT.md) for the product north star — what the engine is, the load-bearing principles
(positions are types; consequences are read from the sequence; every choice is written down with who agrees; an answer
or an error, never a guess), and what is out of scope. Read it before proposing designs or plans. The choices weaver
makes, each with its example, specification, tests and the tools that agree, are in
[`docs/source/choices.md`](docs/source/choices.md).

## Working norms

Operating directives for Claude (and any agent) in this repo; they counteract default model dispositions.

- **Resist the minimal-diff reflex.** Don't reach for the smallest change that hides the symptom (special-casing,
  papering over root causes). Aim for the correct fix at the right complexity level — not the smallest, not gold-plated.
- **Fail loudly and early.** Raise on a missing expected input or precondition; never fall back to a default/placeholder
  to limp along. A placeholder is an explicit caller input, never a code default.
- **Never instruct around a defect — fix the defect.** Don't write prose telling readers to work around broken code —
  "pass it as a string, the converter loses precision". Prose is untested, and callers who didn't read it stay broken.
- **Push back; don't just comply.** When a design, name, or approach seems worse — including a shortcut you're asked to
  take — say so with reasoning, unprompted. The author owns the final call.
- **Offer better alternatives with trade-offs.** When a materially better approach than the proposed one exists, present
  it and the trade-offs — don't just execute the ask.
- **Investigate before producing.** Read the code and verify constraints first. Don't treat a training-pattern
  convention as load-bearing unchecked; don't speculate about what you can read.
- **Explain non-obvious changes first.** For a change whose rationale isn't self-evident, give the why before showing or
  applying the diff.
- **Ask when unsure** rather than assume intent.
- **No intensifiers or emphasis filler.** Drop words and phrases that add emphasis but no information — "that's the
  key", "crucially", "importantly", "the key insight", "it's worth noting". State the point plainly. Applies to all
  prose: chat replies, PR/review comments, commit messages, and docs.

## Code style

@docs/style/general.md

Python follows [`docs/style/python.md`](docs/style/python.md); it loads when a Python file is touched. Rust's mechanical
layer is `cargo fmt` and `cargo clippy`; the judgement layer is `general.md`, and the crate's own conventions — tagged
position types, `thiserror` errors a caller can match, no allocation on a hot path that a stack buffer serves — are read
from the code around the change.

Tests — what one asserts, what it may depend on, what its data may contain — follow
[`docs/style/writing-tests.md`](docs/style/writing-tests.md); it loads when a test file is touched.

## Core and binding

The engine is the `hgvs-weaver` crate under `hgvs-weaver/`; the root crate is the pyo3 binding, and `weaver/_weaver.pyi`
declares what it exposes. A behaviour lives in the core and is tested in Rust; the Python suite tests the binding, the
protocols and the commands. A change to a pyo3 signature changes the stub in the same commit. Build the extension into
the environment with `uv run --no-sync maturin develop` before running pytest; do not rebuild while a validation run is
in progress, since its workers import the shared library.

A change to the core — parsing, normalisation, projection, consequences, equivalence — is gated by the 100,000-variant
ClinVar run: `uv run --group gate weaver-gate` builds the reference store once and then costs the validation alone, and
prints the row-by-row comparison with a baseline output. Report the rows that moved, and check each against the truth
columns; the percentages alone hide a single moved row. A new choice, or a changed one, is an entry in `choices.md` in
its fixed shape.

## Docs

Two audiences, two registers:

- **Instruction files** are prompts and rules — `CLAUDE.md`, `.claude/rules/`, `.claude/skills/`: model-only, only what
  changes behavior, no maintainer notes, no harness mechanics (which rules load when, where files live). A token there
  is paid on every run that loads it; human-facing explanation belongs in `docs/` or code.
- **Everything under `docs/`** is written for a human first — a maintainer who has read
  [`docs/PRODUCT.md`](docs/PRODUCT.md) and the README's Overview but not this area, and has to get the take-aways from
  one read on GitHub. Explain with the clarity and style of Martin Kleppmann — motivation before mechanism, specifics
  out of the argument's way. Detail that restates code — field lists, paths, env vars, test names — stays in the code
  and is linked, never transcribed. A model reads what a human reads. `docs/source/` is the Sphinx site and the design
  record: one living doc per area, rewritten in place; no ADRs — rationale lives in the doc, chronology in git. The
  guide is [`docs/style/design-docs.md`](docs/style/design-docs.md); to write or rewrite one, load the
  `writing-design-docs` skill. `python3 tools/check_links.py` checks every relative link in the tracked Markdown.

## Committing

- **Stage explicit paths**, not `git add -A` / `.`; explicit staging avoids sweeping in an untracked file, and the
  working tree may hold the author's uncommitted edits in files a change also touches.
- **Pre-commit runs lint/format/hygiene** (`.pre-commit-config.yaml`); CI runs the same hooks plus `cargo test` and
  pytest. Ensure hooks are installed (`uv run --no-sync pre-commit install`) — if not, install or ask the author; never
  bypass with `--no-verify`.
- **Lock discipline.** CI syncs from `uv.lock` with `--locked`; a change to `pyproject.toml`'s dependencies comes with
  the relocked file in the same commit. `Cargo.lock` at the root and in `hgvs-weaver/` are tracked for the same reason.
- **Correct a pushed branch with a new commit on top**, not amend + force-push. PRs squash-merge, so `main` history
  stays linear regardless and intermediate fixups vanish on merge. Reserve force-push for rebasing a branch onto `main`.

## Worktrees

Worktrees go in `.claude/worktrees/` (gitignored), never `../` siblings.

- **New branch** → the Claude Code worktree command.
- **Existing branch** → `git worktree add .claude/worktrees/<name> <branch>` (the command only cuts fresh branches).

## CI and review

- **Adversarially review before opening a PR.** For any change with non-trivial code or logic, run adversarial review
  passes in subagents with fresh context — the reviewer sees only the diff, not the authoring conversation — and fix the
  findings autonomously; repeat until a pass surfaces only diminishing findings, then open the PR. Exempt: trivial
  changes, doc-only changes, resource/asset changes.
- **A PR description is written for the human reviewer**: what the change is and why, the take-aways, and where to look
  — the altitude of a design doc's Overview, shorter. The diff carries the detail; don't narrate it. Same style:
  [`docs/style/design-docs.md` § Style](docs/style/design-docs.md#style).
- **Pin third-party GitHub Actions to a commit SHA** with the version in a trailing comment (`@3d3c…  # v7.0.1`), never
  a moving tag: a tag can be moved to malicious code, a SHA cannot. Resolve the SHA from the action's release when
  adding or bumping one, and prefer a tool installed from a checksum-verified release over a third-party action where
  the action would only wrap a download.
