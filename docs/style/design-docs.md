# Design docs

A design doc is the durable record of one area's design, read by people who were not in the conversation that produced
it. Everything in this guide follows from that — who reads the doc, and what they need to get from it. The reader and
the style below apply to every doc under `docs/`; the shape and the policy are for design docs.

## The reader

A design doc is normally reviewed by a second maintainer who was not in the conversation, and it is what every later
reader of the area starts from. So write for a maintainer who has read [`../PRODUCT.md`](../PRODUCT.md) and the README's
Overview, knows nothing about this area — not its adjacent docs, not its code — and has to make a review decision from
one read on GitHub: they should come away with the decisions and the reasons for them without having to ask anyone.

Writing for that human also serves a model reading the doc as context. Writing for the model does not serve the human:
compression strips the antecedents, signposts and consequences a person needs to follow an argument they are seeing for
the first time.

## Style

Explain with the clarity and style of Martin Kleppmann: plain prose built up from first principles, terms defined before
they are used, a concrete example where it aids understanding, and the trade-offs of each option stated honestly rather
than slanted towards the one chosen.

The forms that carry a design here are concrete. Where it has a surface, a mockup of that surface; where it has a flow
across services, a request diagram; a plain statement of what is stored where; and where the doc covers several
interfaces, a subsection per interface in one consistent shape — what it takes, what it produces, what it refuses — so a
reader can compare them. Such a subsection is design, not a restatement of code: what a message means in the system
belongs in the doc, its field list belongs to the type or the stub. A close-to-code sketch — a class or service with
method bodies elided to `...` — is a form on the same footing where the architecture is easiest to see that way: a guess
at writing time, superseded by the contract file or the code once that exists. These forms are additive to the prose,
not carved out of it: a reader takes in a mockup or a diagram in seconds, and it leaves the prose less to say. The
length that has to fit a sitting is the argument's; a doc that runs long tightens its prose or splits along a seam, and
keeps the picture.

The maintainers also prefer a particular presentation, and the reason is again that review decision. Motivation first,
then what was decided, then the mechanism at the level of concepts and interfaces — the specifics after it, in an
appendix where they would break the thread. Interleaved specifics cost the reviewer the thread: having read past a
paragraph of detail to reach the next step of the reasoning, they have to reconstruct where the reasoning was. Spell out
the implications rather than leaving them to be derived, and give enumerations as lists rather than run-on prose. A
point worth making gets its own section, with the reason it matters — never an aside in the middle of an argument about
something else. A mermaid diagram suits a flow or a set of states, though GitHub renders it in the file view and not in
the PR diff, so a reviewer has to open the file.

## Where the low-level detail goes

The doc names decisions and interfaces — "the `TranscriptMapper`", "the `DataProvider` protocol" — and where an
interface has an entry point in code, links it where the interface is introduced: the module, the crate, the stub file.
Anything beyond that restates code and stays out: a per-field paraphrase of a struct, a protocol or a schema, env-var
names, file paths beyond that entry point, function, class or test names, error strings, constants. A paraphrase is a
second copy of something that already has an authoritative one, and it goes stale the next time the code changes.

Cutting such a detail is not deleting it — it moves to the code's own documentation surface: the doc comment on a type
or field, a module or function docstring, a test. What is written there is written for the caller — what the field,
option or message implies for them, what they must do and what they can rely on — not the mechanism behind it. Not an
inline comment beside the implementation; those are what the Comments rule in [`general.md`](general.md) governs, and
its default is no comment. The doc links to that surface instead. A passage that stops making sense once the detail is
cut was pitched at the wrong altitude: rewrite it higher rather than restoring the detail.

The division holds when both surfaces carry the same decision: the doc states the decision and the reason for it, and
the docstring on the function or type states what it implies for the caller — what they must do, what they can rely on —
and at most a clause of why. Where the doc says a transcript-only base projects as the whole run because a one-base
insertion would be a frameshift on every other transcript, the docstring on the projection says that an edit touching
such a run is rewritten over the run. Both name the same promise; that is not a second copy, because the argument sits
in one place — it is the argument, restated in both, that drifts. The pointer runs the other way too: where context
helps the code's reader, the docstring or comment carries a bare one-line pointer to the doc — never a restated
explanation ([`general.md`](general.md), Comments).

The doc also does not narrate what a fixture, a test or a first implementation happened to constrain. Where an
implementation constraint genuinely shapes an interface, state it as a decision and give its reason; otherwise leave it
out. A paragraph that explains an interface by what some implementation does reads to a reviewer as a design smell — the
implementation constraining the interface rather than the other way round.

A doc whose contract is the only code so far is held to the same line: that file — the Python protocol, the stub — is
the spec, and the doc explains it. When nothing is written yet, the doc names the interfaces and says what each one
promises; the field-level detail is written once, into the contract file, when that file exists.

## A default shape

Not a template — the shape that has served, to deviate from where the design reads better another way. A header line
carries `**Related:**` links, each with a few words on what the linked doc covers so a reader can tell whether to open
it. Then:

- **Overview** — the take-aways for a reader who reads nothing else: what this is, what was decided, why.
- **Background** — what the reader needs in order to read the rest: the vocabulary of the area, the problem, and the
  constraints that bound the design, in plain English — the behaviour a third-party API or package was found to have,
  what is non-negotiable from the user's perspective; not how the design got here.
- **Non-goals** — what the design will not do: the deliberate scope exclusions, each with the reason for it. A *not yet
  built* is not a non-goal — the doc describes the design; what of it exists is the code's to show.
- **Design** — the decision, then the mechanism at the level of concepts and interfaces, then the consequences.
- **Alternatives considered** — each option weighed, and the reason it was rejected.
- **Open questions** — unresolved points that need a decision or an input; omit the section when there are none.
- **Appendix** — supporting material that would break the thread: a worked example, a measurement, an evidence table.
  Never a restatement of code.

No section states what is built, what is planned or what has shipped: that changes with every merge, and the code and
git already hold it, so a section stating it rots at once. What such a section would carry belongs elsewhere — a
deliberate deferral and its reason is a design decision (Design, or Alternatives considered), an unresolved gap is an
Open question, and an accepted transitional cost, a deprecated field kept for one release or a two-step retirement, is
stated beside the decision that accepts it.

## Policy

- **One living doc per area**, under `docs/source/`, where Sphinx builds it. Two exist:
  [`choices.md`](../source/choices.md), the record of every choice weaver makes where another tool could reasonably
  answer otherwise, in its own fixed shape (example, agreement, specification, tests); and
  [`equivalence_logic.md`](../source/equivalence_logic.md). A new area gets a doc beside them in the default shape
  above; a new choice is an entry in `choices.md`. No ADRs: rationale lives in the doc's `Alternatives considered`,
  chronology in git.
- **Rewrite in place.** A design change edits the doc — never appends a supersession layer, never spawns a "v2". A
  superseded doc is deleted and its live content folded into its successor, in the same PR.
- **State each decision once** within the doc; restating one three sections later creates a second copy, and the next
  edit will update only one.
- **A doc is followable on its own.** Where the argument turns on another doc's decision, state that decision in a
  clause and link.
- **Short enough to read in a sitting.** A doc well past that usually covers two areas, and splits along its own
  subsection boundaries: "one living doc per area" fixes the number of docs about an area, not the size of an area.
  Where the *how* of a sub-design stands on its own, this doc says why it matters and links to a doc of its own for the
  mechanism.
- **A doc that does not yet meet this guide is brought to it** when next substantially edited, or in a rewrite of its
  own; the PR-time review looks only at the docs a PR changes.
