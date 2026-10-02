# Python style

This is the style guide for Python code in hgvs-weaver: the `weaver` package, its CLIs and its tests. It's the
high-level / human-judgement layer — the mechanical formatting and naming rules are enforced by [ruff] and aren't
restated here. If something contradicts what ruff is configured to enforce, ruff wins, and this doc should be updated.

The content is largely lifted from the [Google Python Style Guide][pyguide] (CC-BY-3.0), edited for the conventions of
this repo. Sections we don't enforce or that ruff already covers have been omitted.

## Imports

**Import modules, not symbols.** Use `import x` to bring in a package or module. Use `from x import y` only when `y` is
itself a module (i.e. `x` is the package prefix). Do not use `from x import y` to bring in a class, function, constant,
or other symbol.

Why: at the call site, a qualified reference (`module.Thing`) makes it clear where `Thing` came from. A bare `Thing` is
ambiguous to anyone reading the file — and to AI assistants navigating the code — without scrolling to the import block.

```python
# Good
import pathlib
import dataclasses

path = pathlib.Path('/etc/hostname')

@dataclasses.dataclass
class Record:
    name: str

# Bad — symbols pulled out of their module
from pathlib import Path
from dataclasses import dataclass

path = Path('/etc/hostname')

@dataclass
class Record:
    name: str
```

The cost is a little extra typing at call sites. It's worth it.

### Carved-out exceptions

The following symbol imports are allowed:

- Names from the `typing` module (`Any`, `Protocol`, `TypeVar`, `NamedTuple`, etc.). These are language vocabulary;
  fully-qualifying them adds noise without clarity.

- Names from `collections.abc` (`Iterator`, `Iterable`, `Mapping`, `Sequence`, etc.). Same rationale.

- `from __future__ import ...` — required syntax for future statements, not a normal import.

Anything not on this list — including `pathlib.Path`, `dataclasses.dataclass`, `contextlib.contextmanager`,
`functools.partial`, third-party classes, and this package's own classes — should be accessed via its module. A private
module is imported the same way (`from weaver.cli import provider`, then `provider.RefSeqDataProvider`); where a
module's name would collide with a local, alias it (`from weaver.cli import validate as validate_cli`). Imports are
absolute (ruff `TID252`); never `from ._sibling import`.

### No `TYPE_CHECKING` blocks

Don't use `from typing import TYPE_CHECKING` to gate type-only imports. The pattern is an ugly two-tier import structure
and almost always papers over a problem that's better fixed elsewhere:

- If the import is type-only because of a *circular import*, the abstractions are usually wrong. Restructure the
  modules.
- If the import is type-only because the library is *expensive to import* at runtime (`hail`, `torch`, `sqlalchemy`),
  the legitimate exceptions are rare in this codebase — handle them case-by-case with an inline comment justifying the
  block.

For everything else (`Iterator`, `Iterable`, `Sequence` from `collections.abc`), import at module level. The runtime
cost is nothing.

```python
# Good
from collections.abc import Iterator

def lines(source: str) -> Iterator[str]:
    ...

# Bad — Iterator costs nothing to import; the block is pure noise
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from collections.abc import Iterator

def lines(source: str) -> Iterator[str]:
    ...
```

## Exception handling

Make use of built-in exception classes when it makes sense. For example, raise a `ValueError` if you were passed a
negative number but were expecting a positive one. Do not use `assert` statements for validating argument values of a
public API. `assert` is used to ensure internal correctness, not to enforce correct usage nor to indicate that some
unexpected event occurred.

Narrow `except` clauses to the specific exceptions you can actually handle. Never use bare `except:` — it catches
`KeyboardInterrupt` and `SystemExit` too. Never silently swallow an exception; at minimum log it, ideally re-raise with
added context using `raise ... from`.

```python
# Good
try:
    return parse_config(path)
except json.JSONDecodeError as e:
    raise ConfigError(f'invalid JSON in {path}') from e

# Bad
try:
    return parse_config(path)
except Exception:
    return {}  # caller has no way to know the file was invalid
```

## Mutable global state

Avoid module-level mutable state. Module globals that get mutated at runtime are effectively a hidden parameter to every
function in the module — they make code hard to test, hard to reason about under concurrency, and surprising when the
module is imported twice.

Constants are fine. Caches, registries, and similar should usually be attached to an explicit object or scoped to a
function.

```python
# Good — caller owns the state
def build_index(records: Iterable[Record]) -> dict[str, Record]:
    return {r.id: r for r in records}

# Bad — every caller now shares this dictionary
_INDEX: dict[str, Record] = {}

def add_to_index(r: Record) -> None:
    _INDEX[r.id] = r
```

## Resource management

Always use the `with` statement for resources that need explicit cleanup — files, sockets, network connections, locks,
database transactions, subprocesses. Don't rely on garbage collection to close them.

```python
# Good
with open(path) as f:
    data = f.read()

# Bad — file stays open until the GC gets to it
f = open(path)
data = f.read()
```

For objects without a context-manager protocol, wrap them with `contextlib.closing` rather than calling `.close()`
manually in a `finally` block.

## Power features

Avoid metaclasses, monkey-patching, dynamic class creation, `eval`, `exec`, custom `__getattr__` / `__setattr__` magic,
and reflective hacks. These features exist for legitimate library-level use cases, but in application code they almost
always make the code harder to read, harder to type-check, and harder for tools (including AI assistants) to reason
about.

If you find yourself reaching for a power feature, the problem is usually that an abstraction is missing or wrong; fix
that instead.

## Docstrings

Use Google-style docstrings on public APIs that warrant explanation. Public APIs *should* have docstrings; this is
policy, not lint — ruff's `D` rules are not selected here, so neither a missing docstring nor its format fails CI.

Write them in the Google convention anyway (one-line summary, blank line, body, closing quotes on their own line).
Beyond format, the *content* matters:

- **Module docstring**: what the module is for; what to use it for.
- **Function/method docstring**: what the caller needs to know to call it correctly. Not a paraphrase of the
  implementation.
- **Class docstring**: the invariant the class maintains, plus any important behaviour the class name doesn't already
  imply.

The standard sections are `Args:`, `Returns:`, `Raises:`, `Yields:`.

```python
# Good
def merge_base(head: str, base: str) -> str:
    """Resolve the merge base of two commits.

    Args:
        head: SHA of the PR head commit.
        base: SHA of the PR base commit (usually `main`).

    Returns:
        The merge-base SHA as a hex string.

    Raises:
        subprocess.CalledProcessError: If `git merge-base` returns
            a non-zero exit code (commits unknown, no common ancestor).
    """
    ...

# Bad — restates the code without adding signal
def merge_base(head: str, base: str) -> str:
    """Call git merge-base with head and base and return the output."""
    ...
```

## Type annotations

Annotations must be present on function and method signatures (ruff's `ANN` rules) and should be *specific*, which is
policy, not lint; the type checker configured in `pyproject.toml` checks what is annotated.

- Prefer concrete types over `typing.Any`. Ruff's `ANN401` forbids it outside the CLI modules and the extension's stubs,
  where it is allowed; where `Any` is the honest type there, say why in a one-line comment, and consider whether a
  `Protocol`, `TypeVar`, or `Union` would carry more signal.
- Prefer `list[Dataset]` over bare `list`. The element type carries meaning.
- Prefer `Iterable[T]` / `Sequence[T]` / `Mapping[K, V]` for function parameters (the most general type the function
  actually uses); use concrete `list` / `dict` for return types.
- Don't annotate the obvious. `def __init__(self, name: str) -> None` is fine; `def f(self) -> Self` is fine;
  `def add(a: int, b: int) -> int: return a + b` is fine. Don't add `# type: Iterable[int]` comments inside a function
  unless the type checker is confused.
- The extension module's surface is declared in `weaver/_weaver.pyi`. A change to a pyo3 signature changes the stub in
  the same commit; the stub is what every annotation in the Python package resolves against.

## Function decomposition

A function should do one thing. Pyguide's loose target is "under 40 lines"; if you're well above that, the function is
probably doing more than one thing. Decompose along natural boundaries — input parsing, the core operation, formatting
the output — rather than chopping arbitrarily.

A function that's structurally simple but long (e.g. a config generator with 100 sequential `setattr` calls) is fine. A
function that's 50 lines of branching, where each branch is a different mode of operation, is not — those branches want
to be separate functions.

## Tooling notes

### `from __future__ import annotations`

Use it. Every Python module in this repo should start with `from __future__ import annotations` (after the module
docstring, if present). It makes all annotations lazily-evaluated strings, which:

- Removes runtime ordering dependencies for forward references.
- Lets us use `list[int]` / `X | Y` annotations without worrying about Python version.
- Means the cost of an unused type-only import is the import cost only, not the resolution cost — supporting the "no
  `TYPE_CHECKING`" stance above.

### pytest

Test files use `pytest`, not `unittest`. Test functions are top-level or grouped under classes named `Test*`. Assertions
use bare `assert` (ruff's `S101` is ignored under `tests/`). The extension is built into the environment with
`uv run --no-sync maturin develop` before the tests import it.

```python
# Good — pytest style
def test_merge_base_returns_a_sha():
    assert _looks_like_a_sha(merge_base('HEAD~1', 'HEAD'))

# Bad — unittest style in this repo
class MergeBaseTest(unittest.TestCase):
    def test_returns_a_sha(self):
        self.assertTrue(...)
```

Use `pytest` fixtures (`tmp_path`, `monkeypatch`, `capsys`, custom-defined) for setup, not class-level setup methods —
for the setup a test does not assert on; a value the assertion depends on is set in the test body
([`writing-tests.md`](writing-tests.md#a-test-reads-without-its-reader-leaving-it)). Parametrise tests with
`@pytest.mark.parametrize` rather than writing N near-identical test functions.

What a test asserts — invariants rather than the current shape, results rather than interaction sequences — and what it
may depend on or contain is language-independent and lives in [`writing-tests.md`](writing-tests.md).

______________________________________________________________________

*Adapted from the [Google Python Style Guide][pyguide], licensed under [CC-BY-3.0]. Sections covered by [ruff]'s lint
rules have been omitted; the linter config in `pyproject.toml` is the source of truth for those.*

[cc-by-3.0]: https://creativecommons.org/licenses/by/3.0/
[pyguide]: https://google.github.io/styleguide/pyguide.html
[ruff]: https://docs.astral.sh/ruff/
