"""The commands pyproject.toml installs name functions that exist."""

from __future__ import annotations

import importlib
import pathlib
import tomllib

import pytest

PYPROJECT = pathlib.Path(__file__).parent.parent / "pyproject.toml"
SCRIPTS = tomllib.loads(PYPROJECT.read_text())["project"]["scripts"]


@pytest.mark.parametrize(("command", "target"), sorted(SCRIPTS.items()))
def test_an_installed_command_names_a_function_that_exists(command: str, target: str) -> None:
    module_name, _, function = target.partition(":")
    module = importlib.import_module(module_name)
    assert callable(getattr(module, function, None)), f"{command}: {target} is not a function"
