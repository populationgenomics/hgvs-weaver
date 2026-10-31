"""The commands an install of hgvs-weaver declares name functions that exist."""

from __future__ import annotations

import importlib.metadata

import pytest

COMMANDS = sorted(
    (entry.name, entry.value)
    for entry in importlib.metadata.distribution("hgvs-weaver").entry_points
    if entry.group == "console_scripts"
)


def test_the_install_declares_commands() -> None:
    assert COMMANDS, "the installed hgvs-weaver declares no console scripts"


@pytest.mark.parametrize(("command", "target"), COMMANDS)
def test_an_installed_command_names_a_function_that_exists(command: str, target: str) -> None:
    module_name, _, function = target.partition(":")
    module = importlib.import_module(module_name)
    assert callable(getattr(module, function, None)), f"{command}: {target} is not a function"
