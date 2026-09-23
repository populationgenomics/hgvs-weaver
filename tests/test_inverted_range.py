"""A range written backwards is refused when parsed, as an HGVSError a caller can catch."""

import pytest

import weaver


@pytest.mark.parametrize(
    "hgvs",
    ["NM_206933.4:c.8559_2A>G", "NM_025137.4:c.6331_6232insG", "NM_025137.4:c.100_50del"],
)
def test_a_range_written_backwards_is_a_parse_error(hgvs: str) -> None:
    with pytest.raises(weaver.ParseError, match="runs backwards") as caught:
        weaver.parse(hgvs)
    assert isinstance(caught.value, weaver.HGVSError)
    assert isinstance(caught.value, Exception)


def test_the_ordered_control_parses() -> None:
    assert str(weaver.parse("NM_025137.4:c.50_100del")) == "NM_025137.4:c.50_100del"
