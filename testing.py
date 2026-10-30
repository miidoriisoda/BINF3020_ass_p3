import pytest
from helper_functions import jukes_cantor, kimura_two_parameter

SEQ1 = "AGGT-GA"
SEQ2 = "AGC-AGA"


def test_jukes_cantor():
    assert jukes_cantor(SEQ1, SEQ2) == pytest.approx(1.163, rel=1e-3)


def test_kimura_two_parameter():
    assert kimura_two_parameter(SEQ1, SEQ2) == pytest.approx(1.196, rel=1e-3)