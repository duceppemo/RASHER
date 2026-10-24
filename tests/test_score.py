"""The validation scoring script (validation/score.py, not part of the package)."""

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("score", Path(__file__).parents[1] / "validation" / "score.py")
score = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score)


def _truth():
    # A and B differ by one SNP outside the repeats, C by one more inside the inverted repeat (two sites)
    def matrix(ab, ac, bc):
        return {"A": {"A": 0, "B": ab, "C": ac}, "B": {"A": ab, "B": 0, "C": bc}, "C": {"A": ac, "B": bc, "C": 0}}
    return {"scenario": "plastid", "distances": {"events": matrix(1, 2, 1), "sites": matrix(1, 3, 2),
                                                  "no_repeat": matrix(1, 1, 0)}}


def test_no_repeat_best_fit_is_named_with_the_exact_pairs():
    truth = _truth()
    blind = {"A": {"A": 0, "B": 1, "C": 1}, "B": {"A": 1, "B": 0, "C": 0}, "C": {"A": 1, "B": 0, "C": 0}}
    events = {"A": {"A": 0, "B": 1, "C": 2}, "B": {"A": 1, "B": 0, "C": 1}, "C": {"A": 2, "B": 1, "C": 0}}
    results = [{"run": "r", "assembly": {}, "snps": {"blind": score.score_snps(blind, truth),
                                                      "counted": score.score_snps(events, truth)}}]
    s = results[0]["snps"]["blind"]
    assert (s["best_convention"], s["exact_pairs"], s["wrong_pairs"]) == ("no_repeat", 3, [])
    md = score.markdown(results, truth)
    assert "| r | blind | 3/3 | no_repeat | 3/3 (no_repeat; events 1/3, sites 1/3) | 0 |  |" in md
    assert "| r | counted | 3/3 | events | 3/3 | 0 |  |" in md
    assert "leaves out SNPs inside repeats" in md
    results[0]["snps"].pop("blind")
    assert "leaves out SNPs inside repeats" not in score.markdown(results, truth)
