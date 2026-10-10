"""Descriptive statistics for framewise predicted cysteine pKa values."""
from __future__ import annotations
import numpy as np
from .io import PKaSeries


def henderson_hasselbalch_deprotonated_fraction(pka: np.ndarray | float, ph: float) -> np.ndarray:
    """Return the idealized thiolate fraction for a monoprotic acid at ``ph``.

    This is an algebraic Henderson–Hasselbalch mapping of a predicted pKa. It is
    not a calibrated estimate of an MD ensemble's actual protonation population.
    """
    if not np.isfinite(ph):
        raise ValueError("pH must be finite")
    values = np.asarray(pka, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("pKa values must be finite")
    x = (values - float(ph)) * np.log(10.0)
    return np.exp(-np.logaddexp(0.0, x))


def summarize_series(series: PKaSeries, ph: float | None = None) -> dict[str, object]:
    """Calculate descriptive, non-inferential metrics for one pKa series."""
    x = series.values
    if x.size == 0:
        raise ValueError("Cannot summarize an empty series")
    result: dict[str, object] = {
        "dataset": series.source.name,
        "source": str(series.source),
        "n_frames": int(x.size),
        "first_frame_label": int(series.frames[0]),
        "last_frame_label": int(series.frames[-1]),
        "mean_pka": float(np.mean(x)),
        "median_pka": float(np.median(x)),
        "std_pka_sample": float(np.std(x, ddof=1)) if x.size > 1 else 0.0,
        "min_pka": float(np.min(x)),
        "q05_pka": float(np.quantile(x, 0.05)),
        "q25_pka": float(np.quantile(x, 0.25)),
        "q75_pka": float(np.quantile(x, 0.75)),
        "q95_pka": float(np.quantile(x, 0.95)),
        "max_pka": float(np.max(x)),
        "lag1_autocorrelation": _lag1_correlation(x),
    }
    if ph is not None:
        fractions = henderson_hasselbalch_deprotonated_fraction(x, ph)
        result.update({
            "reference_pH": float(ph),
            "fraction_frames_with_pka_below_pH": float(np.mean(x < ph)),
            "mean_HH_thiolate_fraction_model_only": float(np.mean(fractions)),
            "median_HH_thiolate_fraction_model_only": float(np.median(fractions)),
        })
    return result


def _lag1_correlation(values: np.ndarray) -> float | str:
    if values.size < 3 or np.std(values) == 0:
        return ""
    a, b = values[:-1], values[1:]
    if np.std(a) == 0 or np.std(b) == 0:
        return ""
    return float(np.corrcoef(a, b)[0, 1])
