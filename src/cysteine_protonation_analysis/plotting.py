"""Plotters for pKa traces and distributions."""
from __future__ import annotations
from pathlib import Path
from .io import PKaSeries


def plot_series(series: PKaSeries, output_prefix: str | Path, bins: int = 25) -> tuple[Path, Path]:
    if bins < 1:
        raise ValueError("bins must be positive")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    trace_path = Path(str(prefix) + "_trace.png")
    hist_path = Path(str(prefix) + "_histogram.png")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(series.frames, series.values, linewidth=0.9)
    ax.set_xlabel("Frame label (not physical time unless mapped)")
    ax.set_ylabel("Predicted cysteine pKa")
    ax.set_title(series.source.name)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(trace_path, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.hist(series.values, bins=bins, edgecolor="black", linewidth=0.4)
    ax.set_xlabel("Predicted cysteine pKa")
    ax.set_ylabel("Frame count")
    ax.set_title(f"Distribution: {series.source.name} (n={len(series.values)})")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(hist_path, dpi=180)
    plt.close(fig)
    return trace_path, hist_path
