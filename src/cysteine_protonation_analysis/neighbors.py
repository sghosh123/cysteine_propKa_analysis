"""Periodic-boundary-aware heavy-atom neighbor analysis around a cysteine atom."""
from __future__ import annotations
from pathlib import Path
import csv
import math

WATER_NAMES = {"HOH", "WAT", "TIP3", "TIP3P", "TIP4", "TIP4P", "SOL", "SPC", "SPCE", "TIP5P"}
COMMON_METALS = {"CU", "CU1", "CU2", "ZN", "FE", "FE2", "FE3", "MG", "MN", "CO", "NI", "CA", "NA", "K", "CD", "HG"}


def _safe_attr(obj, attribute: str, fallback=""):
    try:
        value = getattr(obj, attribute)
    except Exception:
        return fallback
    return fallback if value is None else value


def _chain_id(atom) -> str:
    return str(_safe_attr(atom, "chainID", "") or "")


def _is_hydrogen(atom) -> bool:
    element = str(_safe_attr(atom, "element", "") or "").strip().upper()
    if element:
        return element in {"H", "D", "T"}
    name = "".join(c for c in str(atom.name).upper() if not c.isdigit())
    return name.startswith(("H", "D"))


def _category(atom) -> str:
    resname = str(getattr(atom, "resname", "")).upper()
    element = str(_safe_attr(atom, "element", "") or "").strip().upper()
    name = str(getattr(atom, "name", "")).upper()
    if resname in WATER_NAMES:
        return "water"
    if resname in COMMON_METALS or element in COMMON_METALS or name in COMMON_METALS:
        return "metal_or_ion"
    return "other"


def _valid_box(dimensions) -> bool:
    if dimensions is None or len(dimensions) < 6:
        return False
    try:
        vals = [float(v) for v in dimensions[:6]]
    except (TypeError, ValueError):
        return False
    return all(math.isfinite(v) for v in vals) and all(v > 0 for v in vals[:3]) and all(0 < v < 180 for v in vals[3:6])


def cysteine_atom_contacts(universe, resid: int, atom_name: str = "SG", segid: str | None = None,
                           chain: str | None = None, cutoff: float = 4.0,
                           include_hydrogens: bool = False) -> list[dict[str, object]]:
    """Return long-form per-frame contacts for one uniquely identified atom.

    ``cutoff`` is in Å because MDAnalysis uses Å for standard coordinate readers.
    The target is identified by resid, atom name, and optional segid/chain ID. If
    residue numbers repeat, pass whichever identifier is available in the topology.
    """
    if cutoff <= 0 or not math.isfinite(cutoff):
        raise ValueError("cutoff must be a positive finite distance in Å")
    import numpy as np
    from MDAnalysis.lib.distances import capped_distance

    atoms = universe.atoms
    mask = (atoms.resids == int(resid)) & (atoms.names == str(atom_name))
    if segid is not None:
        mask &= atoms.segids == str(segid)
    if chain is not None:
        try:
            mask &= atoms.chainIDs == str(chain)
        except Exception as exc:
            raise ValueError("A chain was requested, but this topology does not expose chainIDs; use --segid if appropriate") from exc
    target_indices = np.flatnonzero(mask)
    if len(target_indices) != 1:
        qualifier = f", segid={segid!r}" if segid is not None else ""
        qualifier += f", chain={chain!r}" if chain is not None else ""
        raise ValueError(f"Expected exactly one atom for resid={resid}, name={atom_name!r}{qualifier}; found {len(target_indices)}")
    target = atoms[int(target_indices[0])]
    mobile_indices = np.flatnonzero(atoms.resindices != target.resindex)
    if not include_hydrogens:
        mobile_indices = np.asarray([i for i in mobile_indices if not _is_hydrogen(atoms[int(i)])], dtype=int)
    mobile = atoms[mobile_indices]
    rows: list[dict[str, object]] = []
    if len(mobile) == 0:
        return rows
    for ts in universe.trajectory:
        box = ts.dimensions if _valid_box(ts.dimensions) else None
        pairs, distances = capped_distance(
            target.position.reshape(1, 3), mobile.positions, max_cutoff=float(cutoff),
            box=box, return_distances=True,
        )
        time_value = float(ts.time) if np.isfinite(getattr(ts, "time", np.nan)) else ""
        for pair, distance in zip(pairs, distances):
            partner = mobile[int(pair[1])]
            rows.append({
                "frame_index_zero_based": int(ts.frame),
                "time_value_from_reader": time_value,
                "target_segid": str(target.segid),
                "target_chain": _chain_id(target),
                "target_resid": int(target.resid),
                "target_resname": str(target.resname),
                "target_atom": str(target.name),
                "partner_segid": str(partner.segid),
                "partner_chain": _chain_id(partner),
                "partner_resid": int(partner.resid),
                "partner_resindex": int(partner.resindex),
                "partner_resname": str(partner.resname),
                "partner_atom": str(partner.name),
                "partner_element": str(_safe_attr(partner, "element", "") or ""),
                "partner_category": _category(partner),
                "distance_A": float(distance),
            })
    return rows


def write_contacts_csv(rows: list[dict[str, object]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["frame_index_zero_based", "time_value_from_reader", "target_segid", "target_chain", "target_resid", "target_resname", "target_atom",
              "partner_segid", "partner_chain", "partner_resid", "partner_resindex", "partner_resname", "partner_atom", "partner_element",
              "partner_category", "distance_A"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
