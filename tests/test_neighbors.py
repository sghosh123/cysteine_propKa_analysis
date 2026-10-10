import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")
from MDAnalysis.coordinates.memory import MemoryReader
from cysteine_protonation_analysis.neighbors import cysteine_atom_contacts


def make_universe(duplicate_sg=False):
    u = mda.Universe.empty(
        n_atoms=4, n_residues=2, atom_resindex=[0, 0, 1, 1],
        n_segments=2, residue_segindex=[0, 1], trajectory=True,
    )
    u.add_TopologyAttr("names", ["SG", "SG" if duplicate_sg else "CB", "O", "H"])
    u.add_TopologyAttr("resnames", ["CYS", "HOH"])
    u.add_TopologyAttr("resids", [18, 101])
    u.add_TopologyAttr("segids", ["A", "W"])
    u.add_TopologyAttr("elements", ["S", "C", "O", "H"])
    coords = np.array([[[0.2, 0, 0], [0.4, 0, 0], [9.8, 0, 0], [9.7, 0, 0]]], dtype=np.float32)
    dims = np.array([[10, 10, 10, 90, 90, 90]], dtype=np.float32)
    u.load_new(coords, format=MemoryReader, dimensions=dims)
    return u


def test_neighbor_contacts_use_minimum_image_and_exclude_hydrogens():
    rows = cysteine_atom_contacts(make_universe(), resid=18, atom_name="SG", segid="A", cutoff=1.0)
    assert len(rows) == 1
    assert rows[0]["partner_atom"] == "O"
    assert rows[0]["partner_category"] == "water"
    assert rows[0]["distance_A"] == pytest.approx(0.4, abs=1e-3)


def test_neighbor_contacts_can_include_hydrogens():
    rows = cysteine_atom_contacts(make_universe(), resid=18, atom_name="SG", segid="A", cutoff=1.0,
                                 include_hydrogens=True)
    assert len(rows) == 2


def test_ambiguous_or_missing_target_fails():
    u = make_universe(duplicate_sg=True)
    with pytest.raises(ValueError, match="exactly one"):
        cysteine_atom_contacts(u, resid=18, atom_name="SG", segid="A", cutoff=1.0)
    with pytest.raises(ValueError, match="exactly one"):
        cysteine_atom_contacts(u, resid=999, atom_name="SG", cutoff=1.0)
