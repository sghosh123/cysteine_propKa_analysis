# Release checklist

- [ ] Choose a repository license before making the repository public.
- [ ] Verify redistribution permission for the included computed result files.
- [ ] Inspect the actual PDB snapshots for chain IDs, residue numbering, missing atoms, and metal/ligand treatment.
- [ ] Record the exact PROPKA version and command line for any new analysis.
- [ ] Validate one real pKa frame end to end using the local VMD and PROPKA installations.
- [ ] Verify frame spacing and frame-to-time mapping from the original trajectory/configuration.
- [ ] Before comparing data sets, confirm the residues refer to homologous positions and their chain IDs have consistent meanings.
- [ ] Do not claim copper-aware estimates unless the calculation explicitly and correctly models the copper coordination environment.
- [ ] Confirm all tests pass from a clean environment.
