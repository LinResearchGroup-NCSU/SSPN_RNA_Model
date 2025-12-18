# SSPN_RNA_Model
In this repository, we present the RNA modeling capabilities of our Single-Site-Per-Nucleotide (SSPN) model, a modification to OpenABC. The modifications are described and used in `tutorial/single_rna.ipynb`. For an example of simulation of multiple RNAs, see `tutorial/multiple_rnas.ipynb`.

Helpful references:
---
Original paper for OpenABC

OpenABC: Liu, S. et al., **OpenABC enables flexible, simplified, and efficient GPU accelerated simulations of biomolecular condensates.** *PLoS Computational Biology* (2023). https://doi.org/10.1371/journal.pcbi.1011442.

---
---
Our RNA model uses only interactions from the SMOG forcefield, which was originally designed for proteins

SMOG: Noel, J. et al., **SMOG 2: a versatile software package for generating structure-based models.** *PLoS Computational Biology* (2016). https://doi.org/10.1371/journal.pcbi.1004794.

---
---
Our "constant velocity" simulations are pulling the molecule apart at a constant velocity, the same as this optical tweezer experiment. They use the hairpin system P5abc and its truncations, which can be simulated in the tutorial.

Mechanical unfolding of RNA hairpins: Liphardt, J. et al., **Reversible unfolding of single RNA molecules by mechanical force.** *Science*. (2001). https://doi.org/10.1126/science.1058498.

---
---
We can compare data from the small hairpin system P5GA to this paper.

Modeling the P5GA hairpin: Hyeon, C. and Thirumalai, D. **Mechanical unfolding of RNA hairpins.** *PNAS* (2005). https://doi.org/10.1073/pnas.0408314102.

---
---
Here, one thing to of note is that they construct PDB files by rearranging residues in a similar manner to how our PDB files for P5abc$\Delta$

Modeling P5abc: Koculi, E. et al., **Folding path of P5abc RNA involves direct coupling of secondary and tertiary structures.** *Nucleic Acids Research* (2012). https://doi.org/10.1093/nar/gks468

---