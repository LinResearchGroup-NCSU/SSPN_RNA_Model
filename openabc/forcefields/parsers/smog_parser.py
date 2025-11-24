import numpy as np
import pandas as pd
import mdtraj
from openabc.utils import helper_functions
from openabc.utils.shadow_map import find_cg_pairs_from_atomistic_pdb
from openabc.lib import _amino_acids
from openabc.lib import _rna_nucleotides

_smog_amino_acid_mass_dict = dict(ALA=71.0788, ARG=156.1875, ASN=114.1038, ASP=115.0886, CYS=103.1388, 
                                  GLU=129.1155, GLN=128.1307, GLY=57.0519, HIS=137.1411, ILE=113.1594, 
                                  LEU=113.1594, LYS=128.1741, MET=131.1926, PHE=147.1766, PRO=97.1167, 
                                  SER=87.0782, THR=101.1051, TRP=186.2132, TYR=163.1760, VAL=99.1326)

_smog_amino_acid_charge_dict = dict(ALA=0.0, ARG=1.0, ASN=0.0, ASP=-1.0, CYS=0.0, 
                                    GLN=0.0, GLU=-1.0, GLY=0.0, HIS=0.0, ILE=0.0,
                                    LEU=0.0, LYS=1.0, MET=0.0, PHE=0.0, PRO=0.0,
                                    SER=0.0, THR=0.0, TRP=0.0, TYR=0.0, VAL=0.0)

_smog_rna_nucleotide_mass_dict = dict(RA=135.13, RC=111.1, RG=151.13, RU=112.0868)
#_smog_rna_nucleotide_mass_dict = dict(RA=135.13*100, RC=111.1*100, RG=151.13*100, RU=112.0868*100)

_smog_rna_nucleotide_charge_dict = dict(RA=-1.0, RC=-1.0, RG=-1.0, RU=-1.0)

#_smog_rna_nucleotide_mass_dict = dict(A=135.13, C=111.1, G=151.13, U=112.0868)

#_smog_rna_nucleotide_charge_dict = dict(A=-1.0, C=-1.0, G=-1.0, U=-1.0)

class SMOGParser(object):
    """
    SMOG protein/RNA parser.
    """
    def __init__(self, atomistic_pdb, cg_pdb, cut_pdb=None, mol_type='protein', default_parse=True):
        """
        Initialize a protein/rna with SMOG model. 
        
        Parameters
        ----------
        atomistic_pdb : str
            Path for the atomistic pdb file. 
        
        parsed_pdb : str
            path for the output pdb file. 
        
        mol_type : str
            Accepts 'protein' or 'rna'. Whether to parse pdb as protein (CA atoms) or RNA (phosphorus atoms).

        default_parse : bool
            Whether to parse with default settings. 
        
        """
        
        self.mol_type = mol_type

        # Verify the input has correct type parameter (I didn't do this in any other function. . . might be a good idea to do that later.)
        if self.mol_type != 'protein' and self.mol_type != 'rna':
            raise ValueError("Invalid mol_type. Use 'protein' or 'rna'.")
        
        self.atomistic_pdb = atomistic_pdb
        self.pdb = cg_pdb
        self.cut_pdb = cut_pdb
        self.atoms = helper_functions.parse_pdb(cg_pdb, mol_type, write_output=False, unaltered=False)
        if self.mol_type == 'protein':
            # check if all the atoms are protein CA atoms
            assert ((self.atoms['resname'].isin(_amino_acids)).all() and self.atoms['name'].eq('CA').all())
            if default_parse:
                print('Parse configuration with default settings.')
                self.parse_mol()
        elif self.mol_type == 'rna':
            if default_parse:
                print('Parse configuration with default settings.')
                self.parse_mol()

    @classmethod
    def from_atomistic_pdb(cls, atomistic_pdb, cg_pdb, cut_pdb, mol_type='protein', cg_type='phosphorus', write_TER=False, default_parse=True):
        """
        Initialize an SMOG model protein/RNA from atomistic pdb. 
        
        Parameters
        ----------
        atomistic_pdb : str
            Path for the input atomistic pdb file. 
        
        parsed_pdb : str
            Path for the output CA/phosphorus pdb file. 
            
        mol_type : str
            Accepts 'protein' or 'rna'. Whether to parse pdb as protein (CA atoms) or RNA (phosphorus atoms) 
        
        write_TER : bool
            Whether to write TER between two chains. 
        
        default_parse : bool
            Whether to parse with default settings. 
        
        """
        helper_functions.atomistic_pdb_to_cg_pdb(atomistic_pdb, cg_pdb, cut_pdb, mol_type, cg_type, write_TER)
        return cls(atomistic_pdb, cg_pdb, cut_pdb, mol_type, default_parse)
    
    def parse_exclusions(self, exclude12=True, exclude13=True, exclude14=True, exclude_native_pairs=True):
        """
        Parse nonbonded exclusions based on bonds, angles, dihedrals, and native pairs.
        
        Parameters
        ----------
        exclude12 : bool
            Whether to exclude nonbonded interactions between 1-2 atom pairs. 
        
        exclude13 : bool
            Whether to exclude nonbonded interactions between 1-3 atom pairs. 
        
        exclude14 : bool
            Whether to exclude nonbonded interactions between 1-4 atom pairs. 
        
        exclude_native_pairs : bool
            Whether to exclude nonbonded interactions between native pairs. 
        
        """
        exclusions = []
        if exclude12 and hasattr(self, 'sspn_bonds'):
            for i, row in self.sspn_bonds.iterrows():
                exclusions.append((int(row['a1']), int(row['a2'])))
        if exclude13 and hasattr(self, 'sspn_angles'):
            for i, row in self.sspn_angles.iterrows():
                exclusions.append((int(row['a1']), int(row['a3'])))
        if exclude14 and hasattr(self, 'sspn_dihedrals'):
            for i, row in self.sspn_dihedrals.iterrows():
                exclusions.append((int(row['a1']), int(row['a4'])))
        if exclude_native_pairs and hasattr(self, 'native_pairs'):
            for i, row in self.native_pairs.iterrows():
                exclusions.append((int(row['a1']), int(row['a2'])))
        exclusions = np.array(sorted(exclusions))
        self.exclusions = pd.DataFrame(exclusions, columns=['a1', 'a2']).drop_duplicates(ignore_index=True)
        # for using with 3SPN2, also set self.protein_exclusions
        # thomas - 'protein' changed to sspn because this also works with RNA. Both proteins and RNA are single-site-per-residue (sspn) models.
        self.sspn_exclusions = self.exclusions.copy()
        
    def parse_mol(self, mol_type='protein', get_native_pairs=True, frame=0, radius=0.1, bonded_radius=0.05, cutoff=0.6, box=None, 
                  use_pbc=False, exclude12=True, exclude13=True, exclude14=True, exclude_native_pairs=True, 
                  mass_dict=_smog_amino_acid_mass_dict, charge_dict=_smog_amino_acid_charge_dict, 
                  bonded_energy_scale=2.5):
        """
        Parse molecule.
        
        Parameters
        ----------
        get_native_pairs : bool
            Whether to get native pairs from atomistic pdb with shadow algorithm. 
        
        frame : int
            Frame index to compute native configuration parameters for both CG and atomistic models. 
            
        radius : float or int
            Shadow algorithm radius. 
        
        bonded_radius : float or int
            Shadow algorithm bonded radius. 
        
        cutoff : float or int
            Shadow algorithm cutoff. 
        
        box : None or np.ndarray
        Specify box shape. 
        Note `use_pbc` = False is only compatible with orthogonal box. 
        If `use_pbc` = False, then set `box` as None. 
        If `use_pbc` = True, then `box` = np.array([lx, ly, lz, alpha1, alpha2, alpha3]). 
        
        use_pbc : bool
            Whether to use periodic boundary condition (PBC) for searching native pairs. 
        
        exclude12 : bool
            Whether to exclude nonbonded interactions between 1-2 atom pairs. 
        
        exclude13 : bool
            Whether to exclude nonbonded interactions between 1-3 atom pairs. 
        
        exclude14 : bool
            Whether to exclude nonbonded interactions between 1-4 atom pairs. 
        
        exclude_native_pairs : bool
            Whether to exclude nonbonded interactions between native pairs. 
        
        mass_dict : dict
            Mass dictionary. 
        
        charge_dict : dict
            Charge dictionary. 
        
        bonded_energy_scale : float or int
            Bonded energy additional scale factor. 
        """
        if(mol_type=='rna'):
            mass_dict = _smog_rna_nucleotide_mass_dict
            charge_dict = _smog_rna_nucleotide_charge_dict

        # set bonds, angles, and dihedrals
        bonds, angles, dihedrals = [], [], []
        n_atoms = len(self.atoms.index)
        for atom1 in range(n_atoms):
            chain1 = self.atoms.loc[atom1, 'chainID']
            if atom1 < n_atoms - 1:
                atom2 = atom1 + 1
                chain2 = self.atoms.loc[atom2, 'chainID']
                if chain1 == chain2:
                    bonds.append([atom1, atom2])
            if atom1 < n_atoms - 2:
                atom3 = atom1 + 2
                chain3 = self.atoms.loc[atom3, 'chainID']
                if (chain1 == chain2) and (chain1 == chain3):
                    angles.append([atom1, atom2, atom3])
            if atom1 < n_atoms - 3:
                atom4 = atom1 + 3
                chain4 = self.atoms.loc[atom4, 'chainID']
                if (chain1 == chain2) and (chain1 == chain3) and (chain1 == chain4):
                    dihedrals.append([atom1, atom2, atom3, atom4])
        bonds, angles, dihedrals = np.array(bonds), np.array(angles), np.array(dihedrals)
        traj = mdtraj.load_pdb(self.pdb)

        #'protein' replaced with 'sspn', standing for single site per residue. This is to denote that this works for our cg models of proteins and RNA, whereas DNA is 3spn (three site per nucleotide).
        self.sspn_bonds = pd.DataFrame(bonds, columns=['a1', 'a2'])
        self.sspn_bonds['r0'] = mdtraj.compute_distances(traj, bonds, periodic=use_pbc)[frame]
        # # #thomas - 20000 may be incorrect for RNA
        self.sspn_bonds.loc[:, 'k_bond'] = 20000*bonded_energy_scale
        self.sspn_angles = pd.DataFrame(angles, columns=['a1', 'a2', 'a3'])
        self.sspn_angles['theta0'] = mdtraj.compute_angles(traj, angles, periodic=use_pbc)[frame]
        self.sspn_angles.loc[:, 'k_angle'] = 40*bonded_energy_scale
        self.sspn_dihedrals = pd.DataFrame(columns=['a1', 'a2', 'a3', 'a4', 'periodicity', 'phi0', 'k_dihedral'])
        phi = mdtraj.compute_dihedrals(traj, dihedrals, periodic=use_pbc)[frame]
        for i in range(dihedrals.shape[0]):
            row = dihedrals[i].tolist() + [1, phi[i] + np.pi, 1.0*bonded_energy_scale]
            self.sspn_dihedrals.loc[len(self.sspn_dihedrals.index)] = row
            row = dihedrals[i].tolist() + [3, 3*(phi[i] + np.pi), 0.5*bonded_energy_scale]
            self.sspn_dihedrals.loc[len(self.sspn_dihedrals.index)] = row
        # set native pairs
        if get_native_pairs:
            print('Get native pairs with shadow algorithm.')
            self.native_pairs = find_cg_pairs_from_atomistic_pdb(self.atomistic_pdb, self.cut_pdb, frame, radius, bonded_radius, cutoff, mol_type, box, use_pbc)
            self.native_pairs.loc[:, 'epsilon_G'] = 1.0*bonded_energy_scale
            self.native_pairs.loc[:, 'sigma_G'] = 0.05
            self.native_pairs.loc[:, 'alpha_G'] = 1.6777216e-5*bonded_energy_scale
        # set exclusions
        self.parse_exclusions(exclude12, exclude13, exclude14, exclude_native_pairs) 
        # set mass and charge
        for i, row in self.atoms.iterrows():
            self.atoms.loc[i, 'mass'] = mass_dict[row['resname']]
            self.atoms.loc[i, 'charge'] = charge_dict[row['resname']]



