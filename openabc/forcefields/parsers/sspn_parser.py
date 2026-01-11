import numpy as np
import pandas as pd
import mdtraj
from openabc.utils import helper_functions
from openabc.utils.shadow_map import find_cg_pairs_from_atomistic_pdb

_smog_rna_nucleotide_mass_dict = dict(RA=135.13, RC=111.1, RG=151.13, RU=112.0868)

_smog_rna_nucleotide_charge_dict = dict(RA=-1.0, RC=-1.0, RG=-1.0, RU=-1.0)

class SSPNParser(object):
    """
    SSPN RNA parser.
    """
    def __init__(self, atomistic_pdb, cg_pdb, default_parse=True):
        """
        Initialize am RNA with SMOG model. 
        
        Parameters
        ----------
        atomistic_pdb : str
            Path for the atomistic pdb file. 
        
        cg_pdb : str
            path for the output pdb file. 

        default_parse : bool
            Whether to parse with default settings. 
        
        """
        
        self.atomistic_pdb = atomistic_pdb
        self.cg_pdb = cg_pdb
        self.atomistic_dataframe = helper_functions.get_dataframe_from_pdb(cg_pdb)
        if default_parse:
            print('Parse configuration with default settings.')
            self.parse_mol()

    @classmethod
    def from_atomistic_pdb(cls, atomistic_pdb, cg_pdb, cg_type="phosphorus", write_TER=False, default_parse=True):
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
        helper_functions.atomistic_pdb_to_p_pdb(atomistic_pdb, cg_pdb, cg_type, write_TER)
        return cls(atomistic_pdb, cg_pdb)
    
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
        self.sspn_exclusions = self.exclusions.copy()
        
    def parse_mol(self, get_native_pairs=True, frame=0, radius=0.1, bonded_radius=0.05, cutoff=0.6, box=None, 
                  use_pbc=False, exclude12=True, exclude13=True, exclude14=True, exclude_native_pairs=True, 
                  mass_dict=_smog_rna_nucleotide_mass_dict, charge_dict=_smog_rna_nucleotide_charge_dict, 
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

        # set bonds, angles, and dihedrals
        bonds, angles, dihedrals = [], [], []
        n_atoms = len(self.atomistic_dataframe.index)
        for atom1 in range(n_atoms):
            chain1 = self.atomistic_dataframe.loc[atom1, 'chainID']
            if atom1 < n_atoms - 1:
                atom2 = atom1 + 1
                chain2 = self.atomistic_dataframe.loc[atom2, 'chainID']
                if chain1 == chain2:
                    bonds.append([atom1, atom2])
            if atom1 < n_atoms - 2:
                atom3 = atom1 + 2
                chain3 = self.atomistic_dataframe.loc[atom3, 'chainID']
                if (chain1 == chain2) and (chain1 == chain3):
                    angles.append([atom1, atom2, atom3])
            if atom1 < n_atoms - 3:
                atom4 = atom1 + 3
                chain4 = self.atomistic_dataframe.loc[atom4, 'chainID']
                if (chain1 == chain2) and (chain1 == chain3) and (chain1 == chain4):
                    dihedrals.append([atom1, atom2, atom3, atom4])
        bonds, angles, dihedrals = np.array(bonds), np.array(angles), np.array(dihedrals)
        traj = mdtraj.load_pdb(self.cg_pdb)

        self.sspn_bonds = pd.DataFrame(bonds, columns=['a1', 'a2'])
        self.sspn_bonds['r0'] = mdtraj.compute_distances(traj, bonds, periodic=use_pbc)[frame]
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
            self.native_pairs = find_cg_pairs_from_atomistic_pdb(self.atomistic_pdb, frame, radius, bonded_radius, cutoff, 'rna', box, use_pbc)
            self.native_pairs.loc[:, 'epsilon_G'] = 1.0*bonded_energy_scale
            self.native_pairs.loc[:, 'sigma_G'] = 0.05
            self.native_pairs.loc[:, 'alpha_G'] = 1.6777216e-5*bonded_energy_scale
        # set exclusions
        self.parse_exclusions(exclude12, exclude13, exclude14, exclude_native_pairs) 
        # set mass and charge
        for i, row in self.atomistic_dataframe.iterrows():
            self.atomistic_dataframe.loc[i, 'mass'] = mass_dict[row['resname']]
            self.atomistic_dataframe.loc[i, 'charge'] = charge_dict[row['resname']]



