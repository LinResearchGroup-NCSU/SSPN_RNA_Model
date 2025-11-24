import numpy as np
import pandas as pd
try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit
from openabc.forcefields.cg_model import CGModel
from openabc.forcefields import functional_terms
from openabc.forcefields.parameters import Mixin3SPN2ConfigParser
from simtk.openmm import CustomCentroidBondForce
from openmm import CustomBondForce
import warnings
import os

__location__ = os.path.dirname(os.path.abspath(__file__))

class SMOGSSPNModel(CGModel, Mixin3SPN2ConfigParser):
    """
    The class for SMOG+3SPN2 model. 
    To ensure this model works properly, please ensure two neighboring ssDNA chains do not share same chainID. 
    """
    def __init__(self, default_parse_config=True):
        """
        Initialize. 
        
        Parameters
        ----------
        dna_type : str
            DNA type. This is related to force field parameters. 
        
        default_parse_config : bool
            Whether to parse the default 3SPN2 configuration file. 
        
        """
        self.atoms = None
        self.exclusions = None
        self.bonded_attr_names = ['sspn_bonds', 'sspn_angles', 'sspn_dihedrals', 'native_pairs', 'sspn_exclusions']
        if default_parse_config:
            # load parameters
            self.parse_config_file()
    
    def append_mol(self, new_mol, verbose=False):
        """
        The method can append new molecules by concatenating atoms and bonded interaction information saved in dataframes. 
        Please ensure two neighboring chains do not share chainID. 
        
        Parameters
        ----------
        new_mol : a consistent parser object or CG model object
            The object of a new molecule including atom and bonded interaction information. 
        
        verbose : bool
            Whether to report the appended attributes. 
        
        """
        super().append_mol(new_mol, verbose)

    def add_sspn_bonds(self, force_group=1):
        """
        Add sspn bonds.
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        """
        if hasattr(self, 'sspn_bonds'):
            print('Add sspn bonds.')
            force = functional_terms.harmonic_bond_term(self.sspn_bonds, self.use_pbc, force_group)
            self.system.addForce(force)

    def add_sspn_angles(self, threshold=130*np.pi/180, clip=False, force_group=2, verbose=True):
        """
        Add sspn angles. 
        
        Parameters
        ----------
        threshold : float
            The angle for which greater angles will be clipped, if desired. Designed to reduce energy spikes.

        clip : bool
            Whether or not to use angle clipping

        force_group : int
            Force group. 
        
        """
        if hasattr(self, 'sspn_angles'):
            any_theta0_beyond_threshold = False
            for i, row in self.sspn_angles.iterrows():
                a1 = int(row['a1'])
                a2 = int(row['a2'])
                a3 = int(row['a3'])
                theta0 = float(row['theta0'])
                if (theta0 > threshold) and verbose:
                    warnings.warn(f'Warning: angle composed of atom ({a1}, {a2}, {a3}) has theta0 equal to {theta0}, which is larger than the threshold value equal to {threshold}!')
                    any_theta0_beyond_threshold = True
            if clip and any_theta0_beyond_threshold:
                print(f'Decrease all the theta0 values larger than {threshold} to {threshold}.')
                self.sspn_angles.loc[self.sspn_angles['theta0'] > threshold, 'theta0'] = threshold
            force = functional_terms.harmonic_angle_term(self.sspn_angles, self.use_pbc, force_group)
            self.system.addForce(force)

    def add_sspn_dihedrals(self, force_group=3):
        """
        Add sspn dihedrals. 
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        """
        if hasattr(self, 'sspn_dihedrals'):
            print('Add sspn dihedrals.')
            force = functional_terms.periodic_dihedral_term(self.sspn_dihedrals, self.use_pbc, force_group)
            self.system.addForce(force)

    def add_native_pairs(self, force_group=4):
        """
        Add native pairs. 
        
        Parameters
        ----------
        force_group : int
            Force group.
        
        """
        if hasattr(self, 'native_pairs'):
            print('Add native pairs.')
            force = functional_terms.native_pair_gaussian_term(self.native_pairs, self.use_pbc, force_group)
            self.system.addForce(force)
   
    def add_constant_force(self, constant_force, force_group=13):
        """
        Add constant force between ends. 
        
        Parameters
        ----------
        constant_force : float
            Magnitude of force desired

        force_group : int
            Force group.
        
        """
        constant_force_converted = constant_force.value_in_unit(unit.kilojoule / unit.nanometer) * unit.AVOGADRO_CONSTANT_NA
        constant_force_val = constant_force_converted.value_in_unit(constant_force_val.unit)
        first_residue_index, last_residue_index = self.atoms.index[0], self.atoms.index[-1]
        force = CustomCentroidBondForce(2, "force * distance(g1, g2)")
        force.addPerBondParameter("force")
        g1 = force.addGroup([first_residue_index])
        g2 = force.addGroup([last_residue_index])
        force.addBond([g1, g2], [-1*constant_force_val], force_group)
        print('Add constant force.')
        self.system.addForce(force)

    def add_force_extension(self, k, force_group=14):
        """
        Add constant force between ends. 
        
        Parameters
        ----------
        k : float
            Spring constant for simulated optical traps at each end of the system
            
        force_group : int
            Force group.
        
        """
        first_residue_index, last_residue_index = self.atoms.index[0], self.atoms.index[-1]
        k_val = k.in_units_of(unit.kilojoule / unit.nanometer **2) * unit.AVOGADRO_CONSTANT_NA
        force = CustomBondForce("0.5*k*(r-r0)^2")
        force.addPerBondParameter("k")
        force.addGlobalParameter("r0", 0.0)
        force.addBond(first_residue_index, last_residue_index, [k_val])
        force.setUsesPeriodicBoundaryConditions(True)
        force.setForceGroup(force_group)
        print('Add force-extension.')
        self.system.addForce(force)
        
    def add_smog_vdwl(self, cutoff=1.6*unit.nanometer, force_group=11):
        """
        Add SMOG nonbonded Van der Waals interactions. 
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        """
        print('Add all nonbonded contact interactions.')
        force = functional_terms.all_smog_NC_term(self, cutoff, force_group)
        self.system.addForce(force)

    def add_elec_switch(self, salt_conc=150.0*unit.millimolar, temperature=300.0*unit.kelvin, 
                        cutoff1=1.2*unit.nanometer, cutoff2=1.5*unit.nanometer, switch_coeff=[1, 0, 0, -10, 15, -6], 
                        add_native_pair_elec=True, force_group=6):
        """
        Add electrostatic interaction with switch function. 
        
        The switch function switches potential to zero within range cutoff1 < r <= cutoff2. 
        
        Parameters
        ----------
        salt_conc : Quantity
            Monovalent salt concentration. 
        
        temperature : Quantity
            Temperature. 
        
        cutoff1 : Quantity
            Cutoff distance 1. 
        
        cutoff2 : Quantity
            Cutoff distance 2. 
        
        switch_coeff : sequence-like
            Switch function coefficients. 
        
        add_native_pair_elec : bool
            Whether to add electrostatic interactions between native pairs. 
        
        force_group : int
            Force group. 
        
        """
        print('Add electrostatic interactions with distance-dependent dielectric and switch.')
        charges = self.atoms['charge'].tolist()
        force1 = functional_terms.ddd_dh_elec_switch_term(charges, self.exclusions, self.use_pbc, salt_conc, 
                                                          temperature, cutoff1, cutoff2, switch_coeff, force_group)
        self.system.addForce(force1)
        if add_native_pair_elec and hasattr(self, 'native_pairs'):
            print('Add electrostatic interactions between native pair atoms.')
            df_charge_bonds = pd.DataFrame(columns=['a1', 'a2', 'q1', 'q2'])
            for i, row in self.native_pairs.iterrows():
                a1, a2 = int(row['a1']), int(row['a2'])
                q1, q2 = float(charges[a1]), float(charges[a2])
                if (q1 != 0) and (q2 != 0):
                    df_charge_bonds.loc[len(df_charge_bonds.index)] = [a1, a2, q1, q2]
            force2 = functional_terms.ddd_dh_elec_switch_bond_term(df_charge_bonds, self.use_pbc, salt_conc, 
                                                                   temperature, cutoff1, cutoff2, switch_coeff, 
                                                                   force_group)
            self.system.addForce(force2)
        else:
            print('Do not add electrostatic interactions between native pair atoms.')   