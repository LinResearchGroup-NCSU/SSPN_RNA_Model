import os
import math
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import subprocess

which_pdb = "1f9l_cg.pdb"
notebook_path = os.path.abspath('__file__')
notebook_dir = os.path.dirname(notebook_path)

def display_in_VMD(batch_id="batch001", run_id="001"):
    tcl_script = notebook_dir + "/scripts/graphical_settings.tcl"

    batch_dir = os.path.join(notebook_dir, batch_id)
    structure_file_path = os.path.join(batch_dir, run_id, which_pdb)
    trajectory_file_path = os.path.join(batch_dir, run_id, 'trajectory.dcd')
    command = f'vmd.exe -e "{tcl_script}" "{structure_file_path}" "{trajectory_file_path}"'
    output = subprocess.Popen(command)


def calculate_energy_drift(csv_path):
    """
    Given the path to an energy_log.csv file, compute the average relative
    energy drift |delta E|, ignoring the first 0.5 ns entirely, using E0
    as the average total energy from [0.5 ns, 1.0 ns], and measuring drift
    after 1.0 ns.
    """
    df = pd.read_csv(csv_path, skiprows=0)

    # Convert ps to ns
    df["Time (ns)"] = df["Time (ps)"] / 1000.0
    # Define E0 as the average total energy from 0.5 ns to 1.0 ns
    eq_data = df[(df["Time (ns)"] >= 0.5) & (df["Time (ns)"] <= 1.0)]
    if len(eq_data) == 0:
        # If there's no data in that slice, we can't define E0 properly
        return float('nan')
    E0 = eq_data["Total Energy (kJ/mole)"].mean()
    # Only consider frames AFTER 1.0 ns for measuring drift
    production = df[df["Time (ns)"] > 1.0]
    if len(production) == 0:
        # No production data to analyze
        return float('nan')
    total_energy_prod = production["Total Energy (kJ/mole)"]
    # Compute mean of |(Ek - E0)/E0|
    rel_drift = abs((total_energy_prod - E0) / E0).mean()
    return rel_drift

def read_csv_column(file_path, column_index):
    values = []
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # skip header row
        for row in reader:
            if len(row) >= 3:
                values.append(float(row[column_index]))
    return values

def energy_drift_analysis(batch_id, outliers = []):
    # Directory containing run subfolders

    # Lists to store our data
    timesteps_fs = []
    log_drift_vals = []

    batch_parameters_file_path = os.path.join(batch_id,"batch_parameters.csv")
    timesteps_fs = read_csv_column(batch_parameters_file_path, 2)
    
    #this bit of code is incase there are any broken simulations
    for i in sorted(outliers, reverse=True):
        del timesteps_fs[i-1]

    # Loop over each run directory
    run_names = sorted(os.listdir(batch_id))
    for i in range(0,len(run_names)):
        run_name = run_names[i]
        run_path = os.path.join(batch_id, run_name)
        if (i+1) in outliers: continue # this way we don't plot the outliers
        if os.path.isdir(run_path) and run_name.isdigit():

            energy_log_file_path = os.path.join(run_path, "energy_log.csv")
            # Calculate drift value
            drift_value = calculate_energy_drift(energy_log_file_path)

            # Take ln(|delta E|)
            if drift_value > 0:
                log_drift = math.log(drift_value)
            else:
                log_drift = float('nan')

            log_drift_vals.append(log_drift)

    # --- PLOT data and add a trendline ---
    plt.figure()

    # Scatter-plot the data
    plt.plot(timesteps_fs, log_drift_vals, 'o')

    # Fit a straight line y = m x + b
    # (Handle cases where you might have NaN by filtering them out.)
    valid_indices = [i for i in range(len(timesteps_fs)) 
                    if not np.isnan(log_drift_vals[i])]
    if valid_indices:
        x_valid = np.array([timesteps_fs[i] for i in valid_indices])
        y_valid = np.array([log_drift_vals[i] for i in valid_indices])

        # Polyfit for a line
        coeffs = np.polyfit(x_valid, y_valid, deg=1)  # slope, intercept
        poly_fn = np.poly1d(coeffs)  # callable function, e.g. poly_fn(x)

        # Generate a smooth x-range for plotting the line
        x_smooth = np.linspace(min(x_valid), max(x_valid), 100)
        y_smooth = poly_fn(x_smooth)

        plt.plot(x_smooth, y_smooth, '--')

    plt.tick_params(axis='both', labelsize=18)
    plt.xlabel('Timestep (fs)', fontsize=16)
    plt.ylabel('ln(|ΔE|)', fontsize=16)
    plt.title('ln(|ΔE|) vs. Timestep (fs)', fontsize=24)
    plt.legend()
    plt.show()

def EED_analysis(batch_file_path, run_id, temperature, bins,offset):
    structure_pdb_file_path = os.path.join(batch_file_path, run_id, which_pdb)
    trajectory_file_path = os.path.join(batch_file_path, run_id, 'trajectory.dcd')
    traj = mdtraj.load(trajectory_file_path, top=structure_pdb_file_path)

    resid1 = 0
    resid2 = traj.n_residues - 1

    equilibrium_trajectory = traj

    EED_for_each_frame = mdtraj.compute_distances(equilibrium_trajectory, [[resid1, resid2]])[:, 0]
    hist_counts, bin_edges = np.histogram(EED_for_each_frame, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Compute free energy: F(EED) = -kB T ln P(EED)
    probabilities = hist_counts / np.sum(hist_counts)  # Normalize to probabilities
    kB_kJmol = 0.008314  # in kJ/(mol*K)
    kB_kCalmol = kB_kJmol * (1/4.184)
    #free_energy = -kB_kJmol * temperature * np.log(probabilities)
    free_energy = -kB_kCalmol * temperature * np.log(probabilities)
    
    # Handle log(0) cases by shifting to max energy value
    #free_energy -= np.nanmin(free_energy[np.isfinite(free_energy)])  # Normalize lowest energy to 0
    free_energy = free_energy + offset
    return EED_for_each_frame, bin_centers, bin_edges, hist_counts, free_energy

def calculate_Q_N_by_frame(traj, native_positions, native_pairs_df, frame_number, lambda_factor=1.2, r_0=0.5):
    frame_positions = traj.xyz[frame_number]
    r0ij_list = np.linalg.norm(native_positions[native_pairs_df['a1'].astype(int).values] -
                               native_positions[native_pairs_df['a2'].astype(int).values], axis=1)
    d0_list = lambda_factor * r0ij_list
    rij_list = np.linalg.norm(frame_positions[native_pairs_df['a1'].astype(int).values] -
                              frame_positions[native_pairs_df['a2'].astype(int).values], axis=1)
    s_r = np.where(rij_list < d0_list, 1, np.exp(-((rij_list - d0_list) ** 2) / (2 * r_0 ** 2)))
    return np.sum(s_r)/len(native_pairs_df)

def Q_analysis(batch_file_path, run_id, temperature, native_pairs_df, bins, lambda_factor = 1.2, r_0=0.5):
    structure_pdb_file_path = os.path.join(batch_file_path, run_id, which_pdb)
    trajectory_file_path = os.path.join(batch_file_path, run_id, 'trajectory.dcd')

    traj = mdtraj.load(trajectory_file_path, top=structure_pdb_file_path)
    native_traj = mdtraj.load(structure_pdb_file_path)
    native_positions = native_traj.xyz[0]
    Q_for_each_frame = [calculate_Q_N_by_frame(traj, native_positions, native_pairs_df, frame, lambda_factor, r_0 = r_0) for frame in range(traj.n_frames)]

    hist_counts, bin_edges = np.histogram(Q_for_each_frame, bins=bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    kB_kJmol = 0.008314
    probabilities = hist_counts / np.sum(hist_counts)  # Normalize to probabilities
    free_energy = -kB_kJmol * temperature * np.log(probabilities)
    free_energy -= np.nanmin(free_energy[np.isfinite(free_energy)])  # Normalize lowest energy to 0
    return Q_for_each_frame, bin_centers, bin_edges, hist_counts, free_energy

def plot_collective_variable_over_simulation(collective_variable_for_each_frame, collective_variable_type, temperature):
    plt.plot(range(len(collective_variable_for_each_frame)), collective_variable_for_each_frame)
    plt.title(str(temperature) + "K - " + collective_variable_type +" over time")
    plt.xlabel("Frame")
    plt.ylabel(collective_variable_type)

def plot_collective_variable_histogram(bin_centers, hist_counts, bin_edges, collective_variable_type, temperature):
    plt.bar(bin_centers, hist_counts, width=np.diff(bin_edges), color='blue', alpha=0.7, edgecolor='black')
    plt.title(str(temperature) + "K - " + collective_variable_type + " Histogram")
    plt.xlabel(collective_variable_type, fontsize=18)
    plt.ylabel("Probability Density", fontsize=18)

def plot_collective_variable_free_energy(bin_centers, free_energy, collective_variable_type, temperature, color="red", label=""):
    line = plt.plot(bin_centers, free_energy, linestyle='-', color=color, label=label)
    # Remove this line – you are setting the title globally later anyway
    # plt.title(str(temperature) + "K - " + collective_variable_type + " Free Energy")
    plt.xlabel("R [nm]", fontsize=18)
    plt.ylabel("F(R) (kCal/mol)", fontsize=18)
    return line

# Putting these here in case I want to use them later. . .

# def compute_heat_capacity(batch_id, run_id):
#     # Construct the path to the energy log
#     path = os.path.join(batch_id, run_id, 'energy_log.csv')
    
#     # Read the CSV file, ignoring the comment character (#)
#     try:
#         data = pd.read_csv(path, comment='#')
#     except FileNotFoundError:
#         raise FileNotFoundError(f"Could not find energy log at {path}")
    
#     # Extract the necessary columns using indices directly
#     temperatures = data.iloc[:, 5]  # 6th column is Temperature
#     total_energies = data.iloc[:, 2]  # 5th column is Total Energy
    
#     # Compute averages
#     avg_energy = total_energies.mean()
#     avg_energy_squared = (total_energies ** 2).mean()
    
#     # Compute the variance of the energy
#     variance_energy = avg_energy_squared - avg_energy ** 2
    
#     # Boltzmann constant in kJ/mol K
#     k_B = 0.0083145
    
#     # Compute the heat capacity
#     avg_temp = temperatures.mean()
#     heat_capacity = variance_energy / (k_B * avg_temp ** 2)
    
#     print(f"Heat capacity for {batch_id}/{run_id}: {heat_capacity} kJ/mol K")
#     return heat_capacity



# def calculate_total_Q_N_average(run_id, lambda_factor=1.2, r_0=0.05):
#     structure_pdb_file_path = os.path.join(runs_file_path, run_id, '1f9l_cg.pdb')
#     trajectory_file_path = os.path.join(runs_file_path, run_id, 'trajectory.dcd')
#     traj = mdtraj.load(trajectory_file_path, top=structure_pdb_file_path)
#     native_traj = mdtraj.load(structure_pdb_file_path)
#     native_positions = native_traj.xyz[0]
#     native_pairs = rna_model.native_pairs[['a1', 'a2']].to_numpy(dtype=int)

#     r0ij_list = np.linalg.norm(native_positions[native_pairs[:, 0]] - native_positions[native_pairs[:, 1]], axis=1)
#     d0_list = lambda_factor * r0ij_list
#     Q_list = np.zeros(traj.n_frames)

#     for frame_number in range(traj.n_frames):
#         frame_positions = traj.xyz[frame_number]
#         rij_list = np.linalg.norm(frame_positions[native_pairs[:, 0]] - frame_positions[native_pairs[:, 1]], axis=1)
#         s_r = np.where(rij_list < d0_list, 1, np.exp(-((rij_list - d0_list) ** 2) / (2 * r_0 ** 2)))
#         Q_list[frame_number] = np.mean(s_r)
#     return np.mean(Q_list)

# def find_heat_capacity(run_id, temperature):
#     # Load trajectory and native structure
#     structure_pdb_file_path = os.path.join(runs_file_path, run_id, '1f9l_cg.pdb')
#     trajectory_file_path = os.path.join(runs_file_path, run_id, 'trajectory.dcd')

#     traj = mdtraj.load(trajectory_file_path, top=structure_pdb_file_path).xyz
#     native_conformation = mdtraj.load(structure_pdb_file_path).xyz

#     # Combine native conformation with trajectory frames
#     conformations = np.vstack((native_conformation, traj))  # Faster than np.concatenate
#     n_frames = conformations.shape[0]
#     simulation = rna_model.simulation
#     context = simulation.context
#     energies = np.zeros(n_frames)

#     for i in range(n_frames):
#         context.setPositions(conformations[i])
#         state = context.getState(getEnergy=True)
#         energies[i] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

#     k_B = 1  # Assuming temperature is in appropriate units
#     sigma_squared = np.var(energies)  # Faster than np.average calculations
#     C_v = sigma_squared / (k_B * temperature)

#     return C_v

# def find_energies(run_id="001"):
#     trajectory_file_path = os.path.join(runs_file_path,run_id,'trajectory.dcd')
#     structure_pdb_file_path = os.path.join(runs_file_path, run_id, which_pdb)
#     dcd_conformations = mdtraj.load(trajectory_file_path, top=structure_pdb_file_path).xyz
#     native_conformation = mdtraj.load(structure_pdb_file_path).xyz

#     conformations = np.concatenate((native_conformation, dcd_conformations), axis = 0)
#     n_frames = len(conformations)

#     simulation = rna_model.simulation
#     openmm_energies = []
#     for i in range(n_frames):
#         row = [i]
#         simulation.context.setPositions(conformations[i])
#         for j in force_groups:
#             state = rna_model.simulation.context.getState(getEnergy=True, groups={j})
#             energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
#             row.append(energy)
#         state = simulation.context.getState(getEnergy=True)
#         energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
#         row.append(energy)
#         openmm_energies.append(row)

#     openmm_energies = np.array(openmm_energies)
#     columns = ['frame', 'rna bond', 'rna angle', 'rna dihedral', 'native pair', 'elec', 'vdwl', 'sum']
#     df_openmm_energies = pd.DataFrame(openmm_energies, columns=columns).round(6)
#     df_openmm_energies.to_excel(energies_file_path, index=False)

# import os
# import numpy as np
# import matplotlib.pyplot as plt

# plt.figure(figsize=(6, 3), dpi= 80, facecolor='w', edgecolor='k')

# plt.subplot(1, 2, 1)

# def my_le_range(start, end, step):
#     while start <= end:
#         yield start
#         start += step


# data = []
# for i in my_le_range(3000, 3000, 50):
#     temperature  = i / 10
     
#     for name in ['aa', 'ab', 'ac']:
#         folder = "c:/Users/Thomas/Documents/workspace/research/OpenABC_RNA/tests/result//wham1d." + str(name) + '/free'
        
#         os.chdir(folder)
#         if name == 'aa':
#             data = np.loadtxt(str(i))[:,:2].T
#         else:            
#             data_tmp = np.loadtxt(str(i))[:,1]
#             data = np.vstack((data, data_tmp))
    
#     data = data.T
#     print(data)

#     data_x = data[:,0]
#     # Add entropic penalty for cartesian-spherical coordinate Jacobian
#     data_y = np.average(data[:,1:]/(temperature * 0.008314), axis=1)
#     error_y = np.std(data[:,1:]/(temperature * 0.008314), axis=1)  

#     # Shift so that the minimum of y value is zero
#     data_y = data_y - np.amin(data_y)
#     plt.plot(data_x, data_y, 'k-')
#     plt.fill_between(data_x, data_y-error_y, data_y+error_y, alpha=0.3, facecolor='k')

# plt.xlabel('E2E Distance (nm)')
# plt.ylabel('PMF (kT)')  


# plt.tight_layout()
# #plt.savefig('/Users/xl23/Dropbox_work/Ongoing_projects/PRC2_DNA_interaction/figure/binding_test.pdf')


# # Clip Trajectory

# import mdtraj as md

# file_path = "batch043/001/"
# # Load the trajectory
# traj = md.load(file_path+'trajectory.dcd', top=file_path+'P5abc_cg.pdb')  # Or your .psf/.gro file

# # Clip the trajectory to only include frames up to frame N (e.g., frame 100)
# clipped_traj = traj[:7521]  # Includes frames 0 to 100

# # Save the new clipped trajectory
# clipped_traj.save(file_path+'clipped_trajectory.dcd')