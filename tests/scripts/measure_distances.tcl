# Get the filenames and residue IDs from command-line arguments
set pdb_file [lindex $argv 0]
set dcd_file [lindex $argv 1]
set resid1 [lindex $argv 2]
set resid2 [lindex $argv 3]

# Load the PDB and DCD files
mol new $pdb_file
mol addfile $dcd_file type dcd waitfor all

set sel_ref [atomselect top "backbone"]
measure fit $sel_ref $sel_ref

# Open a file to save distances
# Define the absolute path for the output file
set output_dir "C:/Users/Thomas/Documents/workspace/research/OpenABC_RNA/tests/scripts/"
set outfile [open "${output_dir}distances.txt" w]

# Select atoms for the two residues (e.g., first atom in each residue)
set sel1 [atomselect top "resid $resid1 and name P"]
set sel2 [atomselect top "resid $resid2 and name P"]

# Check if the selections are valid
if {[$sel1 num] == 0 || [$sel2 num] == 0} {
    puts "Error: Residue selection failed. Check residue IDs or atom names."
    exit
}

# Initialize an empty list to store distances
set distance_list ""

# Loop through all frames in the trajectory
set num_frames [molinfo top get numframes]
for {set frame 0} {$frame < $num_frames} {incr frame} {
    # Go to the current frame
    animate goto $frame

    # Get the coordinates of the selected atoms
    set coord1 [lindex [$sel1 get {x y z}] 0]
    set coord2 [lindex [$sel2 get {x y z}] 0]

    # Measure the distance
    set distance [veclength [vecsub $coord1 $coord2]]

    # Append the distance to the list
    append distance_list "$distance "
}

# Write the concatenated distances to the file
puts $outfile $distance_list

# Cleanup and close file
close $outfile
$sel1 delete
$sel2 delete

puts "Distances saved to distances.txt"

# Debugging: Print arguments to verify they're passed correctly
puts "PDB file: $pdb_file"
puts "DCD file: $dcd_file"
puts "Residue 1: $resid1"
puts "Residue 2: $resid2"

# Debugging: Print number of atoms in selections
puts "Number of atoms in resid $resid1: [$sel1 num]"
puts "Number of atoms in resid $resid2: [$sel2 num]"

# Debugging: Check frames
puts "Number of frames: $num_frames"