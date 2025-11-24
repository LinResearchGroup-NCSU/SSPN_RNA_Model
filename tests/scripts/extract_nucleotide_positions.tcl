# Load the DCD and PDB files
mol new 1hr2_simulation.pdb
mol addfile ../output.dcd waitfor all

# Select the frame to analyze
set frame 20
set sel [atomselect top "resname RA RU RC RG" frame $frame]

# Get the positions and write them to nucleotide_positions.txt
set positions [$sel get {x y z}]
set outfile [open "nucleotide_positions.txt" w]
puts $outfile [join $positions " "]
close $outfile

# Cleanup
$sel delete
puts "Nucleotide positions for frame $frame saved to nucleotide_positions.txt"
