# Load the structure file
mol new 7qr4_simulation.pdb

# Select RNA nucleotides with custom residue names
set sel [atomselect top "resname RA RU RC RG"]
if {[$sel num] == 0} {
    puts "Error: No RNA nucleotides (RA, RU, RC, RG) found in the selection."
    exit
}

# Open output file
set outfile [open "pdb_positions.txt" w]

# Extract positions and write them to a single line
set positions [$sel get {resid x y z}]
set line ""
foreach pos $positions {
    append line "$pos "
}
puts $outfile $line

# Clean up
close $outfile
$sel delete

puts "Nucleotide positions saved to nucleotide_positions.txt"