set graph_style [lindex $argv 2]
# Set graphical representation
puts "Setting graphical representations..."
if { $graph_style == "VDW" } {
    mol modselect 0 0 all
    mol modstyle 0 0 VDW
} elseif { $graph_style == "Licorice" } {
    mol modselect 0 0 all
    mol modstyle 0 0 Licorice
} else {
    puts "Unknown graphical style: $graph_style. Defaulting to VDW."
    mol modselect 0 0 all
    mol modstyle 0 0 VDW
}

# Apply graphical settings
mol color Name
mol material Transparent
puts "Graphical settings applied."
