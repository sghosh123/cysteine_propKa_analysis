# VMD Tcl template: export selected frames to frame_000000.pdb files.
# EDIT ALL PATHS AND SELECTIONS BEFORE RUNNING.
# Run with VMD's Tcl interpreter, e.g. source scripts/export_frames.tcl
# IMPORTANT: the default "protein" selection excludes Cu and other non-protein
# atoms. PROPKA 3.0 reports in the source archive explicitly warn that this
# version does not handle ligands. Do not interpret these as metal-aware pKa
# estimates for a Cu-bound model.

set topology_file "path/to/system.psf"
set trajectory_file "path/to/trajectory.dcd"
set topology_type "psf"
set trajectory_type "dcd"
set atom_selection "protein"
set output_dir "pka_frames"
set first_frame 0
set last_frame -1     ;# -1 means final frame
set stride 1

if {![file exists $topology_file]} { error "Topology not found: $topology_file" }
if {![file exists $trajectory_file]} { error "Trajectory not found: $trajectory_file" }
if {$stride < 1} { error "stride must be >= 1" }

file mkdir $output_dir
set output_dir [file normalize $output_dir]
set mol_id [mol new $topology_file type $topology_type waitfor all]
mol addfile $trajectory_file type $trajectory_type molid $mol_id waitfor all
set nframes [molinfo $mol_id get numframes]
if {$last_frame < 0} { set last_frame [expr {$nframes - 1}] }
if {$first_frame < 0 || $first_frame >= $nframes} { error "first_frame out of range: $first_frame; trajectory has $nframes frames" }
if {$last_frame < $first_frame || $last_frame >= $nframes} { error "last_frame out of range: $last_frame; trajectory has $nframes frames" }

set manifest_path [file join $output_dir "frames_manifest.csv"]
set manifest [open $manifest_path "w"]
puts $manifest "frame_index_zero_based,frame_label_one_based,pdb_file"
for {set im $first_frame} {$im <= $last_frame} {incr im $stride} {
    set sel [atomselect $mol_id $atom_selection frame $im]
    if {[$sel num] == 0} {
        $sel delete
        close $manifest
        error "Selection '$atom_selection' returned zero atoms at frame $im"
    }
    set basename [format "frame_%06d.pdb" $im]
    set pdb_path [file join $output_dir $basename]
    $sel writepdb $pdb_path
    $sel delete
    puts $manifest "$im,[expr {$im + 1}],$basename"
    puts "Exported frame $im -> $pdb_path"
}
close $manifest
puts "Wrote frame manifest: $manifest_path"
