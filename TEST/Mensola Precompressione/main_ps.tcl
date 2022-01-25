wipe

# some common info
set process_id [getPID]
set is_parallel 0

model basic -ndm 3 -ndf 6
# source definitions
source definitions.tcl
# source materials
source materials_1.tcl
# source sections
source sections_1.tcl
# source node
source nodes.tcl
# source element
source elements.tcl
# source analysis_steps
source analysis_steps.tcl

wipe
