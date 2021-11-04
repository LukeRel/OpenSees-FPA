wipe

# some common info
set process_id [getPID]
set is_parallel 0

# model creation
model basic -ndm 3 -ndf 3

# timeSeries Linear
timeSeries Linear 1

# materials
nDMaterial ElasticIsotropic 1 30000.0 0.2

# nodes tag x y z
node 1 1 1 0
node 2 1 0 0
node 3 1 1 1
node 4 1 0 1
node 5 0 1 0
node 6 0 0 0
node 7 0 1 1
node 8 0 0 1

# brick_elements SSPbrick
element SSPbrick 19 6 5 7 8 2 1 3 4 1

# source analysis_steps
source analysis_steps.tcl

wipe
