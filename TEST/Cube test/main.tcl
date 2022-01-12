wipe

# some common info
set process_id [getPID]
set is_parallel 0

# model creation
model basic -ndm 3 -ndf 3

# timeSeries Linear
timeSeries Linear 1

# materials
nDMaterial PlasticDamage2P 1 37000 0.2 6.6 60.0 25900.0 37.0 9.72e-5 3.0e-4 0.8 5.84e-5 5.0e-4 0.1 1.0

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
