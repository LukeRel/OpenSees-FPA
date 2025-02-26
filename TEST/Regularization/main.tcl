wipe

# some common info
set process_id [getPID]
set is_parallel 0

model basic -ndm 3 -ndf 6

# Serie lineare:
timeSeries Linear 1

# Materials:
uniaxialMaterial Concrete02 1 -40.0 -0.002 -5.0 -0.01 0.05 3.2 35000.0
uniaxialMaterial Steel02 2 450.0 210000.0 0.0 15.0 0.925 0.15
# nDMaterial ElasticIsotropic 3 210000.0 0.3
# nDMaterial ElasticIsotropic 4 210000.0 0.3
# nDMaterial PlasticDamageConcrete3d 5 35000.0 0.2 -40.0 3.0

# Sezioni:
source fibers.tcl

# Nodi:
node 1 0 0 0
node 2 0 0 1000

# Geometric transformation command
geomTransf Linear 1 1.0 0.0 -0.0

# Vettore dei tag di sezione, uno per punto di Gauss:
set secTags "5 5 5 5 5 5 5"

# Elemento trave:
element forceBeamColumn 1 1 2 1 DistReg Lobatto 7 2 50.0 $secTags

# Analisi:
source pushover_analysis.tcl

wipe
