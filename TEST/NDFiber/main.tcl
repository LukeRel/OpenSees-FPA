wipe

puts ""
puts "**************************************************************************"
puts "* Starting NDFiber section shelf with additional prestress strains by LP *"
puts "**************************************************************************"
puts ""
puts "In this model, a 100x10x10cm shelf gets pushed over a 10cm distance."
puts "The model is made of an NDFiber section with prestress effects added"
puts "as additional consistent fiber strains."
puts ""

# some common info
set process_id [getPID]
set is_parallel 0

model basic -ndm 3 -ndf 6

puts "- Model created"
puts ""

# Serie lineare:
timeSeries Linear 1

# Materials:
#nDMaterial ElasticIsotropic 1 30000.0 0.2
#nDMaterial ElasticIsotropic 2 210000.0 0.3
uniaxialMaterial Concrete02 1 -40.0 -0.002 -5.0 -0.05 0.05 3.2 35000.0
uniaxialMaterial Steel02 2 450.0 210000.0 0.0 15.0 0.925 0.15
uniaxialMaterial Steel02 3 1760.0 210000.0 0.0 15.0 0.925 0.15

puts "- Materials added"
puts ""

# Sezioni:
source fibers.tcl

puts "- NDFiber section succesfully created"
puts ""

# Nodi:
node 1 0 0 0
node 2 0 0 1000

# Geometric transformation command
geomTransf Linear 1 1.0 0.0 -0.0

puts "- Nodes created"
puts ""

# Vettore dei tag di sezione, uno per punto di Gauss:
set secTags "5 5 5 5 5 5 5"

# Elemento trave:
#element forceBeamColumn 1 1 2 1 DistReg Lobatto 7 2 50.0 $secTags
element forceBeamColumn 1 1 2 1 Lobatto 5 9

puts "- Elements created"
puts ""

# Analisi:
source pushover_analysis.tcl

wipe
