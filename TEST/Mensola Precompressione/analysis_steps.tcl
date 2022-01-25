
# a list of all monitor and custom function actors to be called by the MonitorFunction
set all_custom_functions {}
set all_monitor_actors {}

# the main custom function caller that will call all actors in $all_monitor_actors and in $all_custom_functions list
proc CustomFunctionCaller {step_id dt T n_iter norm perc process_id is_parallel} {
	global all_monitor_actors
	global all_custom_functions
	# Call monitors: we pass the parameters needed
	foreach p $all_monitor_actors {
		$p $step_id $dt $T $n_iter $norm $perc $process_id $is_parallel
	}
	# Call all other custom functions
	foreach p $all_custom_functions {
		$p
	}
}

recorder mpco "Prova.mpco" \
-N "displacement" "rotation" "reactionForce" "reactionMoment" \
-E "force" "deformation" "section.force" "section.deformation" "material.stress" "material.strain" "section.fiber.stress" "section.fiber.strain"

# Constraints.sp fix
	fix 1 1 1 1 1 1 1

# Patterns.addPattern loadPattern
pattern Plain 3 1 {

# Loads.Force NodeForce
	load 2 0.1 0.0 0.0 0.0 0.0 0.0
}

# analyses command
domainChange
constraints Penalty 1e+16 1e+16
numberer RCM
system UmfPack
test NormDispIncr 0.0001 100  
algorithm Newton
integrator LoadControl 0.0
analysis Static
# ======================================================================================
# NON-ADAPTIVE LOAD CONTROL ANALYSIS
# ======================================================================================

# ======================================================================================
# USER INPUT DATA 
# ======================================================================================

# duration and initial time step
set total_time 1.0
set initial_num_incr 10

set time 0.0
set time_increment [expr $total_time / $initial_num_incr]
integrator LoadControl $time_increment 
for {set increment_counter 1} {$increment_counter <= $initial_num_incr} {incr increment_counter} {
	if {$process_id == 0} {
		puts "Increment: $increment_counter. time_increment = $time_increment. Current time = $time"
	}
	
	set ok [analyze 1 ]
	#barrier
	
	if {$ok == 0} {
		set num_iter [testIter]
		set time [expr $time + $time_increment]
		# print statistics
		set norms [testNorms]
		if {$num_iter > 0} {set last_norm [lindex $norms [expr $num_iter-1]]} else {set last_norm 0.0}
		if {$process_id == 0} {
			puts "Increment: $increment_counter - Iterations: $num_iter - Norm: $last_norm ( [expr $time/$total_time*100.0] % )"
		}
		
		# Call Custom Functions
		set perc [expr $time/$total_time]
		CustomFunctionCaller $increment_counter $time_increment $time $num_iter $last_norm $perc $process_id $is_parallel

	} else {
		error "ERROR: the analysis did not converge"
	}
}

if {$process_id == 0} {
	puts "Target time has been reached. Current time = $time"
	puts "SUCCESS."
}
wipeAnalysis

# Done!
puts "ANALYSIS SUCCESSFULLY FINISHED"
