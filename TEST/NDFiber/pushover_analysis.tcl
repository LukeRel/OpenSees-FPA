
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

recorder mpco "shelf_noprestress.mpco" \
-N "displacement" "rotation" "reactionForce" "reactionMoment" \
-E "force" "deformation" "section.force" "section.deformation" "section.fiber.stress" "section.fiber.strain"

# Constraints.sp fix
	fix 1 1 1 1 1 1 1

# Patterns.addPattern loadPattern
pattern Plain 3 1 {

# Loads.sp sp(prescribedNodalValue)
 # sp node
	sp 2 1 50.0
}

# analyses command
domainChange
constraints Transformation
numberer RCM
system UmfPack
test NormDispIncr 0.0001 1000  
algorithm NewtonLineSearch
integrator LoadControl 0.0
analysis Static
# ======================================================================================
# ADAPTIVE LOAD CONTROL ANALYSIS
# ======================================================================================

# ======================================================================================
# USER INPUT DATA 
# ======================================================================================

# duration and initial time step
set total_time 1.0
set initial_num_incr 1000

# parameters for adaptive time step
set max_factor 1.0
set min_factor 1e-05
set max_factor_increment 1.5
set min_factor_increment 1e-05
set max_iter 1000
set desired_iter 500

set increment_counter 0
set factor 1.0
set old_factor $factor
set time 0.0
set initial_time_increment [expr $total_time / $initial_num_incr]
set time_tolerance [expr abs($initial_time_increment) * 1.0e-8]

while 1 {
	
	incr increment_counter
	if {[expr abs($time)] >= [expr abs($total_time)]} {
		if {$process_id == 0} {
			puts "Target time has been reached. Current time = $time"
			puts "SUCCESS."
		}
		break
	}
	
	set time_increment [expr $initial_time_increment * $factor]
	if {[expr abs($time + $time_increment)] > [expr abs($total_time) - $time_tolerance]} {
		set time_increment [expr $total_time - $time]
	}
	if {$process_id == 0} {
		puts "Increment: $increment_counter. time_increment = $time_increment. Current time = $time"
	}
	
	integrator LoadControl $time_increment 
	set ok [analyze 1]
	#barrier
	
	if {$ok == 0} {
		set num_iter [testIter]
		set factor_increment [expr min($max_factor_increment, [expr double($desired_iter) / double($num_iter)])]
		set factor [expr $factor * $factor_increment]
		if {$factor > $max_factor} {
			set factor $max_factor
		}
		if {$process_id == 0} {
			if {$factor > $old_factor} {
				puts "Increasing increment factor due to faster convergence. Factor = $factor"
			}
		}
		set old_factor $factor
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
		set num_iter $max_iter
		set factor_increment [expr max($min_factor_increment, [expr double($desired_iter) / double($num_iter)])]
		set factor [expr $factor * $factor_increment]
		if {$process_id == 0} {
			puts "Reducing increment factor due to non convergece. Factor = $factor"
		}
		if {$factor < $min_factor} {
			if {$process_id == 0} {
				puts "ERROR: current factor is less then the minimum allowed ($factor < $min_factor)"
				puts "Giving up"
			}
			error "ERROR: the analysis did not converge"
		}
	}
	
}

wipeAnalysis

# Done!
puts "ANALYSIS SUCCESSFULLY FINISHED"
