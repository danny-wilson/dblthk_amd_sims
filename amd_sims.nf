process simulate_dry_run {
	input:
		val taskid
		tuple val(params), val(n)
	output:
		tuple val(params), val(n)
		path("taskid"), emit: taskid
		path("sim.data.RDS"), emit: sim.data
		path("sim.analyses.RDS"), emit: sim.analyses
		path("sim.perf.RDS"), emit: sim.perf
shell:
'''
	echo !{taskid} !{params} !{n}
'''
}

process simulate {
beforeScript "module add R/4.2.1-foss-2022a"
input:
	val taskid
	tuple val(params), val(n)
output:
	tuple val(params), val(n)
	path("taskid"), emit: taskid
	path("sim.data.RDS"), emit: sim.data
	path("sim.analyses.RDS"), emit: sim.analyses
	path("sim.perf.RDS"), emit: sim.perf
shell:
'''
	#!/usr/bin/env Rscript

	# Get working directory
	wd = getwd()
	setwd("/well/bag/wilson/GitHub/dblthk_amd_sims")
	# Load required source code: this in turn sources summary_mvMR_BF.R
	source("biomarker-sim-functions.R")
	# Return to working directory
	setwd(wd)

	# Read arguments
	taskid = as.integer("!{taskid}")
	params = as.numeric("!{params}")
	n = as.integer("!{n}")

	# Check arguments
	stopifnot(!is.na(taskid))
	stopifnot(!any(is.na(params)))
	stopifnot(length(params)==1)
	params = rep(params, 15)
	stopifnot(n>15)

	# Load example AMD data
	full.data = load.data()

	# Perform univariable association tests
	results.univariable.MR = univariable.MR(full.data)

	# Filter example data to the 15 most significant risk factors based on univariable associations
	data = filter.data(full.data, as.integer(15), results.univariable.MR)

	# Set seed to fix the independent variables (x)
	set.seed(0)
	# Simulate the independent variables based on the example AMBdata
	simulate = gen.simulate(data, n)
	
	# Set seed again to unfix the dependent variables (y)
	set.seed(taskid)
	# Simulate the dependent variables
	sim.data = simulate(n, params)
	
	# Perform the standard set of analyses
	sim.analyses = do.analyses(sim.data, nu=full.data@m, mr.bma.nsim=0)

	# Compute performance metrics for the analyses
	sim.perf = calc.performance(sim.analyses, params, freqt.alpha=0.05, bayes.tau=19)

	# Save files
	write(taskid, file="taskid")
	saveRDS(sim.data, file="sim.data.RDS")
	saveRDS(sim.analyses, file="sim.analyses.RDS")
	saveRDS(sim.perf, file="sim.perf.RDS")
	
	cat("Doublethink AMD simulation completed successfully\n")
'''
}

// Default values
params.ntasks = 1000
params.params = 0.0
params.n = 149

// Define the channels
ch_taskid = Channel.of(1..params.ntasks)
ch_params = Channel.of(params.params)
ch_n = Channel.of(params.n)
ch_args = ch_n
				   .combine(ch_params)

// Go
workflow {
	// Perform the simulations
	simulate(ch_taskid, ch_args)

}
