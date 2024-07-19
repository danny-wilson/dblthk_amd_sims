process simulate_dry_run {
input:
	val taskid
	path parameters_filename
output:
	path("${taskid}.sim.data.RDS"), 	emit: sim_data
	path("${taskid}.sim.anal.RDS"),		emit: sim_anal
	path("${taskid}.sim.perf.RDS"), 	emit: sim_perf
shell:
'''
	echo !{taskid} !{params} !{n}
'''
}

// Task batching pattern https://nextflow-io.github.io/patterns/task-batching/
process simulate {
input:
	val taskids
	path parameters_filename
output:
	path("${taskids}.sim.data.RDS"), 		emit: sim_data
	path("${taskids}.sim.anal.RDS"),		emit: sim_anal
	path("${taskids}.sim.perf.RDS"), 		emit: sim_perf
shell:
'''
	#!/usr/bin/env Rscript

	# Get working directory
	wd = getwd()
	cat("repoDir:                !{params.repoDir}\n")
	setwd("!{params.repoDir}")
	# Load required source code: this in turn sources summary_mvMR_BF.R
	source("biomarker-sim-functions.R")
	# Return to working directory
	setwd(wd)

	# Read arguments
	taskids = as.integer(unlist(strsplit("!{taskids}", split=" ")))
	cat("simulate() read arguments:\n")
	cat("taskids:               ", taskids, "\n")
	stopifnot(length(taskids)>0)
	stopifnot(!any(is.na(taskids)))
	stopifnot(all(taskids>0))
	stopifnot(length(unique(taskids))==length(taskids))
	infile = as.character("!{parameters_filename}")
	filenames.sim_data = paste0(taskids, ".sim.data.RDS")
	filenames.sim_anal = paste0(taskids, ".sim.anal.RDS")
	filenames.sim_perf = paste0(taskids, ".sim.perf.RDS")
	cat("infile:                ", infile, "\n")
	cat("filenames.sim_data:    ", filenames.sim_data, "\n")
	cat("filenames.sim_anal:    ", filenames.sim_anal, "\n")
	cat("filenames.sim_perf:    ", filenames.sim_perf, "\n")

	# Implied arguments
	nsim = as.integer("!{params.ntasks}")
	cat("nsim:                  ", nsim, "\n")
	stopifnot(!is.na(nsim))
	stopifnot(nsim>0)
	stopifnot(all(nsim>=taskids))
	n = as.numeric("!{params.n}")
	cat("n:                     ", n, "\n")
	stopifnot(!is.na(n))
	stopifnot(n>15)
	mr_bma_nsim = as.integer("!{params.mr_bma_nsim}")
	cat("mr_bma_nsim:           ", mr_bma_nsim, "\n")
	stopifnot(!is.na(mr_bma_nsim))
	stopifnot(mr_bma_nsim>=0)
	alpha = as.double("!{params.alpha}")
	stopifnot(!is.na(alpha))
	stopifnot(alpha>0)
	cat("alpha:                 ", alpha, "\n")
	tau = as.double("!{params.tau}")
	stopifnot(!is.na(tau))
	stopifnot(tau>0)
	cat("tau:                   ", tau, "\n")
	simulate_independence = as.logical(toupper("!{params.simulate_independence}"))
	stopifnot(!is.na(simulate_independence))
	cat("simulate_independence: ", simulate_independence, "\n")

	# Read parameters
	stopifnot(file.exists(infile))
	all_params = as.matrix(read.delim(infile, header=FALSE, sep=" ", quote=""))
	stopifnot(nrow(all_params)==nsim)
	stopifnot(ncol(all_params)==15)
	stopifnot(!any(is.na(all_params)))

	# Load example AMD data
	setwd("!{params.repoDir}")
	full.data = load.data()
	setwd(wd)

	# Perform univariable association tests
	results.univariable.MR = univariable.MR(full.data)

	# Filter example data to the 15 most significant risk factors based on univariable associations
	data = filter.data(full.data, as.integer(15), results.univariable.MR)

	# Set seed to fix the independent variables (x)
	set.seed(0)
	# Simulate the independent variables based on the example AMBdata
	simulate = gen.simulate(data, n, seed=NA, independence=simulate_independence)
	
	for(subtask in 1:length(taskids)) {
		taskid = taskids[subtask]

		params = all_params[taskid,]
		cat("Beginning taskid:    ", taskid, "\n\n")
		cat("params:              ", params, "\n\n")

		# Set seed again to unfix the dependent variables (y)
		set.seed(taskid)
		# Simulate the dependent variables
		sim.data = simulate(n, params)
		new.data = simulate(n, params)
		
		# Perform the standard set of analyses
		sim.anal = do.analyses(sim.data, params, nu=sim.data@m, mr.bma.nsim=mr_bma_nsim)

		# Compute performance metrics for the analyses
		sim.perf = calc.performance(sim.anal, params, freqt.alpha=alpha, bayes.tau=tau, newdata=new.data)

		# Save files
		saveRDS(sim.data, file=filenames.sim_data[subtask])
		saveRDS(sim.anal, file=filenames.sim_anal[subtask])
		saveRDS(sim.perf, file=filenames.sim_perf[subtask])
	}
	
	cat("simulate: Doublethink AMD simulation completed successfully\n")
'''
}

process combine_performance {
publishDir "${params.publishDir}", mode: "copy"
input:
	path(infiles)
output:
	path("combined.performance.RDS"), emit: combined_perf
shell:
'''
	#!/usr/bin/env Rscript
	
	# Get working directory
	wd = getwd()
	cat("repoDir:               !{params.repoDir}\n")
	setwd("!{params.repoDir}")
	# Load required source code: this in turn sources summary_mvMR_BF.R
	source("biomarker-sim-functions.R")
	# Return to working directory
	setwd(wd)

	# Read arguments
	infiles.all = "!{infiles}"
	outfile = "combined.performance.RDS"
	
	# Check arguments
	infiles = unlist(strsplit(infiles.all, ' '))
	nfiles = length(infiles)
	stopifnot(file.exists(infiles))
	stopifnot(!file.exists(outfile))
	cat("combine_performance: combining", nfiles, "files\n")
	
	# Print arguments
	cat("combine_performance() read arguments:\n")
	cat("infiles:              ", infiles, "\n")
	cat("outfile:              ", outfile, "\n")
	cat("\n")

	results = list()
	for(i in 1:nfiles) {
		# Load simulation performance metrics
		perf = readRDS(infiles[i])
		if(i==1) {
			nanalyses = length(perf)
		} else {
			stopifnot(nanalyses==length(perf))
		}
	
		# Combine analysis results
		for(j in 1:nanalyses) {
			# Combine the performance metrics
			stopifnot(is(perf[[j]], "bmsim_performance"))
			analysis.name = names(perf)[j]
			if(i==1) {
				results[[analysis.name]] = combine.performance.iteratively(NULL, perf[[analysis.name]], i, nfiles)
			} else {
				results[[analysis.name]] = combine.performance.iteratively(results[[analysis.name]], perf[[analysis.name]], i, nfiles)
			}
		}

		# Free memory
		rm("perf")
		if(i %% 100==0) gc()
	}

	# Save results
	saveRDS(results, file="combined.performance.RDS")

	cat("combine_performance: Completed successfully\n")
'''
}

process simulate_parameters {
publishDir "${params.publishDir}", mode: "copy"
input:
output:
	path("${params.model}.parameters.txt"), emit: model_params
shell:
'''
	#!/usr/bin/env Rscript

	# Get working directory
	wd = getwd()
	cat("repoDir:              !{params.repoDir}\n")
	setwd("!{params.repoDir}")
	# Load required source code: this in turn sources summary_mvMR_BF.R
	source("biomarker-sim-functions.R")
	
	# Implied arguments
	nsim = as.integer("!{params.ntasks}")
	stopifnot(!is.na(nsim))
	stopifnot(nsim>0)
	model = as.character("!{params.model}")
	filename.parameters = "!{params.model}.parameters.txt"
	tau = as.numeric("!{params.simulate_parameters_tau}")
	stopifnot(!is.na(tau))
	stopifnot(tau>0)
	h = as.numeric("!{params.simulate_parameters_h}")
	stopifnot(!is.na(h))
	stopifnot(h>0)
	mu = as.numeric("!{params.simulate_parameters_mu}")
	stopifnot(!is.na(mu))
	stopifnot(mu>0)
	n = as.numeric("!{params.n}")
	stopifnot(!is.na(n))
	stopifnot(n>0)
	seed = as.integer("!{params.simulate_parameters_seed}")
	stopifnot(!is.na(seed))

	# Print arguments
	cat("simulate_parameters() read arguments:\n")
	cat("nsim:                ", nsim, "\n")
	cat("model:               ", model, "\n")
	cat("filename.parameters: ", filename.parameters, "\n")
	cat("tau:                 ", tau, "\n")
	cat("h:                   ", h, "\n")
	cat("mu:                  ", mu, "\n")
	cat("n:                   ", n, "\n")
	cat("seed:                ", seed, "\n")
	cat("\n")

	# Load example AMD data
	full.data = load.data()

	# Return to working directory
	setwd(wd)

	# Perform univariable association tests
	results.univariable.MR = univariable.MR(full.data)

	# Filter example data to the 15 most significant risk factors based on univariable associations
	data = filter.data(full.data, as.integer(15), results.univariable.MR)

	# Simulate parameters and write to file
	simulate.parameters(nsim=nsim, m=data@m, model=model, filename=filename.parameters, overwrite=FALSE, data=data, tau=tau, h=h, mu=mu, n=n, seed=seed)
'''
}

// Default values
params.ntasks = 1000
params.n = 145
params.model = 'grand_null'
params.publishDir = '.'
params.simulate_parameters_tau = 9
params.simulate_parameters_h = 1
params.simulate_parameters_mu = 0.1
params.simulate_parameters_seed = 0
params.mr_bma_nsim = 1000
params.alpha = 0.01
params.tau = 9
params.simulate_independence = false
params.repoDir = '/well/bag/wilson/GitHub/dblthk_amd_sims'
params.combine_performance_mem = 14.0
params.task_batch_size = 10

// Print arguments
// Default arguments can be overriden by specifying them in nextflow.config
// Unspecified arguments with no default will throw an error here.
println 'Arguments'
println 'ntasks:			       ' + params.ntasks
println 'n:					       ' + params.n
println 'model:				       ' + params.model
println 'publishDir:		       ' + params.publishDir
println 'simulate_parameters_tau:  ' + params.simulate_parameters_tau
println 'simulate_parameters_h:    ' + params.simulate_parameters_h
println 'simulate_parameters_mu:   ' + params.simulate_parameters_mu
println 'simulate_parameters_seed: ' + params.simulate_parameters_seed
println 'mr_bma_nsim:              ' + params.mr_bma_nsim
println 'alpha:                    ' + params.alpha
println 'tau:                      ' + params.tau
println 'simulate_independence:    ' + params.simulate_independence
println 'repoDir:                  ' + params.repoDir
println 'combine_performance_mem:  ' + params.combine_performance_mem
println 'task_batch_size        :  ' + params.task_batch_size

// Define the channels
ch_taskid = Channel.of(1..params.ntasks) | buffer(size: params.task_batch_size, remainder: true)

// Go
workflow {
	// Simulate the parameters
	simulate_parameters()
	
	// Perform the simulations
	simulate(ch_taskid, simulate_parameters.out.model_params)

	// Combine performance metrics
	combine_performance(simulate.out.sim_perf.collect())
}
