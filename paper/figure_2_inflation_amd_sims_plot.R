#############
# Arguments #
#############
# Working directory
workdir = "~/Downloads"
# Plotting arguments
width_height_inches = 3.5
width_height_inches_label = 0.1
pointsize = 8
# Hyper-parameters
mu = 0.1; h = 1; hsim = h
tau = 9
nu = 15 # Needs to match data@m
# Input filenames and parameters
nsamps = list('A' = 145, 'B' = 14500)
input_filenames = lapply(nsamps, function(nsamp) {
	# Number of simulations
	nsim = ceiling(50000 * sqrt(nsamp/145))
	paste0(workdir, "/nsamp", nsamp, "_nsim", nsim, "_mu", mu, "_h", h, "_", format(Sys.time(), "%Y-%M-%d"))
})
#input_filenames = list('A' = paste0("~/doublethink/more-efficient-sims/nsamp145_nsim50000_mu0.1_h1_rightprior_new20240716.RDS"),
#                       'B' = paste0("~/doublethink/more-efficient-sims/nsamp14500_nsim5e+05_mu0.1_h1_rightprior_new20240716.RDS"))
# Output filename
figure_filename = paste0(workdir, "/figure_2_inflation_amd_sims.pdf")

#############
# Functions #
#############

do_panel = function(panel = 'A') {
	#############
	# Arguments #
	#############
	# Set the sample size independently of the example data
	nsamp = nsamps[[panel]]
	xi = h/(nsamp+h)
	# Input filenames
	input_filename = input_filenames[[panel]]

	#######################
	# Read the simulation #
	#######################
	stopifnot(file.exists(input_filename))
	sim = readRDS(input_filename)

	#######################
	# Do the computations #
	#######################
	fwel.tol = as.numeric(gsub('logPO.truenull.num.', '', colnames(sim)[grep('logPO.truenull.num', colnames(sim))]))
	ptruenull.fp = colMeans(sim[,grep('logPO.truenull.num', colnames(sim))] - sim[,grep('logPO.truenull.den', colnames(sim))] >= log(tau), na.rm=TRUE)
	pasymnull.fp = colMeans(sim[,grep('logPO.asymnull.num', colnames(sim))] - sim[,grep('logPO.asymnull.den', colnames(sim))] >= log(tau), na.rm=TRUE)
	ptheory = pchisq(2*log(tau/nu/mu/sqrt(xi)), 1, low=FALSE)
	nsim.truenull.fp = colSums(!is.na(sim[,grep('logPO.truenull.num', colnames(sim))] - sim[,grep('logPO.truenull.den', colnames(sim))] >= log(tau)))
	nsim.asymnull.fp = colSums(!is.na(sim[,grep('logPO.asymnull.num', colnames(sim))] - sim[,grep('logPO.asymnull.den', colnames(sim))] >= log(tau)))
	environment()
}

#################
# Load the data #
#################
panel_a = do_panel('A')
panel_b = do_panel('B')

########
# Plot #
########

pdf(figure_filename, width=2*width_height_inches, height=width_height_inches, pointsize=pointsize)
par(mar=c(4.5, 4.5, 1.5, 1.5), mfrow=c(1,2))
options(scipen = 999)
# Panel A
with(panel_a, {
	XLIM = range(sqrt(fwel.tol))
	YLIM = c(0, max(c(qbinom(0.975, nsim.asymnull.fp, pasymnull.fp)/nsim.asymnull.fp, qbinom(0.975, nsim.truenull.fp, ptruenull.fp)/nsim.truenull.fp)))
	plot(XLIM, YLIM, xlim=XLIM, ylim=YLIM, type="n", xlab=expression(rho), ylab=expression(BFWER[rho]), lwd=2)
	abline(h=ptheory, col="green3", lwd=2)
	points(sqrt(fwel.tol), pasymnull.fp, col="grey70", lwd=2, pch=19)
	arrows(sqrt(fwel.tol), qbinom(0.025, nsim.asymnull.fp, pasymnull.fp)/nsim.asymnull.fp, sqrt(fwel.tol), qbinom(0.975, nsim.asymnull.fp, pasymnull.fp)/nsim.asymnull.fp, len=0.05, angle=90, code=3, col="grey70", lwd=2)
	points(sqrt(fwel.tol), ptruenull.fp, col="black", lwd=2, pch=19)
	arrows(sqrt(fwel.tol), qbinom(0.025, nsim.truenull.fp, ptruenull.fp)/nsim.truenull.fp, sqrt(fwel.tol), qbinom(0.975, nsim.truenull.fp, ptruenull.fp)/nsim.truenull.fp, len=0, col="black", lwd=2)
	# Add panel label
	x_user = grconvertX(width_height_inches_label, from="inches", to="user")
	y_user = grconvertY(width_height_inches-width_height_inches_label, from="inches", to="user")
	text(x_user, y_user, labels="A", font=2, cex=1.5, xpd=TRUE)
})
# Panel B
with(panel_b, {
	XLIM = range(sqrt(fwel.tol))
	YLIM = c(0, max(c(qbinom(0.975, nsim.asymnull.fp, pasymnull.fp)/nsim.asymnull.fp, qbinom(0.975, nsim.truenull.fp, ptruenull.fp)/nsim.truenull.fp)))
	plot(XLIM, YLIM, xlim=XLIM, ylim=YLIM, type="n", xlab=expression(rho), ylab=expression(BFWER[rho]), lwd=2)
	abline(h=ptheory, col="green3", lwd=2)
	points(sqrt(fwel.tol), pasymnull.fp, col="grey70", lwd=2, pch=19)
	arrows(sqrt(fwel.tol), qbinom(0.025, nsim.asymnull.fp, pasymnull.fp)/nsim.asymnull.fp, sqrt(fwel.tol), qbinom(0.975, nsim.asymnull.fp, pasymnull.fp)/nsim.asymnull.fp, len=0.05, angle=90, code=3, col="grey70", lwd=2)
	points(sqrt(fwel.tol), ptruenull.fp, col="black", lwd=2, pch=19)
	arrows(sqrt(fwel.tol), qbinom(0.025, nsim.truenull.fp, ptruenull.fp)/nsim.truenull.fp, sqrt(fwel.tol), qbinom(0.975, nsim.truenull.fp, ptruenull.fp)/nsim.truenull.fp, len=0, col="black", lwd=2)
	# Add panel label
	x_user = grconvertX(width_height_inches+width_height_inches_label, from="inches", to="user")
	y_user = grconvertY(width_height_inches-width_height_inches_label, from="inches", to="user")
	text(x_user, y_user, labels="B", font=2, cex=1.5, xpd=TRUE)
})
dev.off()
