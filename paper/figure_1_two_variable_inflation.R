######################
# Plotting arguments #
######################
filename = "~/Downloads/figure_1_two_variable_inflation.pdf"
width_height_inches = 3.5
width_height_inches_label = 0.1
pointsize = 8

###########
# PANEL A #
###########
# Simplified two-variable model: inflation under the grand null
panel_a = {function(){
	# Simulation arguments
	nsim = 1e7
	nu = 2
	mu = 1
	h = 1
	n = 145
	n.grid = 10^(0:5)
	tau = 9
	rho = 0

	# Functions
	logsum = function(x, y) { mx = pmax(x,y); mx + log(exp(x-mx) + exp(y-mx)) }
	logsum.mat = function(mat) { mx = apply(mat,1,max); mx + log(rowSums(exp(mat-mx))) }
	theorem2 = function(tau, V.len, n=get('n', parent.env(environment())), mu=get('mu', parent.env(environment())), h=get('h', parent.env(environment()))) pchisq(2*log(tau/V.len/mu/sqrt(h/(n+h))), 1, low=FALSE)

	# Simulate
	U1 = rnorm(nsim, 0, 1)
	U2 = rnorm(nsim, 0, 1)

	# Calculate the FPR as a function of n
	calc.log_PO12 = function(n) {
		cn = mu*sqrt(h/(n+h))
		log_PO11 = 2*log(cn) + 0.5*(U1^2 + U2^2)
		log_PO10 = log(cn) + 0.5*(sqrt(1-rho^2)*U1 + rho*U2)^2
		log_PO01 = log(cn) + 0.5*U2^2
		log_PO00 = rep(log(1), nsim)

		logsum.mat(cbind(log_PO10, log_PO01, log_PO11)) - log_PO00
	}
	fpr_PO12 = sapply(n.grid, function(n) mean(calc.log_PO12(n)>=log(tau)))
	e.fpr_PO12 = sapply(n.grid, function(N) theorem2(tau, 2, N))
	environment()
}}()

###########
# PANEL B #
###########
# Simplified two-variable model: inflation under an elementary null
panel_b = {function() {
	# Simulation arguments
	nsim = 1e7
	nu = 2
	mu = 0.1
	h = 1
	n = 145
	beta1 = 0
	beta2.grid = seq(0, 1, len=21)
	rho.grid = seq(0, 1, len=6)
	sigma = 1
	tau = 9

	# Functions
	logsum = function(x, y) { mx = pmax(x,y); mx + log(exp(x-mx) + exp(y-mx)) }
	logsum.mat = function(mat) { mx = apply(mat,1,max); mx + log(rowSums(exp(mat-mx))) }
	theorem2 = function(tau, V.len, n=get('n', parent.env(environment())), mu=get('mu', parent.env(environment())), h=get('h', parent.env(environment()))) pchisq(2*log(tau/V.len/mu/sqrt(h/(n+h))), 1, low=FALSE)

	# Simulate
	Z1 = rnorm(nsim, 0, 1)
	Z2 = rnorm(nsim, 0, 1)

	calc.fpr_PO1 = Vectorize(function(BETA2, RHO) {
		cn = mu*sqrt(h/(n+h))
		U1 = Z1 + beta1*sqrt(n)/sigma
		U2 = Z2 + (RHO*beta1 + sqrt(1-RHO^2)*BETA2)*sqrt(n)/sigma

		log_PO11 = 2*log(cn) + 0.5*(U1^2 + U2^2)
		log_PO10 = log(cn) + 0.5*(sqrt(1-RHO^2)*U1 + RHO*U2)^2
		log_PO01 = log(cn) + 0.5*U2^2
		log_PO00 = rep(log(1), nsim)

		#log_PO1 = (logsum(log_PO10, log_PO11) - logsum(log_PO00, log_PO01))
		mean((logsum(log_PO10, log_PO11) - logsum(log_PO00, log_PO01)) >= log(tau))
	})

	log10_fpr_PO1 = log10(outer(beta2.grid, rho.grid, calc.fpr_PO1))
	e.fpr_PO1 = theorem2(tau, 1)
	environment()
}}()


# Generate the figure
pdf(filename, width=2*width_height_inches, height=width_height_inches, pointsize=pointsize)
par(mar=c(4.5,4.5,1.5,1.5), mfrow=c(1,2))
# Plot panel A
with(panel_a, {
	plot(log10(n.grid), log10(fpr_PO12), xlab=expression(log[10]*" Sample size"), ylab=expression(log[10]*' '*Pr(PO[A[v]:O[v]]*' > '*tau*' ; '*theta)), xaxs="i", yaxs="i", type="o", lwd=2, pch=19)
	curve(log10(theorem2(tau, 2, 10^x)), log10(min(n.grid)), log10(max(n.grid)), add=TRUE, col="green3", lwd=2)
	# Add panel label
	x_user = grconvertX(width_height_inches_label, from="inches", to="user")
	y_user = grconvertY(width_height_inches-width_height_inches_label, from="inches", to="user")
	text(x_user, y_user, labels="A", font=2, cex=1.5, xpd=TRUE)
})
# Plot panel B
with(panel_b, {
	COL = colorRampPalette(c("black","grey80"))(length(rho.grid))
	plot(range(beta2.grid), range(log10_fpr_PO1[is.finite(log10_fpr_PO1)]), type="n", xlab=expression(beta[2]), ylab=expression(log[10]*' '*Pr(PO[A[v]:O[v]]*' > '*tau*' ; '*theta)), xaxs="i",  ylim=c(min(log10_fpr_PO1[is.finite(log10_fpr_PO1)]), -1.8))
	for(i in 1:length(rho.grid)) {
		lines(beta2.grid, log10_fpr_PO1[, i], col=COL[i], lwd=2)
	}
	abline(h = log10(e.fpr_PO1), col="green3", lty=2, lwd=2)
	text(c(0.8, 0.57, 0.45, 0.33, 0.25), -c(2.0, 2.6, 3.0, 3.35, 3.6), sprintf("%.1f", rev(rho.grid)[-1]))
	text(0.15, -5.40, expression(rho=='1.0'))
	# Add panel label
	x_user = grconvertX(width_height_inches+width_height_inches_label, from="inches", to="user")
	y_user = grconvertY(width_height_inches-width_height_inches_label, from="inches", to="user")
	text(x_user, y_user, labels="B", font=2, cex=1.5, xpd=TRUE)
})
dev.off()
