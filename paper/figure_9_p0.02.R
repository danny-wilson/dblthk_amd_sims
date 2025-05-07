#############
# Arguments #
#############
filename = "~/Downloads/figure_9_p_0.02.pdf"
width_height_inches = 3.5
pointsize = 8
L = 10
nrep = 1e7
df = 1

#############
# Functions #
#############
ptail.1 = function(x) pchisq(2*log(x), df, low=FALSE)
ptail.approx = function(x, L) L*ptail.1(L*x)
ptail.2 = Vectorize(function(y) 1-integrate(function(x) (1-ptail.1(2*y-exp(0.5*x))) * dchisq(x, df), 0, 2*log(2*y))$value)
f = function(y) ptail.1(y) - ptail.2(y)

############
# Simulate #
############
L.grid = 1:L
system.time((Xbar = sapply(L.grid,function(L.GRID) {
    X = matrix(exp(0.5*rchisq(nrep*L.GRID, df)),nrow=nrep)
    sort(rowMeans(X))
})))
#   user  system elapsed 
# 80.378   3.938  84.320 
Xbar.range = range(Xbar)
Xbar.length.range = c(1,nrep)

########
# Plot #
########
quartz(file=filename, width=width_height_inches, height=width_height_inches, pointsize=pointsize, type="pdf")
closeup = FALSE
par(mar=c(4.5, 4.5, 1.5, 1.5))
COL = rainbow(length(L.grid))
XXLIM = c(1, 200)
XLIM = c(1, 100)
YYLIM = c(ptail.approx(XXLIM[2], max(L.grid)), 1)
if(closeup) { XXLIM = c(11, 13); YYLIM = c(ptail.1(XXLIM[2]), ptail.1(XXLIM[1])) }
plot(XXLIM, YYLIM, log="xy", type="n", xlab="x", ylab=expression(Pr(bar(X)[k] > x)))
for(i in length(L.grid):1) {
	curve(1-ecdf(Xbar[,i])(x), XLIM[1], XLIM[2], add=TRUE, col=COL[i])
	if(!closeup) curve(ptail.approx(x, L.grid[i]), XLIM[2], XXLIM[2], add=TRUE, col=COL[i], lty=3)
}

# The crossover value of Rbar and p
(Rbar.root = uniroot(f, c(1.1, 100), tol=1e-12)$root)	# 11.92362
(ptail.1.root = ptail.1(Rbar.root))			# 0.0259846
abline(h=ptail.1.root, col=8)
abline(v=Rbar.root, col=8)

legend("topright", paste("k =", L.grid), col=COL, lwd=1, bty="n")
text(Rbar.root, YYLIM[1], sprintf("%.3g", Rbar.root), adj=c(-0.2, 1), col=8)
text(XXLIM[1], ptail.1.root, sprintf("%.3g", ptail.1.root), adj=c(0, -1), col=8)
if(names(dev.cur())!="quartz") dev.off()
