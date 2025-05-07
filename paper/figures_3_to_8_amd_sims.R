# Based on 'explore combined performance-chromethroat.R'
setwd("biomarker MR/amd-sim")
combined_performance_file = "~/Downloads/n145_fwer_rho/combined.performance.RDS"
outdir = "~/Downloads"

#############
# Functions #
#############
require(grDevices)
se = function(alpha) c(qbinom(.025, nsim, alpha)/nsim, alpha, qbinom(.975, nsim, alpha)/nsim)
se.approx = function(alpha) alpha - 2*sqrt(alpha/nsim)
len02na = function(x) ifelse(length(x)==0, NA, x)
rk = function(x) ifelse(is.na(x), NA, sapply(x, function(X) sum(setdiff(unique(x), NA)<=X, na.rm=TRUE)))

################
# Load results #
################

res = readRDS(combined_performance_file)
# From the relevant nextflow.config, manually set the thresholds and # simulations
tau = 9; alpha = 0.01; nsim = 10000
fdr = 1/(1+tau)

# Set the reference analysis
dba = "doublethink bma mu = 0.1; h = 1"

# Remove 'bma preciser', 'bma 1df preciser', 'mr-bma bma permute'
for(i in names(res)[grep("bma preciser", names(res))]) res[[i]] <- NULL
for(i in names(res)[grep("bma 1df preciser", names(res))]) res[[i]] <- NULL
res[["mr-bma bma permute"]] <- NULL

# Names of remaining analyses
ana.pos = c(
"ridge regression"=2,
"grand null"=1,
"backward elimination"=3,
"mr-bma bma"=4,
"mr-bma model selection"=4,
"elastic net model selection"=5,
"elastic net"=5,
"lasso model selection"=6,
"lasso"=6,
"doublethink model selection mu = 0.2; h = 4"=7,
"doublethink model selection mu = 0.1; h = 4"=8,
"doublethink model selection mu = 0.05; h = 4"=9,
"doublethink model selection mu = 0.2; h = 1"=10,
"doublethink model selection mu = 0.1; h = 1"=11,
"doublethink model selection mu = 0.05; h = 1"=12,
"doublethink model selection mu = 0.2; h = 0.25"=13,
"doublethink model selection mu = 0.1; h = 0.25"=14,
"doublethink model selection mu = 0.05; h = 0.25"=15,
"grand alternative"=16,
"doublethink bma mu = 0.2; h = 4"=7,
"doublethink bma mu = 0.1; h = 4"=8,
"doublethink bma mu = 0.05; h = 4"=9,
"doublethink bma mu = 0.2; h = 1"=10,
"doublethink bma mu = 0.1; h = 1"=11,
"doublethink bma mu = 0.05; h = 1"=12,
"doublethink bma mu = 0.2; h = 0.25"=13,
"doublethink bma mu = 0.1; h = 0.25"=14,
"doublethink bma mu = 0.05; h = 0.25"=15,
"oracle"=17
)
ana = names(ana.pos)
ana.lab = c(
"Ridge regression",
"Grand null",
"Backward elimination",
"MR-BMA",
"MR-BMA",
"Elastic net",
"Elastic net",
"LASSO",
"LASSO",
"Doublethink mu = 0.2; h = 4",
"Doublethink mu = 0.1; h = 4",
"Doublethink mu = 0.05; h = 4",
"Doublethink mu = 0.2; h = 1",
"Doublethink mu = 0.1; h = 1",
"Doublethink mu = 0.05; h = 1",
"Doublethink mu = 0.2; h = 0.25",
"Doublethink mu = 0.1; h = 0.25",
"Doublethink mu = 0.05; h = 0.25",
"Grand alternative",
"Doublethink mu = 0.2; h = 4",
"Doublethink mu = 0.1; h = 4",
"Doublethink mu = 0.05; h = 4",
"Doublethink mu = 0.2; h = 1",
"Doublethink mu = 0.1; h = 1",
"Doublethink mu = 0.05; h = 1",
"Doublethink mu = 0.2; h = 0.25",
"Doublethink mu = 0.1; h = 0.25",
"Doublethink mu = 0.05; h = 0.25",
"Oracle"
)
ana.lab = gsub("mu", "μ", ana.lab)
# For each, sub-analyses
subana = c("", ".Bonf", ".Simes", ".Hommel", ".HMP")
# Combination
ana.sub = outer(ana, subana, paste0); rownames(ana.sub) <- ana; colnames(ana.sub) <- subana
# Pch for sub-analyses
sub.pch = c(1, rep(124, length(subana)-1)); sub.cex = c(1, seq(1, by=-0.2, len=length(subana)-1)); sub.lwd = 1/sub.cex

# Col for main analyses
col = rep("grey", length(ana)); names(col) <- ana
for(i in 1:length(col)) {
	if(grepl("doublethink bma mu", names(col)[i])) col[i] = "blue2"
	if(names(col)[i]=="doublethink bma mu = 0.1; h = 1") col[i] = "blue2"
	if(grepl("doublethink bma 1df", names(col)[i])) col[i] = "skyblue2"
	if(grepl("doublethink model selection", names(col)[i])) col[i] = "skyblue"
	if(grepl("mr-bma bma", names(col)[i])) col[i] = "red3"
	if(grepl("mr-bma model selection", names(col)[i])) col[i] = "red"
	if(grepl("oracle", names(col)[i])) col[i] = "purple"
	if(grepl("grand null", names(col)[i])) col[i] = "black"
	if(grepl("grand alternative", names(col)[i])) col[i] = "grey50"
	if(grepl("backward elimination", names(col)[i])) col[i] = "grey40"
	if(names(col)[i]=="lasso") col[i] = "green3"
	if(grepl("lasso model selection", names(col)[i])) col[i] = "green2"
	if(names(col)[i]=="elastic net") col[i] = "yellow3"
	if(grepl("elastic net model selection", names(col)[i])) col[i] = "yellow"
	if(names(col)[i]=="ridge regression") col[i] = "orange"
}
subana.lab = gsub(".", "", ifelse(subana=="", "Main", subana), fixed=TRUE)

###################
# Plotting params #
###################
plot.ready = function(filename=NULL, width_height_inches=3.5, pointsize=8) {
	if(!is.null(filename)) {
		filename = file.path(outdir, filename)
		quartz(file=filename, type="pdf", width=2*width_height_inches, height=width_height_inches, pointsize=pointsize)
	}
	par(mar=c(5,12.5,1,1), mfrow=c(1,2))
}
plot.done = function() {
	if(names(dev.cur())!="quartz") dev.off()
}
plot.panel_lab = function(lab='A', width_height_inches=3.5, width_height_inches_label=0.1) {
	if(lab=='A') {
		x_user = grconvertX(width_height_inches_label, from="inches", to="user")
		y_user = grconvertY(width_height_inches-width_height_inches_label, from="inches", to="user")
		text(x_user, y_user, labels="A", font=2, cex=1.5, xpd=TRUE)
	} else if(lab=='B') {
		x_user = grconvertX(width_height_inches+width_height_inches_label, from="inches", to="user")
		y_user = grconvertY(width_height_inches-width_height_inches_label, from="inches", to="user")
		text(x_user, y_user, labels="B", font=2, cex=1.5, xpd=TRUE)
	} else {
		stop("plot.panel_lab(): lab must be 'A' or 'B'")
	}
}

#########################################
# Figure 3: Prediction error and legend #
#########################################
##############
# Prediction #
##############
# Compare prediction L2
test = unlist(lapply(res, function(RES) mean(RES@prediction.L2)))
test = test[!is.na(test)]
XLAB = "Prediction L2 error"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

plot.ready("figure_3_prediction_error_and_legend.pdf")
XLIM = c(0, max(test))
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="", xaxs="r")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])
plot.panel_lab("A")

##################
# Methods legend #
##################

legend.text = c(
"Model selection",
"Oracle",
"Grand alternative",
"Doublethink",
"LASSO",
"Elastic net",
"MR-BMA",
"Backward elimination",
"Grand null",

"Machine learning",
NA,
"Ridge regression",
NA,
"LASSO",
"Elastic net",
NA,
NA,
NA,

"BMA",
NA,
NA,
"Doublethink",
NA,
NA,
"MR-BMA",
NA,
NA)

legend.col = col[c(
NA,
"oracle",
"grand alternative",
"doublethink model selection mu = 0.05; h = 1",
"lasso model selection",
"elastic net model selection",
"mr-bma model selection",
"backward elimination",
"grand null",

NA,
NA,
"ridge regression",
NA,
"lasso",
"elastic net",
NA,
NA,
NA,

NA,
NA,
NA,
"doublethink bma mu = 0.05; h = 1",
NA,
NA,
"mr-bma bma",
NA,
NA)]
legend.pch = ifelse(is.na(legend.text), NA, 22)
legend.pch[match(c("Model selection", "Machine learning", "BMA"), legend.text)] = NA
legend.font = 1; legend.font[match(c("Model selection", "Machine learning", "BMA"), legend.text)] = 2

wh = which(!is.na(legend.text))
addspace = c(1:9, NA, 10:13, NA, 14:16)
par(mar=rep(0.1,4))
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
suppressWarnings(legend("center", legend.text[wh][addspace], ncol=1, bg="white", fill=legend.col[wh][addspace], border=1-is.na(legend.col[wh][addspace]), text.font=legend.font[wh][addspace], bty="n", xpd=TRUE, cex=1.4))
plot.panel_lab("B")
plot.done()

################################
# Figure 4                     #
# Point and interval estimates #
################################
# Compare estimator L2
test = unlist(lapply(res, function(RES) RES@estimate.L2))
test = test[!is.na(test)]
XLAB = "Estimate mean L2 error"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

plot.ready("figure_4_estimator_error_and_stderr_coverage.pdf")
XLIM = c(0, max(test))
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="", xaxs="r")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])
plot.panel_lab("A")

# Compare coverage
test = unlist(lapply(res, function(RES) mean(RES@stderr.coverage)))
test = test[!is.na(test)]
XLAB = paste0(floor(100*(1-alpha)), "% Standard error coverage")

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

XLIM = c(0.9, 1)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="", xaxs="r")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
abline(v=se(1-alpha), lty=1, col=c("grey80", "black", "grey80"))
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])
plot.panel_lab("B")
plot.done()

##############################################
# Figure 5: marginal FWER and strikeout rate #
##############################################
#################
# Marginal FWER #
#################
test = c(unlist(lapply(res, function(RES) RES@freqt.marginal@familywiseI)),
        unlist(lapply(res, function(RES) lapply(RES@pvalueTests.marginal, function(PVT) PVT@familywiseI))))
test = test[!is.na(test)]
names(test) <- gsub("\\.1$", "", names(test))
XLAB = "Marginal type I FWER"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

plot.ready("figure_5_marginal_fwerI_and_strikeoutII.pdf")
XLIM = c(alpha/10, 1)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="i")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
abline(v=se(alpha), lty=1, col=c("grey80", "black", "grey80"))
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
plot.panel_lab("A")

#######################
# Marginal Strikeouts #
#######################
test = c(unlist(lapply(res, function(RES) RES@freqt.marginal@strikeoutII)),
		unlist(lapply(res, function(RES) lapply(RES@pvalueTests.marginal, function(PVT) PVT@strikeoutII))))
test = test[!is.na(test)]
names(test) <- gsub("\\.1$", "", names(test))
XLAB = "Marginal type II strikeout"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

XLIM = range(test[!is.na(wh.ana)])
XLIM = c(0.02, 0.2)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="r")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
plot.panel_lab("B")
plot.done()


##############################################
# Figure 6: pairwise FWER and strikeout rate #
##############################################
#################
# Pairwise FWER #
#################
test = c(unlist(lapply(res, function(RES) RES@freqt.pairwise@familywiseI)),
        unlist(lapply(res, function(RES) lapply(RES@pvalueTests.pairwise, function(PVT) PVT@familywiseI))))
test = test[!is.na(test)]
names(test) <- gsub("\\.1$", "", names(test))
XLAB = "Pairwise type I FWER"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

plot.ready("figure_6_pairwise_fwerI_and_strikeoutII.pdf")
XLIM = c(alpha/10, 1)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="i")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
abline(v=se(alpha), lty=1, col=c("grey80", "black", "grey80"))
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
plot.panel_lab("A")

#######################
# Pairwise Strikeouts #
#######################
test = c(unlist(lapply(res, function(RES) RES@freqt.pairwise@strikeoutII)),
		unlist(lapply(res, function(RES) lapply(RES@pvalueTests.pairwise, function(PVT) PVT@strikeoutII))))
test = test[!is.na(test)]
names(test) <- gsub("\\.1$", "", names(test))
XLAB = "Pairwise type II strikeout"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

XLIM = range(test[!is.na(wh.ana)])
XLIM = c(0.02, 0.2)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="r")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
plot.panel_lab("B")
plot.done()


##############################################
# Figure 7: headline FWER and strikeout rate #
##############################################
#################
# Headline FWER #
#################
test = c(unlist(lapply(res, function(RES) RES@freqt.headline@familywiseI)),
        unlist(lapply(res, function(RES) lapply(RES@pvalueTests.headline, function(PVT) PVT@familywiseI))))
test = test[!is.na(test)]
names(test) <- gsub("\\.1$", "", names(test))
XLAB = "Headline type I FWER"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

plot.ready("figure_7_headline_fwerI_and_strikeoutII.pdf")
XLIM = c(alpha/10, 1)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="i")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
abline(v=se(alpha), lty=1, col=c("grey80", "black", "grey80"))
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
plot.panel_lab("A")


#######################
# Headline Strikeouts #
#######################
test = c(unlist(lapply(res, function(RES) RES@freqt.headline@strikeoutII)),
		unlist(lapply(res, function(RES) lapply(RES@pvalueTests.headline, function(PVT) PVT@strikeoutII))))
test = test[!is.na(test)]
names(test) <- gsub("\\.1$", "", names(test))
XLAB = "Headline type II strikeout"

wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
wh.ana.pos = rk(ana.pos[wh.ana])
wh.ana.lab = ana.lab[wh.ana]

XLIM = range(test[!is.na(wh.ana)])
XLIM = c(0.02, 0.2)
plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="r")
box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
#abline(v=se(alpha*tau), lty=1, col=c("grey80", "black", "grey80"))
points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
plot.panel_lab("B")
plot.done()


###########################################
# Figure 8: truenull FWER_1 and FWER_0.43 #
###########################################
#################
# Truenull FWER #
#################

fwer.rho.gradient = unname(gsub("fwer_", "", sapply(names(res[[1]]@pvalueTests.truenull$Bonf@typeI), function(s) unlist(strsplit(s, " "))[1])))

plot.ready("figure_8_fwerI_rho_1_vs_0.43_.pdf")
for(fwer.rho in c("1", "0.43")) {
	test = c(unlist(lapply(res, function(RES) unname(RES@freqt.truenull@typeI[grep(paste0("fwer_", fwer.rho), names(RES@freqt.truenull@typeI))]))),
			unlist(lapply(res, function(RES) lapply(RES@pvalueTests.truenull, function(PVT) unname(PVT@typeI[grep(paste0("fwer_", fwer.rho), names(PVT@typeI))])))))
	test = test[!is.na(test)]
	names(test) <- gsub("\\.1$", "", names(test))
	XLAB = paste("Type I FWER ρ", fwer.rho)

	wh.ana = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"row"]))
	wh.sub = sapply(names(test), function(TEST) len02na(which(TEST==ana.sub, arr=TRUE)[,"col"]))
	wh.ana.pos = rk(ana.pos[wh.ana])
	wh.ana.lab = ana.lab[wh.ana]
	# tp = sort(test[!is.na(wh.sub)]); tp[grep("oracle", names(tp))]

	XLIM = c(alpha/10, 1)
	plot(XLIM, range(wh.ana.pos, na.rm=TRUE), type="n", xlab=XLAB, ylab="", axes=FALSE, log="x", xaxs="i")
	box(); axis(1); axis(2, wh.ana.pos[!duplicated(wh.ana.pos)], wh.ana.lab[!duplicated(wh.ana.pos)], las=2)
	abline(v=se(alpha), lty=1, col=c("grey80", "black", "grey80"))
	points(pmin(XLIM[2], pmax(XLIM[1], test)), wh.ana.pos, pch=sub.pch[wh.sub], col=col[wh.ana], cex=sub.cex[wh.sub], lwd=sub.lwd[wh.sub])

	legend("right", subana.lab, pch=sub.pch, bg="white", pt.cex=sub.cex, pt.lwd=sub.lwd)
	if(fwer.rho=="1") plot.panel_lab("A")
	if(fwer.rho=="0.43") plot.panel_lab("B")
}
plot.done()
