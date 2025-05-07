#############
# Arguments #
#############
# Local clone of https://github.com/danny-wilson/dblthk_amd_sims
dblthk_amd_sims_repo = "~/Documents/GitHub/dblthk_amd_sims"
# Adjust for the 49 or 15 biomarkers?
number_of_tests = c("full_data", "partial_data")[1]
# Disallow the inclusion of variables in high correlation with others already included
#rsq.max.thresh = 0.2   # Three pairs covering 37/49 variables
rsq.max.thresh = 0.66  # Reduced gap in the distribution of rsq below the two most correlated pairs
output_filename = paste0("~/Downloads/table_2_mr_amd_subanalysis_rsqmax", rsq.max.thresh, ".csv")
overwrite = FALSE

######################
# Load code and data #
######################
# Get working directory to revert to after sourcing the GitHub files
wd = getwd()
setwd(dblthk_amd_sims_repo)
# Load required source code: this in turn sources summary_mvMR_BF.R
source("biomarker-sim-functions.R")
# Load example AMD data
full.data = load.data()
# Revert to user working directory
setwd(wd)

####################################
# Filter variables by correlations #
####################################
# Perform univariable association tests
results.univariable.MR = univariable.MR(full.data)

# Filter example data to the 15 most significant risk factors based on univariable associations
top.data = filter.data(full.data, as.integer(15), results.univariable.MR)
top.data@x.names
# [1] "ApoA1"      "HDL.C"      "HDL.D"      "IDL.TG"     "L.HDL.C"    "L.VLDL.C"  
# [7] "L.VLDL.TG"  "M.VLDL.C"   "M.VLDL.TG"  "S.HDL.TG"   "S.VLDL.C"   "S.VLDL.TG" 
#[13] "Serum.TG"   "XL.HDL.C"   "XS.VLDL.TG"

# Visualize the correlations between variables
dst = dist(t(full.data@x))
hc = hclust(dst)
od = hc$order
rsq = cor(full.data@x)^2
image(1:full.data@m, 1:full.data@m, rsq[od,od]<0.8, axes=FALSE, xlab="", ylab="")
axis(1, 1:full.data@m, full.data@x.names[od], las=2)
axis(2, 1:full.data@m, full.data@x.names[od], las=2)

# Introduce variables in descending order of univariate association
# Disallow the inclusion of variables in high correlation with others already included
data = filter.correlated.data(full.data, m=as.integer(15), rsq.max.thresh = rsq.max.thresh, results.univariable.MR=results.univariable.MR)
stopifnot(data@m==15)
keep = match(data@x.names, full.data@x.names)

# Strength of univariable associations (-log10 p-values adjusted for 49 tests)
tp = attributes(results.univariable.MR)$signif.neglog10padj; names(tp) <- attributes(results.univariable.MR)$names; sort(tp[keep], decreasing=TRUE)
#  XL.HDL.C    L.HDL.C   S.HDL.TG  S.VLDL.TG      ApoA1     IDL.TG      LDL.D 
# 4.3689996  4.1160842  2.0526599  1.2706733  0.6903025  0.0247417 -0.2116284 
#       Ace  XL.HDL.TG       ApoB     VLDL.D    M.HDL.C        His        Ala 
#-0.2960249 -0.4656011 -0.5870546 -0.6219905 -0.6548148 -0.6668876 -0.6686777 
#       Gln 
#-0.8895398 

# What proportion of variation is explained by the first 15 variables?
sum(colMeans((data@x %*% diag(data@m))^2))/
sum(colMeans((full.data@x %*% diag(full.data@m))^2))
# [1] 0.3406715

# Visualize squared correlations
hist(rsq[keep, keep][lower.tri(rsq[keep, keep])], 50)
# Check there are three pairs of highly correlated variables
(tp = rowSums(rsq[keep, keep] > rsq.max.thresh))
# Error if size of correlated clusters exceeds 2
stopifnot(max(tp)<=2)

# List the filtered variable names
cat(data@x.names)
# 0.20:  XL.HDL.C L.HDL.C IDL.TG LDL.D Ace ApoB His Ala Gln Gly Ile Cit Tyr AcAce Leu
# 0.66:  XL.HDL.C L.HDL.C S.HDL.TG S.VLDL.TG ApoA1 IDL.TG LDL.D Ace XL.HDL.TG ApoB VLDL.D M.HDL.C His Ala Gln

# Display the top correlation coefficients
r.corr = cor(full.data@x)
tp = which(rsq[keep, keep]>=rsq.max.thresh, arr=TRUE); cbind(tp, colnames(rsq)[keep][tp[,2]], r.corr[keep,keep][tp])[tp[,1]<tp[,2],]
# rsq.max.thresh = 0.2:
#XL.HDL.C "1"  "2"  "L.HDL.C" "0.815187508469619"
#IDL.TG   "3"  "6"  "ApoB"    "0.913316563074477"
#Ile      "11" "15" "Leu"     "0.853081394495008"
#
# rsq.max.thresh = 0.66
#          row col
# XL.HDL.C "1" "2"  "L.HDL.C"   "0.815187508469618"
# S.HDL.TG "3" "4"  "S.VLDL.TG" "0.921914794378598"
# IDL.TG   "6" "10" "ApoB"      "0.913316563074477"
# Calculate the maximum absolute correlation coefficient below the threshold
sqrt(max(rsq[keep,keep][rsq[keep,keep]<rsq.max.thresh]))
# [1] 0.7981375

# Top squared correlation coefficients (3 highly correlated pairs):
head(sort(rsq[keep, keep][lower.tri(rsq[keep, keep])], decreasing=TRUE))
# Top absolute correlation coefficients (3 highly correlated pairs):
sqrt(head(sort(rsq[keep, keep][lower.tri(rsq[keep, keep])], decreasing=TRUE)))

#########################
# Run the data analysis #
#########################
nu = full.data@m
if(number_of_tests == "partial_data") nu = data@m
system.time((db = doublethink.x(data, h=1, mu=0.1, nu=nu)))
#   user  system elapsed 
# 46.224   0.902  47.510 
# Extract the results of the Doublethink analysis
dba = db$doublethink.bma$"mu = 0.1; h = 1"

#######################
# Inspect the results #
#######################
# NB: when analysing just a subset see Theorem 5 for valid interpretation of sub-analyses.
#############################
# Headline test: grand null #
#############################
tau = 9; mu = 0.1; h = 1; n = data@n; xi = h/(n+h)
### Sub-analysis ###
# Posterior odds against any of the 15 variables having an effect
10^dba@headline.signif.log10po
# Prior odds against any of the 15 variables having an effect
PrO = (1-(1-mu)^15)/(1-mu)^15
# Corresponding Bayes factor
10^dba@headline.signif.log10po/PrO
# Unadjusted p-value (Theorem 2)
pchisq(2*log(10^dba@headline.signif.log10po/((1+mu*sqrt(xi))^15-1)), 1, low=F) # preciser denom
pchisq(2*log(10^dba@headline.signif.log10po/(15*mu*sqrt(xi))), 1, low=F)       # theorem 2
### Full analysis ###
# Prior odds against any of the 49 variables having an effect
PrO = (1-(1-mu)^49)/(1-mu)^49
# Corresponding Bayes factor
10^dba@headline.signif.log10po/PrO
# Unadjusted p-value (Theorem 2)
pchisq(2*log(10^dba@headline.signif.log10po/((1+mu*sqrt(xi))^49-1)), 1, low=F) # preciser denom
pchisq(2*log(10^dba@headline.signif.log10po/(49*mu*sqrt(xi))), 1, low=F)       # theorem 2
# Adjusted p-value (Theorem 3)
pchisq(2*log(10^dba@headline.signif.log10po/((1+mu*sqrt(xi))^nu-1)), 1, low=F) # preciser denom
pchisq(2*log(10^dba@headline.signif.log10po/(nu*mu*sqrt(xi))), 1, low=F)       # theorem 3
# Cross-check adjusted p-value with doublethink output
10^-dba@headline.signif.neglog10padj

##################
# Marginal tests #
##################
### Sub-analysis ###
# Posterior odds against each of the 15 variables having an effect
10^dba@signif.log10po
# Corresponding Bayes factor
10^dba@signif.log10po/mu
# Unadjusted p-value (Theorem 2)
pchisq(2*log(10^dba@signif.log10po/(1*mu*sqrt(xi))), 1, low=F)       # theorem 2
### Full analysis ###
# Prior odds against any of the 35 variables having an effect
PrO = (1-(1-mu)^35)/(1-mu)^35
# Corresponding Bayes factor
10^dba@signif.log10po/PrO
# Unadjusted p-value (Theorem 2)
pchisq(2*log(10^dba@signif.log10po/(35*mu*sqrt(xi))), 1, low=F)       # theorem 2
# Adjusted p-value (Theorem 3)
pchisq(2*log(10^dba@signif.log10po/(49*mu*sqrt(xi))), 1, low=F)       # theorem 3
# Cross-check adjusted p-value with doublethink output
10^-dba@signif.neglog10padj

##################
# Pairwise tests #
##################
### Sub-analysis ###
# Posterior odds against each of the 15 variables having an effect
head(sort(10^dba@pairwise.signif.log10po, decreasing=TRUE))
# Prior odds against either of the variables having an effect
PrO = (1-(1-mu)^2)/(1-mu)^2
# Corresponding Bayes factor
head(sort(10^dba@pairwise.signif.log10po, decreasing=TRUE))/PrO
# Unadjusted p-value (Theorem 2)
pchisq(2*log(head(sort(10^dba@pairwise.signif.log10po, decreasing=TRUE))/(2*mu*sqrt(xi))), 1, low=F)       # theorem 2
### Full analysis ###
# Prior odds against any of the 36 variables having an effect
PrO = (1-(1-mu)^36)/(1-mu)^36
# Corresponding Bayes factor
head(sort(10^dba@pairwise.signif.log10po, decreasing=TRUE))/PrO
# Unadjusted p-value (Theorem 2)
pchisq(2*log(head(sort(10^dba@pairwise.signif.log10po, decreasing=TRUE))/(36*mu*sqrt(xi))), 1, low=F)       # theorem 2
# Adjusted p-value (Theorem 3)
pchisq(2*log(head(sort(10^dba@pairwise.signif.log10po, decreasing=TRUE))/(49*mu*sqrt(xi))), 1, low=F)       # theorem 3
# Cross-check adjusted p-value with doublethink output
head(sort(10^-dba@pairwise.signif.neglog10padj, decreasing=FALSE))

# Compare univariable and multivariable p-values
log10(drop1(lm(amd ~ 0 + XL.HDL.C, data=data.frame(full.data@y, full.data@x)), test="Chisq")$"Pr(>Chi)")
#[1]        NA -6.115903
log10(drop1(lm(amd ~ 0 + L.HDL.C, data=data.frame(full.data@y, full.data@x)), test="Chisq")$"Pr(>Chi)")
#[1]        NA -5.860463
log10(drop1(lm(amd ~ 0 + XL.HDL.C + L.HDL.C, data=data.frame(full.data@y, full.data@x)), test="Chisq")$"Pr(>Chi)")
#[1]         NA -1.1488410 -0.8397098

#########################
# Save results to Table #
#########################
# Output unsorted results to plain text file
# Column names:
#   po       Posterior odds of including each variable
#   bf       Corresponding Bayes factor
#   pp       Corresponding posterior probability
#   p.adj    p-value adjusted for the number of tests (Theorem 3). p > 0.02 --> p = 1
#   p.unadj  p-value unadjusted for the number of tests (Theorem 2)
#   p.Bonf   Bonferroni correction for the unadjusted p-value
#   pe       Point estimate (posterior mean)
#   stderr   Standard error (root posterior variance)
if(!file.exists(output_filename) | overwrite==TRUE) write.csv(data.frame('po'=10^dba@signif.log10po, 'bf'=10^dba@signif.log10po/mu, 'pp'=10^dba@signif.log10po/(1+10^dba@signif.log10po), 'p.adj'=10^-dba@signif.neglog10padj, 'p.unadj'=pchisq(2*log(10^dba@signif.log10po/mu/sqrt(xi)), 1, low=FALSE), 'p.Bonf'=pchisq(2*log(10^dba@signif.log10po/mu/sqrt(xi)), 1, low=FALSE)*nu, 'pe'=dba@estimate, 'stderr'=dba@stderror), output_filename)

#############################
# p-value combination tests #
#############################
# Marginal-based combined tests, based on Theorem 2.
# Do any of them yield significance, headline or pairwise (L.HDL.C | XL.HDL.C) ?
# NB Theorem 5: interpretation of sub-analyses

# Bonferroni - NEITHER HEADLINE NOR PAIRWISE
10^-dba@pvalueTests$Bonf@marginal.neglog10padj
#   XL.HDL.C    L.HDL.C
# 0.05242055 0.17402875
head(sort(10^-dba@pvalueTests$Bonf@pairwise.neglog10padj, decreasing=FALSE), n=20)
# L.HDL.C | XL.HDL.C  S.HDL.TG | XL.HDL.C S.VLDL.TG | XL.HDL.C     ApoA1 | XL.HDL.C    IDL.TG | XL.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.05242055
#   LDL.D | XL.HDL.C       Ace | XL.HDL.C XL.HDL.TG | XL.HDL.C      ApoB | XL.HDL.C    VLDL.D | XL.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.05242055
# M.HDL.C | XL.HDL.C       His | XL.HDL.C       Ala | XL.HDL.C       Gln | XL.HDL.C   S.HDL.TG | L.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.17402875
10^-dba@pvalueTests$Bonf@headline.neglog10padj
# 0.05242055

# HMP - NEITHER
10^-dba@pvalueTests$HMP@marginal.neglog10padj
#   XL.HDL.C    L.HDL.C
# 0.07289146 0.38426250
head(sort(10^-dba@pvalueTests$HMP@pairwise.neglog10padj, decreasing=FALSE), n=20)
# L.HDL.C | XL.HDL.C  S.HDL.TG | XL.HDL.C S.VLDL.TG | XL.HDL.C     ApoA1 | XL.HDL.C
#         0.05243803           0.07289146           0.07289146           0.07289146
10^-dba@pvalueTests$HMP@headline.neglog10padj
# 0.05243803

# Simes - NEITHER
pComb = calc.Simes(p.th2, nu)
10^-dba@pvalueTests$Simes@marginal.neglog10padj
#   XL.HDL.C     L.HDL.C
# 0.05242055  0.17402875
head(sort(10^-dba@pvalueTests$Simes@pairwise.neglog10padj, decreasing=FALSE), n=20)
# L.HDL.C | XL.HDL.C  S.HDL.TG | XL.HDL.C S.VLDL.TG | XL.HDL.C     ApoA1 | XL.HDL.C    IDL.TG | XL.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.05242055
#   LDL.D | XL.HDL.C       Ace | XL.HDL.C XL.HDL.TG | XL.HDL.C      ApoB | XL.HDL.C    VLDL.D | XL.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.05242055
# M.HDL.C | XL.HDL.C       His | XL.HDL.C       Ala | XL.HDL.C       Gln | XL.HDL.C   S.HDL.TG | L.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.17402875
10^-dba@pvalueTests$Simes@headline.neglog10padj
# 0.05242055

# Hommel - NEITHER
10^-dba@pvalueTests$Hommel@marginal.neglog10padj
#   XL.HDL.C     L.HDL.C
# 0.05242055  0.17047714
head(sort(10^-dba@pvalueTests$Hommel@pairwise.neglog10padj, decreasing=FALSE), n=20)
# L.HDL.C | XL.HDL.C  S.HDL.TG | XL.HDL.C S.VLDL.TG | XL.HDL.C     ApoA1 | XL.HDL.C    IDL.TG | XL.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.05242055
#   LDL.D | XL.HDL.C       Ace | XL.HDL.C XL.HDL.TG | XL.HDL.C      ApoB | XL.HDL.C    VLDL.D | XL.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.05242055
# M.HDL.C | XL.HDL.C       His | XL.HDL.C       Ala | XL.HDL.C       Gln | XL.HDL.C   S.HDL.TG | L.HDL.C
#         0.05242055           0.05242055           0.05242055           0.05242055           0.17047714
10^-dba@pvalueTests$Hommel@headline.neglog10padj
# 0.05242055
