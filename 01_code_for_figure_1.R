## code to produce Figure 1 and the corresponding numerical study

## BEGIN SETUP ##

## load libraries
require(foreach)
require(doParallel)
require(doSNOW)

## set parameters
K <- c(seq(0.01, 1.01, 0.1), seq(2,5)) ## ratio of variability in group 2 vs. group 1
rho_lam <- seq(-0.975, 0.975, 0.05) ## correlation between two delta components

S1 <- c(seq(0.01, 1.01, 0.1), seq(2,5)) ## signal-to-noise ratio for group 1
S2 <- c(seq(0.01, 1.01, 0.1), seq(2,5)) ## signal-to-noise ratio for group 2
rho_sig <- seq(-0.975, 0.975, 0.05) ## correlation between two error components

## create a matrix with all possible combinations of parameters
params <- expand.grid(K = K, rho_lam = rho_lam, S1 = S1, S2 = S2, rho_sig = rho_sig)

## set up parallelization
cores=detectCores()
cl <- makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb <- txtProgressBar(max = 0.1*nrow(params), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## split results into 10 matrices to reduce memory
for (j in 1:10){
  print(j)
  results <- foreach(k=1:(0.1*nrow(params)), .combine=rbind, .packages = c("MASS", "mvtnorm"),
                     .options.snow=opts) %dopar% {
                       
                       ## code to extract from the correct row of the matrix
                       i <- (j-1)*0.1*nrow(params) + k
                       params_i <- as.numeric(unlist(params[i,]))
                       
                       ## extract 5 components from corresponding row
                       K_i <- params_i[1]; rho_lam_i <- params_i[2]
                       S1_i <- params_i[3]; S2_i <- params_i[4]; rho_sig_i <- params_i[5]
                       
                       ## compute Lambda and Sigma matrices corresponding to these inputs
                       Lambda_i <- matrix(c(1, sqrt(K_i)*rho_lam_i,
                                            sqrt(K_i)*rho_lam_i, K_i), nrow = 2)
                       
                       Sigma_i <- matrix(c(1/S1_i, sqrt(K_i/(S1_i*S2_i))*rho_sig_i,
                                           sqrt(K_i/(S1_i*S2_i))*rho_sig_i, K_i/S2_i), nrow = 2)
                       
                       ## compute univariate and bivariate variances using the formulas
                       ## from the paper
                       var_uni_i <- sum(Lambda_i)*sum(Sigma_i)/(sum(Lambda_i) + sum(Sigma_i))
                       
                       var_bi_i <- sum(solve(solve(as.matrix(Lambda_i)) + solve(Sigma_i)))
                       
                       ## add univariate and bivariate variances to the row and return
                       c(params_i, var_uni_i, var_bi_i)
                     }
  ## output results to a .csv file
  write.csv(results,paste0("two_components_",j,".csv"), row.names = FALSE)
}

## compile results from 10 .csv files into res
res <- NULL
for (j in 1:10){
  print(j)
  res <- rbind(res, read.csv(paste0("two_components_", j, ".csv")))
}

## column 8 is the ratio of variances (univariate / bivariate >= 1)
res[,8] <- res[,6]/res[,7]

## univariate summaries based on segmentation

## we segment into 8 categories based on lambda_rho (res[,2]), lambda_Sigma (res[,5])
## along with the difference between S_1 and S_2 (res[,3] - res[,4])
sub1 <- subset(res, abs(res[,2]) > 0.8 & abs(res[,5]) > 0.8 & abs(res[,3] - res[,4]) > 2)
sub2 <- subset(res, abs(res[,2]) > 0.8 & abs(res[,5]) > 0.8 & abs(res[,3] - res[,4]) < 0.1)
sub3 <- subset(res, abs(res[,2]) > 0.8 & abs(res[,5]) < 0.1 & abs(res[,3] - res[,4]) > 2)
sub4 <- subset(res, abs(res[,2]) > 0.8 & abs(res[,5]) < 0.1 & abs(res[,3] - res[,4]) < 0.1)
sub5 <- subset(res, abs(res[,2]) < 0.1 & abs(res[,5]) > 0.8 & abs(res[,3] - res[,4]) > 2)
sub6 <- subset(res, abs(res[,2]) < 0.1 & abs(res[,5]) > 0.8 & abs(res[,3] - res[,4]) < 0.1)
sub7 <- subset(res, abs(res[,2]) < 0.1 & abs(res[,5]) < 0.1 & abs(res[,3] - res[,4]) > 2)
sub8 <- subset(res, abs(res[,2]) < 0.1 & abs(res[,5]) < 0.1 & abs(res[,3] - res[,4]) < 0.1)

## we first consider combinations where S_1 and S_2 are mismatched (i.e., differ by more than 2
## in absolute value). These are combinations 1, 3, 5, and 7 from earlier. The variate of interest
## is the ratio of univariate to bivariate variances (given by the 8th column).
df.plot <- data.frame(factor = c(sub1[,8], sub3[,8], sub5[,8], sub7[,8]),
                      group = c(rep(1, length(sub1[,8])), rep(2, length(sub3[,8])),
                                rep(3, length(sub5[,8])), rep(4, length(sub7[,8]))),
                      med = round(c(rep(median(sub1[,8]), length(sub1[,8])), rep(median(sub3[,8]), length(sub3[,8])),
                                    rep(median(sub5[,8]), length(sub5[,8])), rep(median(sub7[,8]), length(sub7[,8]))), 2))

## create more informative labels for the subplots
df.plot$group <- factor(df.plot$group, levels = c(1, 2, 3, 4),
                        labels = c(expression("|"*rho[lambda]~"| > 0.8, |"*rho[Sigma]~"| > 0.8"), 
                                   expression("|"*rho[lambda]~"| > 0.8, |"*rho[Sigma]~"| < 0.1"), 
                                   expression("|"*rho[lambda]~"| < 0.1, |"*rho[Sigma]~"| > 0.8"), 
                                   expression("|"*rho[lambda]~"| < 0.1, |"*rho[Sigma]~"| < 0.1"))
)

## generate Figure 1
p <- ggplot(data = df.plot, aes(x=factor)) + geom_density(fill="firebrick", alpha = 0.5) + theme_bw() + xlim(0,11)
p <- p + facet_wrap( ~ group, labeller = label_parsed) + xlab("\nVariance Reduction Factor") +
  ylab("Density\n") + 
  geom_segment(aes(x = med, xend = med, y =0, yend = Inf, group = group), colour = 'black', lty = 2) +
  geom_text(
    data    = df.plot,
    mapping = aes(x = med + 0.1, y = 3, label = med),
    hjust   = -0.1,
    vjust   = -1
  ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(strip.text.x = element_text(size = 10))

## output to pdf
pdf(file = "SNRMisMatch.pdf",   # The directory you want to save the file in
    width = 6, 
    height = 4.5) 
p
dev.off()

## we also consider combinations where S_1 and S_2 are not mismatched (i.e., differ by less than 0.1
## in absolute value). These are combinations 2, 4, 6, and 8 from earlier. We ignore combination 8 since 
## its density curve is essentially a point mass at 1. The variate of interest
## is the ratio of univariate to bivariate variances (given by the 8th column).
df.plot <- data.frame(factor = c(sub6[,8], sub4[,8], sub2[,8]),
                      group = c(rep(5, length(sub6[,8])), rep(6, length(sub4[,8])),
                                rep(7, length(sub2[,8]))),
                      med = round(c(rep(median(sub6[,8]), length(sub6[,8])), rep(median(sub4[,8]), length(sub4[,8])),
                                    rep(median(sub2[,8]), length(sub2[,8]))), 2))

## create more informative labels for the subplots
df.plot$group <- factor(df.plot$group, levels = c(5,6,7),
                        labels = c(expression("|"*rho[lambda]~"| < 0.1, |"*rho[Sigma]~"| > 0.8"), 
                                   expression("|"*rho[lambda]~"| > 0.8, |"*rho[Sigma]~"| < 0.1"), 
                                   expression("|"*rho[lambda]~"| > 0.8, |"*rho[Sigma]~"| > 0.8"))
)

## create additional figure that is referenced but not provided in the paper.
p <- ggplot(data = df.plot, aes(x=factor)) + geom_density(fill="firebrick", alpha = 0.5) + theme_bw() + xlim(0,5)
p <- p + facet_wrap( ~ group, labeller = label_parsed) + xlab("\nVariance Reduction Factor") +
  ylab("Density\n") + 
  geom_segment(aes(x = med, xend = med, y =0, yend = Inf, group = group), colour = 'black', lty = 2) +
  geom_text(
    data    = df.plot,
    mapping = aes(x = med + 0.1, y = 6.75, label = med),
    hjust   = -0.1,
    vjust   = -1
  ) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(strip.text.x = element_text(size = 10))
 
## output to pdf
pdf(file = "SNRMatch.pdf",   # The directory you want to save the file in
    width = 6, 
    height = 3) 
p
dev.off()