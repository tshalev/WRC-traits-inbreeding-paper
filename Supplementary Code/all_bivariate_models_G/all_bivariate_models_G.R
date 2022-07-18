rm(list = ls())

library(asreml)
library(extrafont)
library(asremlPlus)
library(synbreed)
library(sommer)
library(nadiv)
library(MCMCglmm)
library(data.table)
library(matrixcalc)
library(Matrix)
library(rrBLUP)
library(data.table)
library(plotrix)
library(lme4)
library(reshape2)
library(ggsci)
library(tidyverse)

source("publication_theme.r")

font_import()
loadfonts(device = "win")
windowsFonts()

varcomps <- function(x) {
  varcomp.df <- data.frame(summary(x)$varcomp)
  varcomp.df <- tibble::rownames_to_column(varcomp.df, "effect")
  varcomp.df <- varcomp.df %>%  # Add a column for percentage of all variance components
    mutate(vc.percent = (component/sum(component))*100)
  print(varcomp.df)
}

wrc.pheno.dat <- read.table("master_phenotype_dataset.csv", na.strings = "NA", header = T, sep = ",")

# Some basic plots
monoterpenes <- wrc.pheno.dat %>% 
  select(fath, fbth, fsab, fmyr, flim, fapi)

melted_monoterpenes <- melt(monoterpenes)


col <- pal_nejm()(8)

ggplot(melted_monoterpenes, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.4) +
  theme_Publication() +
  xlab(expression(bold(paste("Concentration (",mu,"g/g DW")))) +
  ylab("Density") +
  scale_fill_manual(name="Monoterpene",
                    # breaks=c("fath", "fbth", "fsab", "fmyr", "flim", "fapi"),
                    labels=c(expression(paste(alpha,"-thujone")), expression(paste(beta,"-thujone")), "Sabinene", "Myrcene",
                             expression(paste(italic(r),"-limonene")), expression(paste("(+)-",alpha,"-pinene"))),
                    values = col) +
  theme(legend.position = c(0.9, 0.9), legend.background = element_rect(colour = "black")) +
  scale_x_continuous(limits = c(0, 60000))

# ggsave("monoterpene_density_plot.svg", dpi = 300, height = 8, width = 10)


# Read in RRM
# Read in 012 file
alleles <- read.csv("training_samples.012",
                    header = F, sep = "\t", na.strings = -1)
str(alleles)
alleles <- alleles[,-1]

# Read in SNP ids
snps <- read.csv("training_samples.012.pos",
                 header = F, sep  = "\t")
snps$V1 <- as.character(snps$V1)
snps$V2 <- as.character(snps$V2)

# Merge columns and transpose
snps <- snps %>% 
  unite(col = snp, V1:V2, sep = '_')
t_snps <- t(snps)

# Read sample ids
indivs <- read.csv("training_samples.012.indv", header = F, sep = "\t")

# Put row and column names into 012 file
row.names(alleles) <- indivs$V1
colnames(alleles) <- t_snps

alleles <- alleles[order(as.numeric(row.names(alleles))),]
rownames(alleles[1:40,])

matchingSamples <- intersect(wrc.pheno.dat$tree, indivs$V1)

# Convert to matrix
g.mat <- as.matrix(alleles)

# Subtract to make -101
g.mat <- g.mat - 1
g.mat[1:10,1:2]

# Get imputed relationships
G <- A.mat(g.mat, return.imputed = T)
G$imputed[1:10,1:10]
G012 <- G$imputed + 1

# Make G matrix
G <- A.mat(g.mat)

# Make G near positive definite
require(Matrix)
RealizedPD <- nearPD(G, keepDiag = T)
G <- matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
G <- G + diag(0.01, nrow(G))
attr(G, "dimnames") = RealizedPD[[1]]@Dimnames
class(G) = "relationshipMatrix"

# Make G a relationship matrix
wrc.giv <- write.relationshipMatrix(G, sorting = "ASReml", type = "ginv", file = NULL)
head(attr(wrc.giv, "rowNames"))
names(wrc.giv) <- c("row", "column", "coefficient")
head(wrc.giv)
attr(wrc.giv, "INVERSE") <- TRUE

#Filter out non-genotyped trees
wrc.pheno.dat <- wrc.pheno.dat %>% 
  filter(tree %in% matchingSamples)

# Convert to factor, numeric
wrc.pheno.dat[1:14] <- lapply(wrc.pheno.dat[1:14], as.factor)
wrc.pheno.dat[15:74] <- lapply(wrc.pheno.dat[15:74], as.numeric)

### All bivariate models ######
# a-thujone - b-thujone
athuj_bthuj_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, fbth) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_bthuj_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, fbth) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_bthuj_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, fbth) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_bthuj_bivar_model2, athuj_bthuj_bivar_model1)
REMLRT.asreml(athuj_bthuj_bivar_model3, athuj_bthuj_bivar_model1)
REMLRT.asreml(athuj_bthuj_bivar_model2, athuj_bthuj_bivar_model3)


varcomps(athuj_bthuj_bivar_model1)
asreml::vpredict(athuj_bthuj_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.6502 +/-0.09800418
asreml::vpredict(athuj_bthuj_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.6032738 +/-0.01848261


# a-thujone - sabinene
athuj_sab_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ id(units):us(trait),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_sab_bivar_model2, athuj_sab_bivar_model1)
REMLRT.asreml(athuj_sab_bivar_model3, athuj_sab_bivar_model1)
REMLRT.asreml(athuj_sab_bivar_model2, athuj_sab_bivar_model3)
REMLRT.asreml(athuj_sab_bivar_model3, athuj_sab_bivar_model4)
REMLRT.asreml(athuj_sab_bivar_model3, athuj_sab_bivar_model5)


varcomps(athuj_sab_bivar_model3)
asreml::vpredict(athuj_sab_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.3464871 +/-0.1065419
asreml::vpredict(athuj_sab_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.655284 +/-0.01667956

# b-thujone - sabinene
bthuj_sab_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ id(units):us(trait),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_sab_bivar_model2, bthuj_sab_bivar_model1)
REMLRT.asreml(bthuj_sab_bivar_model3, bthuj_sab_bivar_model1)
REMLRT.asreml(bthuj_sab_bivar_model2, bthuj_sab_bivar_model3)
REMLRT.asreml(bthuj_sab_bivar_model3, bthuj_sab_bivar_model4)
REMLRT.asreml(bthuj_sab_bivar_model3, bthuj_sab_bivar_model5)


varcomps(bthuj_sab_bivar_model1)
asreml::vpredict(bthuj_sab_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2949608 +/-0.1371582
asreml::vpredict(bthuj_sab_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.3193808 +/-0.02764154

# a-thujone - a-thujaplicin
athuj_athuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_athuja_bivar_model2, athuj_athuja_bivar_model1)
REMLRT.asreml(athuj_athuja_bivar_model3, athuj_athuja_bivar_model1)
REMLRT.asreml(athuj_athuja_bivar_model2, athuj_athuja_bivar_model3)
REMLRT.asreml(athuj_athuja_bivar_model3, athuj_athuja_bivar_model4)
REMLRT.asreml(athuj_athuja_bivar_model3, athuj_athuja_bivar_model5)


varcomps(athuj_athuja_bivar_model3)
asreml::vpredict(athuj_athuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.108717 +/-0.1354607
asreml::vpredict(athuj_athuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.02666598 +/-0.02847733

# a-thujone - b-thujaplicin
athuj_bthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)#, G.param = athuj_bthuja_bivar_model1$G.param, R.param = athuj_bthuja_bivar_model1$R.param)

athuj_bthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_bthuja_bivar_model2, athuj_bthuja_bivar_model1)
REMLRT.asreml(athuj_bthuja_bivar_model3, athuj_bthuja_bivar_model1)
REMLRT.asreml(athuj_bthuja_bivar_model2, athuj_bthuja_bivar_model3)
REMLRT.asreml(athuj_bthuja_bivar_model3, athuj_bthuja_bivar_model4)
REMLRT.asreml(athuj_bthuja_bivar_model3, athuj_bthuja_bivar_model5)


varcomps(athuj_bthuja_bivar_model3)
asreml::vpredict(athuj_bthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1744438 +/-0.1298524
asreml::vpredict(athuj_bthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.01975459 +/-0.02779304

# a-thujone - g-thujaplicin
athuj_gthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06, G.param = athuj_gthuja_bivar_model1$G.param, R.param = athuj_gthuja_bivar_model1$R.param)

athuj_gthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_gthuja_bivar_model2, athuj_gthuja_bivar_model1)
REMLRT.asreml(athuj_gthuja_bivar_model3, athuj_gthuja_bivar_model1)
REMLRT.asreml(athuj_gthuja_bivar_model2, athuj_gthuja_bivar_model3)
REMLRT.asreml(athuj_gthuja_bivar_model3, athuj_gthuja_bivar_model4)
REMLRT.asreml(athuj_gthuja_bivar_model3, athuj_gthuja_bivar_model5)


varcomps(athuj_gthuja_bivar_model3)
asreml::vpredict(athuj_gthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.07063596 +/-0.1391504
asreml::vpredict(athuj_gthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.03290031 +/-0.02763803

# a-thujone - b-thujaplicinol
athuj_bthuja_ol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ id(units):us(trait),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_bthuja_ol_bivar_model2, athuj_bthuja_ol_bivar_model1)
REMLRT.asreml(athuj_bthuja_ol_bivar_model3, athuj_bthuja_ol_bivar_model1)
REMLRT.asreml(athuj_bthuja_ol_bivar_model2, athuj_bthuja_ol_bivar_model3)
REMLRT.asreml(athuj_bthuja_ol_bivar_model3, athuj_bthuja_ol_bivar_model4)
REMLRT.asreml(athuj_bthuja_ol_bivar_model3, athuj_bthuja_ol_bivar_model5)


varcomps(athuj_bthuja_ol_bivar_model3)
asreml::vpredict(athuj_bthuja_ol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.1988805 +/-0.158779
asreml::vpredict(athuj_bthuja_ol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.06059499 +/-0.02708283

# b-thujone - a-thujaplicin
bthuj_athuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06, G.param = bthuj_athuja_bivar_model3$G.param, R.param = bthuj_athuja_bivar_model3$R.param)

bthuj_athuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_athuja_bivar_model2, bthuj_athuja_bivar_model1)
REMLRT.asreml(bthuj_athuja_bivar_model3, bthuj_athuja_bivar_model1)
REMLRT.asreml(bthuj_athuja_bivar_model2, bthuj_athuja_bivar_model3)
REMLRT.asreml(bthuj_athuja_bivar_model3, bthuj_athuja_bivar_model4)
REMLRT.asreml(bthuj_athuja_bivar_model3, bthuj_athuja_bivar_model5)


varcomps(bthuj_athuja_bivar_model3)
asreml::vpredict(bthuj_athuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1240375 +/-0.1681627
asreml::vpredict(bthuj_athuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.03450621 +/-0.02839637

# b-thujone - b-thujaplicin
bthuj_bthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_bthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_bthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_bthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_bthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_bthuja_bivar_model2, bthuj_bthuja_bivar_model1)
REMLRT.asreml(bthuj_bthuja_bivar_model3, bthuj_bthuja_bivar_model1)
REMLRT.asreml(bthuj_bthuja_bivar_model2, bthuj_bthuja_bivar_model3)
REMLRT.asreml(bthuj_bthuja_bivar_model3, bthuj_bthuja_bivar_model4)
REMLRT.asreml(bthuj_bthuja_bivar_model3, bthuj_bthuja_bivar_model5)


varcomps(bthuj_bthuja_bivar_model3)
asreml::vpredict(bthuj_bthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1007333 +/-0.163554
asreml::vpredict(bthuj_bthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.0185997 +/-0.02725272

# b-thujone - g-thujaplicin
bthuj_gthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06, G.param = bthuj_gthuja_bivar_model1$G.param, R.param = bthuj_gthuja_bivar_model1$R.param)

bthuj_gthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_gthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_gthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_gthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_gthuja_bivar_model2, bthuj_gthuja_bivar_model1)
REMLRT.asreml(bthuj_gthuja_bivar_model3, bthuj_gthuja_bivar_model1)
REMLRT.asreml(bthuj_gthuja_bivar_model2, bthuj_gthuja_bivar_model3)
REMLRT.asreml(bthuj_gthuja_bivar_model3, bthuj_gthuja_bivar_model4)
REMLRT.asreml(bthuj_gthuja_bivar_model3, bthuj_gthuja_bivar_model5)


varcomps(bthuj_gthuja_bivar_model3)
asreml::vpredict(bthuj_gthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.02289319 +/-0.1738303
asreml::vpredict(bthuj_gthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.005025206 +/-0.02632923

# b-thujone - b-thujaplicinol
bthuj_bthuja_ol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06, G.param = bthuj_bthuja_ol_bivar_model1$G.param, R.param = bthuj_bthuja_ol_bivar_model1$R.param)

bthuj_bthuja_ol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ id(units):us(trait),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_bthuja_ol_bivar_model2, bthuj_bthuja_ol_bivar_model1)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model3, bthuj_bthuja_ol_bivar_model1)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model2, bthuj_bthuja_ol_bivar_model3)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model3, bthuj_bthuja_ol_bivar_model4)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model3, bthuj_bthuja_ol_bivar_model5)


varcomps(bthuj_bthuja_ol_bivar_model3)
asreml::vpredict(bthuj_bthuja_ol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.196501 +/-0.1905226
asreml::vpredict(bthuj_bthuja_ol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.02610444 +/-0.02516525


# sabinene - a-thujaplicin
sab_athuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_athuja_bivar_model2, sab_athuja_bivar_model1)
REMLRT.asreml(sab_athuja_bivar_model3, sab_athuja_bivar_model1)
REMLRT.asreml(sab_athuja_bivar_model2, sab_athuja_bivar_model3)
REMLRT.asreml(sab_athuja_bivar_model3, sab_athuja_bivar_model4)
REMLRT.asreml(sab_athuja_bivar_model3, sab_athuja_bivar_model5)


varcomps(sab_athuja_bivar_model3)
asreml::vpredict(sab_athuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.3097327 +/-0.1321742
asreml::vpredict(sab_athuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.02199255 +/-0.02901415

# sabinene - b-thujaplicin
sab_bthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_bthuja_bivar_model2, sab_bthuja_bivar_model1)
REMLRT.asreml(sab_bthuja_bivar_model3, sab_bthuja_bivar_model1)
REMLRT.asreml(sab_bthuja_bivar_model2, sab_bthuja_bivar_model3)
REMLRT.asreml(sab_bthuja_bivar_model3, sab_bthuja_bivar_model4)
REMLRT.asreml(sab_bthuja_bivar_model3, sab_bthuja_bivar_model5)


varcomps(sab_bthuja_bivar_model3)
asreml::vpredict(sab_bthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.05816861 +/-0.1231022
asreml::vpredict(sab_bthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.009707441 +/-0.02875327

# sabinene - g-thujaplicin
sab_gthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_gthuja_bivar_model2, sab_gthuja_bivar_model1)
REMLRT.asreml(sab_gthuja_bivar_model3, sab_gthuja_bivar_model1)
REMLRT.asreml(sab_gthuja_bivar_model2, sab_gthuja_bivar_model3)
REMLRT.asreml(sab_gthuja_bivar_model3, sab_gthuja_bivar_model4)
REMLRT.asreml(sab_gthuja_bivar_model3, sab_gthuja_bivar_model5)


varcomps(sab_gthuja_bivar_model3)
asreml::vpredict(sab_gthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.07420525 +/-0.1334569
asreml::vpredict(sab_gthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.02291048 +/-0.02773944

# sabinene - b-thujaplicinol
sab_bthuja_ol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                     residual = ~ id(units):us(trait),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_bthuja_ol_bivar_model2, sab_bthuja_ol_bivar_model1)
REMLRT.asreml(sab_bthuja_ol_bivar_model3, sab_bthuja_ol_bivar_model1)
REMLRT.asreml(sab_bthuja_ol_bivar_model2, sab_bthuja_ol_bivar_model3)
REMLRT.asreml(sab_bthuja_ol_bivar_model3, sab_bthuja_ol_bivar_model4)
REMLRT.asreml(sab_bthuja_ol_bivar_model3, sab_bthuja_ol_bivar_model5)


varcomps(sab_bthuja_ol_bivar_model3)
asreml::vpredict(sab_bthuja_ol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.009028014 +/-0.1516016
asreml::vpredict(sab_bthuja_ol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.02648834 +/-0.02676502

# a-thujone - ht15
athuj_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06, G.param = athuj_ht15_bivar_model3$G.param, R.param = athuj_ht15_bivar_model3$R.param)

athuj_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_ht15_bivar_model2, athuj_ht15_bivar_model1)
REMLRT.asreml(athuj_ht15_bivar_model3, athuj_ht15_bivar_model1)
REMLRT.asreml(athuj_ht15_bivar_model2, athuj_ht15_bivar_model3)
REMLRT.asreml(athuj_ht15_bivar_model3, athuj_ht15_bivar_model4)
REMLRT.asreml(athuj_ht15_bivar_model3, athuj_ht15_bivar_model5)


varcomps(athuj_ht15_bivar_model1)
asreml::vpredict(athuj_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.0310714 +/-0.160137
asreml::vpredict(athuj_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.05457004 +/-0.02907407

# a-thujone - dbh15
athuj_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_dbh15_bivar_model2, athuj_dbh15_bivar_model1)
REMLRT.asreml(athuj_dbh15_bivar_model3, athuj_dbh15_bivar_model1)
REMLRT.asreml(athuj_dbh15_bivar_model2, athuj_dbh15_bivar_model3)
REMLRT.asreml(athuj_dbh15_bivar_model3, athuj_dbh15_bivar_model4)
REMLRT.asreml(athuj_dbh15_bivar_model3, athuj_dbh15_bivar_model5)


varcomps(athuj_dbh15_bivar_model3)
asreml::vpredict(athuj_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1924222 +/-0.1747985
asreml::vpredict(athuj_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.097931359 +/-0.02843033

# b-thujone - ht15
bthuj_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_ht15_bivar_model2, bthuj_ht15_bivar_model1)
REMLRT.asreml(bthuj_ht15_bivar_model3, bthuj_ht15_bivar_model1)
REMLRT.asreml(bthuj_ht15_bivar_model2, bthuj_ht15_bivar_model3)
REMLRT.asreml(bthuj_ht15_bivar_model3, bthuj_ht15_bivar_model4)
REMLRT.asreml(bthuj_ht15_bivar_model3, bthuj_ht15_bivar_model5)


varcomps(bthuj_ht15_bivar_model1)
asreml::vpredict(bthuj_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2004948 +/-0.191928
asreml::vpredict(bthuj_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.04303263 +/-0.02749378

# b-thujone - dbh15
bthuj_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06, G.param = bthuj_dbh15_bivar_model1$G.param, R.param = bthuj_dbh15_bivar_model1$R.param)

bthuj_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

bthuj_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06, G.param = bthuj_dbh15_bivar_model3$G.param, R.param = bthuj_dbh15_bivar_model3$R.param)

bthuj_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

bthuj_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)



REMLRT.asreml(bthuj_dbh15_bivar_model2, bthuj_dbh15_bivar_model1)
REMLRT.asreml(bthuj_dbh15_bivar_model3, bthuj_dbh15_bivar_model1)
REMLRT.asreml(bthuj_dbh15_bivar_model2, bthuj_dbh15_bivar_model3)
REMLRT.asreml(bthuj_dbh15_bivar_model3, bthuj_dbh15_bivar_model4)
REMLRT.asreml(bthuj_dbh15_bivar_model3, bthuj_dbh15_bivar_model5)


varcomps(bthuj_dbh15_bivar_model1)
asreml::vpredict(bthuj_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.3897746 +/-0.2047885
asreml::vpredict(bthuj_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.03552683 +/-0.02765822

# sabinene - ht15
sab_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

sab_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                residual = ~ id(units):us(trait),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

sab_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

sab_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

sab_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_ht15_bivar_model2, sab_ht15_bivar_model1)
REMLRT.asreml(sab_ht15_bivar_model3, sab_ht15_bivar_model1)
REMLRT.asreml(sab_ht15_bivar_model2, sab_ht15_bivar_model3)
REMLRT.asreml(sab_ht15_bivar_model3, sab_ht15_bivar_model4)
REMLRT.asreml(sab_ht15_bivar_model3, sab_ht15_bivar_model5)


varcomps(sab_ht15_bivar_model1)
asreml::vpredict(sab_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # -0.07585974 +/-0.1516457
asreml::vpredict(sab_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.01085857 +/-0.02846426

# sabinene - dbh15
sab_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ id(units):us(trait),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_dbh15_bivar_model2, sab_dbh15_bivar_model1)
REMLRT.asreml(sab_dbh15_bivar_model3, sab_dbh15_bivar_model1)
REMLRT.asreml(sab_dbh15_bivar_model2, sab_dbh15_bivar_model3)
REMLRT.asreml(sab_dbh15_bivar_model3, sab_dbh15_bivar_model4)
REMLRT.asreml(sab_dbh15_bivar_model3, sab_dbh15_bivar_model5)


varcomps(sab_dbh15_bivar_model3)
asreml::vpredict(sab_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.02982608 +/-0.1698015
asreml::vpredict(sab_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.0328575 +/-0.0282203

# athujaplicin - bthujaplicin
alpha_thujaplicin_beta_thujaplicin_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicin_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicin_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06, #G.param = alpha_thujaplicin_beta_thujaplicin_bivar_model3$G.param, R.param = alpha_thujaplicin_beta_thujaplicin_bivar_model3$R.param)

alpha_thujaplicin_beta_thujaplicin_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicin_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model2, alpha_thujaplicin_beta_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model3, alpha_thujaplicin_beta_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model2, alpha_thujaplicin_beta_thujaplicin_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model3, alpha_thujaplicin_beta_thujaplicin_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model3, alpha_thujaplicin_beta_thujaplicin_bivar_model5)


varcomps(alpha_thujaplicin_beta_thujaplicin_bivar_model3)
asreml::vpredict(alpha_thujaplicin_beta_thujaplicin_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.4598395 +/-0.09668725
asreml::vpredict(alpha_thujaplicin_beta_thujaplicin_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.5384901 +/-0.02064679

# athujaplicin - gthujaplicin
alpha_thujaplicin_gamma_thujaplicin_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                          residual = ~ id(units):us(trait),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model2, alpha_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, alpha_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model2, alpha_thujaplicin_gamma_thujaplicin_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, alpha_thujaplicin_gamma_thujaplicin_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, alpha_thujaplicin_gamma_thujaplicin_bivar_model5)


varcomps(alpha_thujaplicin_gamma_thujaplicin_bivar_model3)
asreml::vpredict(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.5286983 +/-0.1001498
asreml::vpredict(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.5435931 +/-0.02005425

# athujaplicin - bthujaplicinol
alpha_thujaplicin_beta_thujaplicinol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model1$G.param, R.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model1$R.param)

alpha_thujaplicin_beta_thujaplicinol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model3$G.param, R.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model3$R.param)

alpha_thujaplicin_beta_thujaplicinol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model3$G.param, R.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model3$R.param)

alpha_thujaplicin_beta_thujaplicinol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicinol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model2, alpha_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, alpha_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model2, alpha_thujaplicin_beta_thujaplicinol_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, alpha_thujaplicin_beta_thujaplicinol_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, alpha_thujaplicin_beta_thujaplicinol_bivar_model5)

varcomps(alpha_thujaplicin_beta_thujaplicinol_bivar_model3)
asreml::vpredict(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.3540227 +/-0.1353299
asreml::vpredict(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.4663703 +/-0.02141168

# bthujaplicin - gthujaplicin
beta_thujaplicin_gamma_thujaplicin_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = beta_thujaplicin_gamma_thujaplicin_bivar_model1$G.param, R.param = beta_thujaplicin_gamma_thujaplicin_bivar_model1$R.param)

beta_thujaplicin_gamma_thujaplicin_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model2, beta_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model3, beta_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model2, beta_thujaplicin_gamma_thujaplicin_bivar_model3)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model3, beta_thujaplicin_gamma_thujaplicin_bivar_model4)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model3, beta_thujaplicin_gamma_thujaplicin_bivar_model5)


varcomps(beta_thujaplicin_gamma_thujaplicin_bivar_model3)
asreml::vpredict(beta_thujaplicin_gamma_thujaplicin_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2210959 +/-0.1197793
asreml::vpredict(beta_thujaplicin_gamma_thujaplicin_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.4071734 +/-0.02379809

# bthujaplicin - bthujaplicinol
beta_thujaplicin_beta_thujaplicinol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = beta_thujaplicin_beta_thujaplicinol_bivar_model1$G.param, R.param = beta_thujaplicin_beta_thujaplicinol_bivar_model1$R.param)

beta_thujaplicin_beta_thujaplicinol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model2, beta_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model3, beta_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model2, beta_thujaplicin_beta_thujaplicinol_bivar_model3)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model3, beta_thujaplicin_beta_thujaplicinol_bivar_model4)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model3, beta_thujaplicin_beta_thujaplicinol_bivar_model5)


varcomps(beta_thujaplicin_beta_thujaplicinol_bivar_model3)
asreml::vpredict(beta_thujaplicin_beta_thujaplicinol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1982669 +/-0.1395766
asreml::vpredict(beta_thujaplicin_beta_thujaplicinol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.3572472 +/-0.02393615

# gthujaplicin - bthujaplicinol
gamma_thujaplicin_beta_thujaplicinol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model1$G.param,
                                                           R.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model1$R.param)

gamma_thujaplicin_beta_thujaplicinol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

gamma_thujaplicin_beta_thujaplicinol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model3$G.param,
                                                           R.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model3$R.param)

gamma_thujaplicin_beta_thujaplicinol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

gamma_thujaplicin_beta_thujaplicinol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model2, gamma_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gamma_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model2, gamma_thujaplicin_beta_thujaplicinol_bivar_model3)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gamma_thujaplicin_beta_thujaplicinol_bivar_model4)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gamma_thujaplicin_beta_thujaplicinol_bivar_model5)

varcomps(gamma_thujaplicin_beta_thujaplicinol_bivar_model1)
asreml::vpredict(gamma_thujaplicin_beta_thujaplicinol_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.7642685 +/-0.06775346
asreml::vpredict(gamma_thujaplicin_beta_thujaplicinol_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.7596794 +/-0.012039

# athujaplicin - ht15
alpha_thujaplicin_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model2, alpha_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model3, alpha_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model2, alpha_thujaplicin_ht15_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model3, alpha_thujaplicin_ht15_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model3, alpha_thujaplicin_ht15_bivar_model5)

varcomps(alpha_thujaplicin_ht15_bivar_model1)
asreml::vpredict(alpha_thujaplicin_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.08181686 +/-0.1592407
asreml::vpredict(alpha_thujaplicin_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.03236499 +/-0.02885599
 

# athujaplicin - dbh15
alpha_thujaplicin_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ id(units):us(trait),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model2, alpha_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model3, alpha_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model2, alpha_thujaplicin_dbh15_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model3, alpha_thujaplicin_dbh15_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model3, alpha_thujaplicin_dbh15_bivar_model5)

varcomps(alpha_thujaplicin_dbh15_bivar_model3)
asreml::vpredict(alpha_thujaplicin_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.05101451 +/-0.1798843
asreml::vpredict(alpha_thujaplicin_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.101446 +/-0.02793766

# bthujaplicin - ht15
beta_thujaplicin_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06, G.param = beta_thujaplicin_ht15_bivar_model1$G.param, 
                                              R.param = beta_thujaplicin_ht15_bivar_model1$R.param)

beta_thujaplicin_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ id(units):us(trait),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06, G.param = beta_thujaplicin_ht15_bivar_model3$G.param, 
                                             R.param = beta_thujaplicin_ht15_bivar_model3$R.param)

beta_thujaplicin_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_ht15_bivar_model2, beta_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model3, beta_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model2, beta_thujaplicin_ht15_bivar_model3)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model3, beta_thujaplicin_ht15_bivar_model4)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model3, beta_thujaplicin_ht15_bivar_model5)

varcomps(beta_thujaplicin_ht15_bivar_model1)
asreml::vpredict(beta_thujaplicin_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.02665735 +/-0.15001
asreml::vpredict(beta_thujaplicin_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.1084534 +/-0.02829591

# bthujaplicin - dbh15
beta_thujaplicin_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ id(units):us(trait),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model2, beta_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model3, beta_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model2, beta_thujaplicin_dbh15_bivar_model3)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model3, beta_thujaplicin_dbh15_bivar_model4)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model3, beta_thujaplicin_dbh15_bivar_model5)


varcomps(beta_thujaplicin_dbh15_bivar_model3)
asreml::vpredict(beta_thujaplicin_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1186306 +/- 0.1680799
asreml::vpredict(beta_thujaplicin_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.1474809 +/-0.02767835

# gthujaplicin - ht15
gamma_thujaplicin_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06, G.param = gamma_thujaplicin_ht15_bivar_model1$G.param, 
                                             R.param = gamma_thujaplicin_ht15_bivar_model1$R.param)

gamma_thujaplicin_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ id(units):us(trait),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06, G.param = gamma_thujaplicin_ht15_bivar_model3$G.param, 
                                             R.param = gamma_thujaplicin_ht15_bivar_model3$R.param)

gamma_thujaplicin_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model2, gamma_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model3, gamma_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model2, gamma_thujaplicin_ht15_bivar_model3)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model3, gamma_thujaplicin_ht15_bivar_model4)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model3, gamma_thujaplicin_ht15_bivar_model5)

varcomps(gamma_thujaplicin_ht15_bivar_model1)
asreml::vpredict(gamma_thujaplicin_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.03026305 +/-0.161835
asreml::vpredict(gamma_thujaplicin_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.01114055 +/-0.02927066

# gthujaplicin - dbh15
gamma_thujaplicin_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_dbh15_bivar_model1$G.param, R.param = gamma_thujaplicin_dbh15_bivar_model1$R.param)

gamma_thujaplicin_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ id(units):us(trait),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_dbh15_bivar_model3$G.param, R.param = gamma_thujaplicin_dbh15_bivar_model3$R.param)

gamma_thujaplicin_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model2, gamma_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model3, gamma_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model2, gamma_thujaplicin_dbh15_bivar_model3)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model3, gamma_thujaplicin_dbh15_bivar_model4)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model3, gamma_thujaplicin_dbh15_bivar_model5)

varcomps(gamma_thujaplicin_dbh15_bivar_model3)
asreml::vpredict(gamma_thujaplicin_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.02798345 +/-0.1821306
asreml::vpredict(gamma_thujaplicin_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.06682985 +/-0.02858625

# bthujaplicinol - ht15
beta_thujaplicinol_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06, G.param = beta_thujaplicinol_ht15_bivar_model1$G.param, 
                                             R.param = beta_thujaplicinol_ht15_bivar_model1$R.param)

beta_thujaplicinol_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ id(units):us(trait),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

beta_thujaplicinol_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicinol_ht15_bivar_model3$G.param, R.param = beta_thujaplicinol_ht15_bivar_model3$R.param)

beta_thujaplicinol_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

beta_thujaplicinol_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model2, beta_thujaplicinol_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model3, beta_thujaplicinol_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model2, beta_thujaplicinol_ht15_bivar_model3)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model3, beta_thujaplicinol_ht15_bivar_model4)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model3, beta_thujaplicinol_ht15_bivar_model5)

varcomps(beta_thujaplicinol_ht15_bivar_model1)
asreml::vpredict(beta_thujaplicinol_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # -0.1040754 +/-0.1798054
asreml::vpredict(beta_thujaplicinol_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.1233018 +/-0.02908282

# bthujaplicinol - dbh15
beta_thujaplicinol_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06, G.param = beta_thujaplicinol_dbh15_bivar_model1$G.param, 
R.param = beta_thujaplicinol_dbh15_bivar_model1$R.param)

beta_thujaplicinol_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                               residual = ~ id(units):us(trait),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)

beta_thujaplicinol_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicinol_dbh15_bivar_model3$G.param, R.param = beta_thujaplicinol_dbh15_bivar_model3$R.param)

beta_thujaplicinol_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)

beta_thujaplicinol_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model2, beta_thujaplicinol_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model3, beta_thujaplicinol_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model2, beta_thujaplicinol_dbh15_bivar_model3)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model3, beta_thujaplicinol_dbh15_bivar_model4)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model3, beta_thujaplicinol_dbh15_bivar_model5)

varcomps(beta_thujaplicinol_dbh15_bivar_model1)
asreml::vpredict(beta_thujaplicinol_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # -0.06887698 +/-0.2022927
asreml::vpredict(beta_thujaplicinol_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.1787882 +/-0.02814062

# ht15 - dbh15
ht15_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)#, G.param = ht15_dbh15_bivar_model1$G.param, R.param = ht15_dbh15_bivar_model1$R.param)

ht15_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                residual = ~ id(units):us(trait),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)

ht15_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)#, G.param = ht15_dbh15_bivar_model3$G.param, R.param = ht15_dbh15_bivar_model3$R.param)

ht15_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.giv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)

ht15_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.giv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)

REMLRT.asreml(ht15_dbh15_bivar_model2, ht15_dbh15_bivar_model1)
REMLRT.asreml(ht15_dbh15_bivar_model3, ht15_dbh15_bivar_model1)
REMLRT.asreml(ht15_dbh15_bivar_model2, ht15_dbh15_bivar_model3)
REMLRT.asreml(ht15_dbh15_bivar_model3, ht15_dbh15_bivar_model4)
REMLRT.asreml(ht15_dbh15_bivar_model3, ht15_dbh15_bivar_model5)

varcomps(ht15_dbh15_bivar_model3)
asreml::vpredict(ht15_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.7490649 +/-0.08400558
asreml::vpredict(ht15_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.8476982 +/-0.008267215
