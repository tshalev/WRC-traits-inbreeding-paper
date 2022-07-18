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

# Some plots
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

ggsave("monoterpene_density_plot.svg", dpi = 300, height = 8, width = 10)

# Pedigree
wrc.pheno.ped <- read.csv("wrc_training_ped.csv", header = F, na.strings = "NA") 
head(wrc.pheno.ped)

wrc.ainv <- ainverse(wrc.pheno.ped)

matching_ped_dat <- as.integer(setdiff(wrc.pheno.dat$tree, wrc.pheno.ped$V1))
wrc.pheno.dat <- wrc.pheno.dat %>% 
  filter(!(tree %in% matching_ped_dat))

cor(wrc.pheno.dat$fath, wrc.pheno.dat$fmono, use = "na.or.complete")

plot(fath ~ fbth, wrc.pheno.dat)

# Read sample ids
indivs <- read.csv("training_samples_genome_paper_snps_maf01.012.indv", header = F, sep = "\t")

matchingSamples <- intersect(wrc.pheno.dat$tree, indivs$V1)


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
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_bthuj_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, fbth) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_bthuj_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, fbth) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_bthuj_bivar_model2, athuj_bthuj_bivar_model1)
REMLRT.asreml(athuj_bthuj_bivar_model3, athuj_bthuj_bivar_model1)
REMLRT.asreml(athuj_bthuj_bivar_model2, athuj_bthuj_bivar_model3)


varcomps(athuj_bthuj_bivar_model1)
asreml::vpredict(athuj_bthuj_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.6489049 +/-0.09785095
asreml::vpredict(athuj_bthuj_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.6035576 +/-0.0186063


# a-thujone - sabinene
athuj_sab_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ id(units):us(trait),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

athuj_sab_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fath, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_sab_bivar_model2, athuj_sab_bivar_model1)
REMLRT.asreml(athuj_sab_bivar_model3, athuj_sab_bivar_model1)
REMLRT.asreml(athuj_sab_bivar_model2, athuj_sab_bivar_model3)
REMLRT.asreml(athuj_sab_bivar_model3, athuj_sab_bivar_model4)
REMLRT.asreml(athuj_sab_bivar_model3, athuj_sab_bivar_model5)


varcomps(athuj_sab_bivar_model3)
asreml::vpredict(athuj_sab_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.3545644 +/-0.1055315
asreml::vpredict(athuj_sab_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.6510925 +/-0.01717888

# b-thujone - sabinene
bthuj_sab_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ id(units):us(trait),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

bthuj_sab_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fbth, fsab) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_sab_bivar_model2, bthuj_sab_bivar_model1)
REMLRT.asreml(bthuj_sab_bivar_model3, bthuj_sab_bivar_model1)
REMLRT.asreml(bthuj_sab_bivar_model2, bthuj_sab_bivar_model3)
REMLRT.asreml(bthuj_sab_bivar_model3, bthuj_sab_bivar_model4)
REMLRT.asreml(bthuj_sab_bivar_model3, bthuj_sab_bivar_model5)


varcomps(bthuj_sab_bivar_model1)
asreml::vpredict(bthuj_sab_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.298448 +/-0.1359854
asreml::vpredict(bthuj_sab_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.3187039 +/-0.02795357

# a-thujone - a-thujaplicin
athuj_athuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_athuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_athuja_bivar_model2, athuj_athuja_bivar_model1)
REMLRT.asreml(athuj_athuja_bivar_model3, athuj_athuja_bivar_model1)
REMLRT.asreml(athuj_athuja_bivar_model2, athuj_athuja_bivar_model3)
REMLRT.asreml(athuj_athuja_bivar_model3, athuj_athuja_bivar_model4)
REMLRT.asreml(athuj_athuja_bivar_model3, athuj_athuja_bivar_model5)


varcomps(athuj_athuja_bivar_model3)
asreml::vpredict(athuj_athuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1625098 +/-0.1299578
asreml::vpredict(athuj_athuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3)))

# a-thujone - b-thujaplicin
athuj_bthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_bthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_bthuja_bivar_model2, athuj_bthuja_bivar_model1)
REMLRT.asreml(athuj_bthuja_bivar_model3, athuj_bthuja_bivar_model1)
REMLRT.asreml(athuj_bthuja_bivar_model2, athuj_bthuja_bivar_model3)
REMLRT.asreml(athuj_bthuja_bivar_model3, athuj_bthuja_bivar_model4)
REMLRT.asreml(athuj_bthuja_bivar_model3, athuj_bthuja_bivar_model5)


varcomps(athuj_bthuja_bivar_model3)
asreml::vpredict(athuj_bthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1625098 +/-0.1299578
asreml::vpredict(athuj_bthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.01611803 +/-0.02830258

# a-thujone - g-thujaplicin
athuj_gthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

athuj_gthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fath, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_gthuja_bivar_model2, athuj_gthuja_bivar_model1)
REMLRT.asreml(athuj_gthuja_bivar_model3, athuj_gthuja_bivar_model1)
REMLRT.asreml(athuj_gthuja_bivar_model2, athuj_gthuja_bivar_model3)
REMLRT.asreml(athuj_gthuja_bivar_model3, athuj_gthuja_bivar_model4)
REMLRT.asreml(athuj_gthuja_bivar_model3, athuj_gthuja_bivar_model5)


varcomps(athuj_gthuja_bivar_model3)
asreml::vpredict(athuj_gthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.08088585 +/-0.1377588
asreml::vpredict(athuj_gthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.03085237 +/-0.02807806

# a-thujone - b-thujaplicinol
athuj_bthuja_ol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ id(units):us(trait),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

athuj_bthuja_ol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fath, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_bthuja_ol_bivar_model2, athuj_bthuja_ol_bivar_model1)
REMLRT.asreml(athuj_bthuja_ol_bivar_model3, athuj_bthuja_ol_bivar_model1)
REMLRT.asreml(athuj_bthuja_ol_bivar_model2, athuj_bthuja_ol_bivar_model3)
REMLRT.asreml(athuj_bthuja_ol_bivar_model3, athuj_bthuja_ol_bivar_model4)
REMLRT.asreml(athuj_bthuja_ol_bivar_model3, athuj_bthuja_ol_bivar_model5)


varcomps(athuj_bthuja_ol_bivar_model3)
asreml::vpredict(athuj_bthuja_ol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.1769193 +/-0.1582504 
asreml::vpredict(athuj_bthuja_ol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.06149127 +/-0.02738515

# b-thujone - a-thujaplicin
bthuj_athuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_athuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, alpha_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_athuja_bivar_model2, bthuj_athuja_bivar_model1)
REMLRT.asreml(bthuj_athuja_bivar_model3, bthuj_athuja_bivar_model1)
REMLRT.asreml(bthuj_athuja_bivar_model2, bthuj_athuja_bivar_model3)
REMLRT.asreml(bthuj_athuja_bivar_model3, bthuj_athuja_bivar_model4)
REMLRT.asreml(bthuj_athuja_bivar_model3, bthuj_athuja_bivar_model5)


varcomps(bthuj_athuja_bivar_model3)
asreml::vpredict(bthuj_athuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1625098 +/-0.1299578
asreml::vpredict(bthuj_athuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3)))

# b-thujone - b-thujaplicin
bthuj_bthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
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
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_bthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_bthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, beta_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_bthuja_bivar_model2, bthuj_bthuja_bivar_model1)
REMLRT.asreml(bthuj_bthuja_bivar_model3, bthuj_bthuja_bivar_model1)
REMLRT.asreml(bthuj_bthuja_bivar_model2, bthuj_bthuja_bivar_model3)
REMLRT.asreml(bthuj_bthuja_bivar_model3, bthuj_bthuja_bivar_model4)
REMLRT.asreml(bthuj_bthuja_bivar_model3, bthuj_bthuja_bivar_model5)


varcomps(bthuj_bthuja_bivar_model1)
asreml::vpredict(bthuj_bthuja_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(bthuj_bthuja_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714

varcomps(bthuj_bthuja_bivar_model3)
asreml::vpredict(bthuj_bthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.1625098 +/-0.1299578
asreml::vpredict(bthuj_bthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3)))

# b-thujone - g-thujaplicin
bthuj_gthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06, G.param = bthuj_gthuja_bivar_model1$G.param, R.param = bthuj_gthuja_bivar_model1$R.param)

bthuj_gthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ id(units):us(trait),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_gthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_gthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

bthuj_gthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                    fixed = cbind(fbth, gamma_thujaplicin) ~ trait + trait:site,
                                    random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                    residual = ~ dsum(~ id(units):us(trait) | site),
                                    na.action = na.method(x = "include", y = "include"),
                                    maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_gthuja_bivar_model2, bthuj_gthuja_bivar_model1)
REMLRT.asreml(bthuj_gthuja_bivar_model3, bthuj_gthuja_bivar_model1)
REMLRT.asreml(bthuj_gthuja_bivar_model2, bthuj_gthuja_bivar_model3)
REMLRT.asreml(bthuj_gthuja_bivar_model3, bthuj_gthuja_bivar_model4)
REMLRT.asreml(bthuj_gthuja_bivar_model3, bthuj_gthuja_bivar_model5)

varcomps(bthuj_gthuja_bivar_model1)
asreml::vpredict(bthuj_gthuja_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(bthuj_gthuja_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714

varcomps(bthuj_gthuja_bivar_model3)
asreml::vpredict(bthuj_gthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.01219317 +/-0.1720537
asreml::vpredict(bthuj_gthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.00477361 +/-0.02653703

# b-thujone - b-thujaplicinol
bthuj_bthuja_ol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ id(units):us(trait),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

bthuj_bthuja_ol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                       fixed = cbind(fbth, beta_thujaplicinol) ~ trait + trait:site,
                                       random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                       residual = ~ dsum(~ id(units):us(trait) | site),
                                       na.action = na.method(x = "include", y = "include"),
                                       maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_bthuja_ol_bivar_model2, bthuj_bthuja_ol_bivar_model1)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model3, bthuj_bthuja_ol_bivar_model1)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model2, bthuj_bthuja_ol_bivar_model3)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model3, bthuj_bthuja_ol_bivar_model4)
REMLRT.asreml(bthuj_bthuja_ol_bivar_model3, bthuj_bthuja_ol_bivar_model5)


varcomps(bthuj_bthuja_ol_bivar_model1)
asreml::vpredict(bthuj_bthuja_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(bthuj_bthuja_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714


varcomps(bthuj_bthuja_ol_bivar_model3)
asreml::vpredict(bthuj_bthuja_ol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # -0.1535861 +/-0.1912298
asreml::vpredict(bthuj_bthuja_ol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.02725478 +/-0.02533936


# sabinene - a-thujaplicin
sab_athuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_athuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, alpha_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_athuja_bivar_model2, sab_athuja_bivar_model1)
REMLRT.asreml(sab_athuja_bivar_model3, sab_athuja_bivar_model1)
REMLRT.asreml(sab_athuja_bivar_model2, sab_athuja_bivar_model3)
REMLRT.asreml(sab_athuja_bivar_model3, sab_athuja_bivar_model4)
REMLRT.asreml(sab_athuja_bivar_model3, sab_athuja_bivar_model5)


varcomps(sab_athuja_bivar_model1)
asreml::vpredict(sab_athuja_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(sab_athuja_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714

varcomps(sab_athuja_bivar_model3)
asreml::vpredict(sab_athuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2900154 +/-0.1329092
asreml::vpredict(sab_athuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.02575366 +/-0.02935965

# sabinene - b-thujaplicin
sab_bthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_bthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, beta_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_bthuja_bivar_model2, sab_bthuja_bivar_model1)
REMLRT.asreml(sab_bthuja_bivar_model3, sab_bthuja_bivar_model1)
REMLRT.asreml(sab_bthuja_bivar_model2, sab_bthuja_bivar_model3)
REMLRT.asreml(sab_bthuja_bivar_model3, sab_bthuja_bivar_model4)
REMLRT.asreml(sab_bthuja_bivar_model3, sab_bthuja_bivar_model5)


varcomps(sab_bthuja_bivar_model1)
asreml::vpredict(sab_bthuja_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(sab_bthuja_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714


varcomps(sab_bthuja_bivar_model3)
asreml::vpredict(sab_bthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.04055208 +/-0.1227806
asreml::vpredict(sab_bthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.009094749 +/-0.02922578

# sabinene - g-thujaplicin
sab_gthuja_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

sab_gthuja_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fsab, gamma_thujaplicin) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_gthuja_bivar_model2, sab_gthuja_bivar_model1)
REMLRT.asreml(sab_gthuja_bivar_model3, sab_gthuja_bivar_model1)
REMLRT.asreml(sab_gthuja_bivar_model2, sab_gthuja_bivar_model3)
REMLRT.asreml(sab_gthuja_bivar_model3, sab_gthuja_bivar_model4)
REMLRT.asreml(sab_gthuja_bivar_model3, sab_gthuja_bivar_model5)


varcomps(sab_gthuja_bivar_model1)
asreml::vpredict(sab_gthuja_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(sab_gthuja_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714

varcomps(sab_gthuja_bivar_model3)
asreml::vpredict(sab_gthuja_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.08226586 +/-0.1319484
asreml::vpredict(sab_gthuja_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.07781355 +/-0.03328632

# sabinene - b-thujaplicinol
sab_bthuja_ol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                     residual = ~ id(units):us(trait),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

sab_bthuja_ol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                     fixed = cbind(fsab, beta_thujaplicinol) ~ trait + trait:site,
                                     random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                     residual = ~ dsum(~ id(units):us(trait) | site),
                                     na.action = na.method(x = "include", y = "include"),
                                     maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_bthuja_ol_bivar_model2, sab_bthuja_ol_bivar_model1)
REMLRT.asreml(sab_bthuja_ol_bivar_model3, sab_bthuja_ol_bivar_model1)
REMLRT.asreml(sab_bthuja_ol_bivar_model2, sab_bthuja_ol_bivar_model3)
REMLRT.asreml(sab_bthuja_ol_bivar_model3, sab_bthuja_ol_bivar_model4)
REMLRT.asreml(sab_bthuja_ol_bivar_model5, sab_bthuja_ol_bivar_model1)


varcomps(sab_bthuja_ol_bivar_model5)
asreml::vpredict(sab_bthuja_ol_bivar_model5, gc ~ V6 / sqrt(V5 * V7)) # 0.00711989 +/-0.150459
asreml::vpredict(sab_bthuja_ol_bivar_model5, pc ~ (V6 + (V10 + V14 + V18)/3)/ sqrt((V5 +(V9 + V13 + V17)/3) * (V7 + (V11 + V15 + V19)/3))) # 0.00259667 +/-0.002057941

# a-thujone - ht15
athuj_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

athuj_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fath, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_ht15_bivar_model2, athuj_ht15_bivar_model1)
REMLRT.asreml(athuj_ht15_bivar_model3, athuj_ht15_bivar_model1)
REMLRT.asreml(athuj_ht15_bivar_model2, athuj_ht15_bivar_model3)
REMLRT.asreml(athuj_ht15_bivar_model3, athuj_ht15_bivar_model4)
REMLRT.asreml(athuj_ht15_bivar_model3, athuj_ht15_bivar_model5)


varcomps(athuj_ht15_bivar_model1)
asreml::vpredict(athuj_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.04847007 +/-0.1591267
asreml::vpredict(athuj_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.05448618 +/-0.02936389

varcomps(athuj_ht15_bivar_model3)
asreml::vpredict(athuj_ht15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.08226586 +/-0.1319484
asreml::vpredict(athuj_ht15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # -0.07781355 +/-0.03328632

# a-thujone - dbh15
athuj_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

athuj_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fath, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

REMLRT.asreml(athuj_dbh15_bivar_model2, athuj_dbh15_bivar_model1)
REMLRT.asreml(athuj_dbh15_bivar_model3, athuj_dbh15_bivar_model1)
REMLRT.asreml(athuj_dbh15_bivar_model2, athuj_dbh15_bivar_model3)
REMLRT.asreml(athuj_dbh15_bivar_model3, athuj_dbh15_bivar_model4)
REMLRT.asreml(athuj_dbh15_bivar_model3, athuj_dbh15_bivar_model5)


varcomps(athuj_dbh15_bivar_model3)
asreml::vpredict(athuj_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2119233 +/-0.1734458
asreml::vpredict(athuj_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.09892769 +/-0.0286578

# b-thujone - ht15
bthuj_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ id(units):us(trait),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

bthuj_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                  fixed = cbind(fbth, ht15) ~ trait + trait:site,
                                  random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                  residual = ~ dsum(~ id(units):us(trait) | site),
                                  na.action = na.method(x = "include", y = "include"),
                                  maxit = 100, workspace = 80e+06)

REMLRT.asreml(bthuj_ht15_bivar_model2, bthuj_ht15_bivar_model1)
REMLRT.asreml(bthuj_ht15_bivar_model3, bthuj_ht15_bivar_model1)
REMLRT.asreml(bthuj_ht15_bivar_model2, bthuj_ht15_bivar_model3)
REMLRT.asreml(bthuj_ht15_bivar_model3, bthuj_ht15_bivar_model4)
REMLRT.asreml(bthuj_ht15_bivar_model3, bthuj_ht15_bivar_model5)


varcomps(bthuj_ht15_bivar_model1)
asreml::vpredict(bthuj_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2166937 +/-0.1888924
asreml::vpredict(bthuj_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.04135103 +/-0.02762858

# b-thujone - dbh15
bthuj_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06, G.param = bthuj_dbh15_bivar_model1$G.param, R.param = bthuj_dbh15_bivar_model1$R.param)

bthuj_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

bthuj_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06, G.param = bthuj_dbh15_bivar_model3$G.param, R.param = bthuj_dbh15_bivar_model3$R.param)

bthuj_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

bthuj_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(fbth, dbh15) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)



REMLRT.asreml(bthuj_dbh15_bivar_model2, bthuj_dbh15_bivar_model1)
REMLRT.asreml(bthuj_dbh15_bivar_model3, bthuj_dbh15_bivar_model1)
REMLRT.asreml(bthuj_dbh15_bivar_model2, bthuj_dbh15_bivar_model3)
REMLRT.asreml(bthuj_dbh15_bivar_model3, bthuj_dbh15_bivar_model4)
REMLRT.asreml(bthuj_dbh15_bivar_model3, bthuj_dbh15_bivar_model5)


varcomps(bthuj_dbh15_bivar_model1)
asreml::vpredict(bthuj_dbh15_bivar_model2, gc ~ V8 / sqrt(V7 * V9)) # 0.3397626 +/-0.2140569
asreml::vpredict(bthuj_dbh15_bivar_model2, pc ~ (V8 + V12)/ sqrt((V7 + V11) * (V9 + V13))) # -0.03983509 +/-0.02863514

varcomps(bthuj_dbh15_bivar_model3)
asreml::vpredict(bthuj_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2119233 +/-0.1734458
asreml::vpredict(bthuj_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.09892769 +/-0.0286578

# sabinene - ht15
sab_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06, G.param = sab_ht15_bivar_model1$G.param, R.param = sab_ht15_bivar_model1$R.param)

sab_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                residual = ~ id(units):us(trait),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

sab_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06, G.param = sab_ht15_bivar_model3$G.param, R.param = sab_ht15_bivar_model3$R.param)

sab_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

sab_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                fixed = cbind(fsab, ht15) ~ trait + trait:site,
                                random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                residual = ~ dsum(~ id(units):us(trait) | site),
                                na.action = na.method(x = "include", y = "include"),
                                maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_ht15_bivar_model2, sab_ht15_bivar_model1)
REMLRT.asreml(sab_ht15_bivar_model3, sab_ht15_bivar_model1)
REMLRT.asreml(sab_ht15_bivar_model2, sab_ht15_bivar_model3)
REMLRT.asreml(sab_ht15_bivar_model3, sab_ht15_bivar_model4)
REMLRT.asreml(sab_ht15_bivar_model3, sab_ht15_bivar_model5)


varcomps(sab_ht15_bivar_model1)
asreml::vpredict(sab_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # -0.05827582 +/-0.1511344
asreml::vpredict(sab_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # -0.01184089 +/-0.02878424

# sabinene - dbh15
sab_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ id(units):us(trait),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06, G.param = sab_dbh15_bivar_model3$G.param, R.param = sab_dbh15_bivar_model3$R.param)

sab_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

sab_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                 fixed = cbind(fsab, dbh15) ~ trait + trait:site,
                                 random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                 residual = ~ dsum(~ id(units):us(trait) | site),
                                 na.action = na.method(x = "include", y = "include"),
                                 maxit = 100, workspace = 80e+06)

REMLRT.asreml(sab_dbh15_bivar_model2, sab_dbh15_bivar_model1)
REMLRT.asreml(sab_dbh15_bivar_model3, sab_dbh15_bivar_model1)
REMLRT.asreml(sab_dbh15_bivar_model2, sab_dbh15_bivar_model3)
REMLRT.asreml(sab_dbh15_bivar_model3, sab_dbh15_bivar_model4)
REMLRT.asreml(sab_dbh15_bivar_model3, sab_dbh15_bivar_model5)


varcomps(sab_dbh15_bivar_model1)
asreml::vpredict(sab_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # -0.01173373 +/-0.1697236
asreml::vpredict(sab_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.02982982 +/-0.0284169

varcomps(sab_dbh15_bivar_model3)
asreml::vpredict(sab_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2119233 +/-0.1734458
asreml::vpredict(sab_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.09892769 +/-0.0286578

# athujaplicin - bthujaplicin
alpha_thujaplicin_beta_thujaplicin_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicin_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ id(units):us(trait),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicin_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)#, G.param = alpha_thujaplicin_beta_thujaplicin_bivar_model3$G.param, R.param = alpha_thujaplicin_beta_thujaplicin_bivar_model3$R.param)

alpha_thujaplicin_beta_thujaplicin_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicin_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                   fixed = cbind(alpha_thujaplicin, beta_thujaplicin) ~ trait + trait:site,
                                   random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                   residual = ~ dsum(~ id(units):us(trait) | site),
                                   na.action = na.method(x = "include", y = "include"),
                                   maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model2, alpha_thujaplicin_beta_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model3, alpha_thujaplicin_beta_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model2, alpha_thujaplicin_beta_thujaplicin_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model3, alpha_thujaplicin_beta_thujaplicin_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicin_bivar_model3, alpha_thujaplicin_beta_thujaplicin_bivar_model5)

varcomps(alpha_thujaplicin_beta_thujaplicin_bivar_model1)
asreml::vpredict(alpha_thujaplicin_beta_thujaplicin_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # -0.01173373 +/-0.1697236
asreml::vpredict(alpha_thujaplicin_beta_thujaplicin_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.02982982 +/-0.0284169

varcomps(alpha_thujaplicin_beta_thujaplicin_bivar_model3)
asreml::vpredict(alpha_thujaplicin_beta_thujaplicin_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_beta_thujaplicin_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# athujaplicin - gthujaplicin
alpha_thujaplicin_gamma_thujaplicin_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                          residual = ~ id(units):us(trait),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

alpha_thujaplicin_gamma_thujaplicin_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                          fixed = cbind(alpha_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                          random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                          residual = ~ dsum(~ id(units):us(trait) | site),
                                                          na.action = na.method(x = "include", y = "include"),
                                                          maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model2, alpha_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, alpha_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model2, alpha_thujaplicin_gamma_thujaplicin_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, alpha_thujaplicin_gamma_thujaplicin_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, alpha_thujaplicin_gamma_thujaplicin_bivar_model5)


varcomps(alpha_thujaplicin_gamma_thujaplicin_bivar_model3)
asreml::vpredict(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_gamma_thujaplicin_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# athujaplicin - bthujaplicinol
alpha_thujaplicin_beta_thujaplicinol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model1$G.param, R.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model1$R.param)

alpha_thujaplicin_beta_thujaplicinol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicinol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model3$G.param, R.param = alpha_thujaplicin_beta_thujaplicinol_bivar_model3$R.param)

alpha_thujaplicin_beta_thujaplicinol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_beta_thujaplicinol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model2, alpha_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, alpha_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model2, alpha_thujaplicin_beta_thujaplicinol_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, alpha_thujaplicin_beta_thujaplicinol_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, alpha_thujaplicin_beta_thujaplicinol_bivar_model5)

varcomps(alpha_thujaplicin_beta_thujaplicinol_bivar_model1)
asreml::vpredict(alpha_thujaplicin_beta_thujaplicinol_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_beta_thujaplicinol_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(alpha_thujaplicin_beta_thujaplicinol_bivar_model3)
asreml::vpredict(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_beta_thujaplicinol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# bthujaplicin - gthujaplicin
beta_thujaplicin_gamma_thujaplicin_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_gamma_thujaplicin_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, gamma_thujaplicin) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model2, beta_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model3, beta_thujaplicin_gamma_thujaplicin_bivar_model1)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model2, beta_thujaplicin_gamma_thujaplicin_bivar_model3)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model3, beta_thujaplicin_gamma_thujaplicin_bivar_model4)
REMLRT.asreml(beta_thujaplicin_gamma_thujaplicin_bivar_model3, beta_thujaplicin_gamma_thujaplicin_bivar_model5)


varcomps(beta_thujaplicin_gamma_thujaplicin_bivar_model3)
asreml::vpredict(beta_thujaplicin_gamma_thujaplicin_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicin_gamma_thujaplicin_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# bthujaplicin - bthujaplicinol
beta_thujaplicin_beta_thujaplicinol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

beta_thujaplicin_beta_thujaplicinol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(beta_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model2, beta_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model3, beta_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model2, beta_thujaplicin_beta_thujaplicinol_bivar_model3)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model3, beta_thujaplicin_beta_thujaplicinol_bivar_model4)
REMLRT.asreml(beta_thujaplicin_beta_thujaplicinol_bivar_model3, beta_thujaplicin_beta_thujaplicinol_bivar_model5)


varcomps(beta_thujaplicin_beta_thujaplicinol_bivar_model3)
asreml::vpredict(beta_thujaplicin_beta_thujaplicinol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicin_beta_thujaplicinol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# gthujaplicin - bthujaplicinol
gamma_thujaplicin_beta_thujaplicinol_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model1$G.param,
                                                           R.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model1$R.param)

gamma_thujaplicin_beta_thujaplicinol_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

gamma_thujaplicin_beta_thujaplicinol_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06, G.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model3$G.param,
                                                           R.param = gamma_thujaplicin_beta_thujaplicinol_bivar_model3$R.param)

gamma_thujaplicin_beta_thujaplicinol_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

gamma_thujaplicin_beta_thujaplicinol_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(gamma_thujaplicin, beta_thujaplicinol) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model2, gamma_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gamma_thujaplicin_beta_thujaplicinol_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model2, gamma_thujaplicin_beta_thujaplicinol_bivar_model3)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gamma_thujaplicin_beta_thujaplicinol_bivar_model4)
REMLRT.asreml(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gamma_thujaplicin_beta_thujaplicinol_bivar_model5)

varcomps(gamma_thujaplicin_beta_thujaplicinol_bivar_model1)
asreml::vpredict(gamma_thujaplicin_beta_thujaplicinol_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.07216421 +/-0.1622209
asreml::vpredict(gamma_thujaplicin_beta_thujaplicinol_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.01968432 +/-0.02781714

varcomps(gamma_thujaplicin_beta_thujaplicinol_bivar_model3)
asreml::vpredict(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(gamma_thujaplicin_beta_thujaplicinol_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# athujaplicin - ht15
alpha_thujaplicin_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ id(units):us(trait),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

alpha_thujaplicin_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                           fixed = cbind(alpha_thujaplicin, ht15) ~ trait + trait:site,
                                                           random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                           residual = ~ dsum(~ id(units):us(trait) | site),
                                                           na.action = na.method(x = "include", y = "include"),
                                                           maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model2, alpha_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model3, alpha_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model2, alpha_thujaplicin_ht15_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model3, alpha_thujaplicin_ht15_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_ht15_bivar_model3, alpha_thujaplicin_ht15_bivar_model5)

varcomps(alpha_thujaplicin_ht15_bivar_model1)
asreml::vpredict(alpha_thujaplicin_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(alpha_thujaplicin_ht15_bivar_model3)
asreml::vpredict(alpha_thujaplicin_ht15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_ht15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# athujaplicin - dbh15
alpha_thujaplicin_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ id(units):us(trait),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

alpha_thujaplicin_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(alpha_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model2, alpha_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model3, alpha_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model2, alpha_thujaplicin_dbh15_bivar_model3)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model3, alpha_thujaplicin_dbh15_bivar_model4)
REMLRT.asreml(alpha_thujaplicin_dbh15_bivar_model3, alpha_thujaplicin_dbh15_bivar_model5)

varcomps(alpha_thujaplicin_dbh15_bivar_model1)
asreml::vpredict(alpha_thujaplicin_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(alpha_thujaplicin_dbh15_bivar_model3)
asreml::vpredict(alpha_thujaplicin_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(alpha_thujaplicin_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# bthujaplicin - ht15
beta_thujaplicin_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicin_ht15_bivar_model1$G.param, 
                                              R.param = beta_thujaplicin_ht15_bivar_model1$R.param)

beta_thujaplicin_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ id(units):us(trait),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicin_ht15_bivar_model3$G.param, 
                                             R.param = beta_thujaplicin_ht15_bivar_model3$R.param)

beta_thujaplicin_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, ht15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_ht15_bivar_model2, beta_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model3, beta_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model2, beta_thujaplicin_ht15_bivar_model3)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model3, beta_thujaplicin_ht15_bivar_model4)
REMLRT.asreml(beta_thujaplicin_ht15_bivar_model3, beta_thujaplicin_ht15_bivar_model5)

varcomps(beta_thujaplicin_ht15_bivar_model1)
asreml::vpredict(beta_thujaplicin_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicin_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(beta_thujaplicin_ht15_bivar_model3)
asreml::vpredict(beta_thujaplicin_ht15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicin_ht15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# bthujaplicin - dbh15
beta_thujaplicin_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ id(units):us(trait),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

beta_thujaplicin_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                              fixed = cbind(beta_thujaplicin, dbh15) ~ trait + trait:site,
                                              random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                              residual = ~ dsum(~ id(units):us(trait) | site),
                                              na.action = na.method(x = "include", y = "include"),
                                              maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model2, beta_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model3, beta_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model2, beta_thujaplicin_dbh15_bivar_model3)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model3, beta_thujaplicin_dbh15_bivar_model4)
REMLRT.asreml(beta_thujaplicin_dbh15_bivar_model3, beta_thujaplicin_dbh15_bivar_model5)


varcomps(beta_thujaplicin_dbh15_bivar_model3)
asreml::vpredict(beta_thujaplicin_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicin_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# gthujaplicin - ht15
gamma_thujaplicin_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_ht15_bivar_model1$G.param, 
                                             R.param = gamma_thujaplicin_ht15_bivar_model1$R.param)

gamma_thujaplicin_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ id(units):us(trait),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_ht15_bivar_model3$G.param, 
                                             R.param = gamma_thujaplicin_ht15_bivar_model3$R.param)

gamma_thujaplicin_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model2, gamma_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model3, gamma_thujaplicin_ht15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model2, gamma_thujaplicin_ht15_bivar_model3)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model3, gamma_thujaplicin_ht15_bivar_model4)
REMLRT.asreml(gamma_thujaplicin_ht15_bivar_model3, gamma_thujaplicin_ht15_bivar_model5)

varcomps(gamma_thujaplicin_ht15_bivar_model1)
asreml::vpredict(gamma_thujaplicin_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(gamma_thujaplicin_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(gamma_thujaplicin_ht15_bivar_model3)
asreml::vpredict(gamma_thujaplicin_ht15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(gamma_thujaplicin_ht15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# gthujaplicin - dbh15
gamma_thujaplicin_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_dbh15_bivar_model1$G.param, 
                                             R.param = gamma_thujaplicin_dbh15_bivar_model1$R.param)

gamma_thujaplicin_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ id(units):us(trait),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = gamma_thujaplicin_dbh15_bivar_model3$G.param, 
                                             R.param = gamma_thujaplicin_dbh15_bivar_model3$R.param)

gamma_thujaplicin_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

gamma_thujaplicin_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(gamma_thujaplicin, dbh15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model2, gamma_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model3, gamma_thujaplicin_dbh15_bivar_model1)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model2, gamma_thujaplicin_dbh15_bivar_model3)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model3, gamma_thujaplicin_dbh15_bivar_model4)
REMLRT.asreml(gamma_thujaplicin_dbh15_bivar_model3, gamma_thujaplicin_dbh15_bivar_model5)

varcomps(gamma_thujaplicin_dbh15_bivar_model1)
asreml::vpredict(gamma_thujaplicin_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(gamma_thujaplicin_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(gamma_thujaplicin_dbh15_bivar_model3)
asreml::vpredict(gamma_thujaplicin_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(gamma_thujaplicin_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# bthujaplicinol - ht15
beta_thujaplicinol_ht15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06, G.param = beta_thujaplicinol_ht15_bivar_model1$G.param, 
                                             R.param = beta_thujaplicinol_ht15_bivar_model1$R.param)

beta_thujaplicinol_ht15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ id(units):us(trait),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

beta_thujaplicinol_ht15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicinol_ht15_bivar_model3$G.param, 
                                             R.param = beta_thujaplicinol_ht15_bivar_model3$R.param)

beta_thujaplicinol_ht15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

beta_thujaplicinol_ht15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                             fixed = cbind(beta_thujaplicinol, ht15) ~ trait + trait:site,
                                             random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                             residual = ~ dsum(~ id(units):us(trait) | site),
                                             na.action = na.method(x = "include", y = "include"),
                                             maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model2, beta_thujaplicinol_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model3, beta_thujaplicinol_ht15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model2, beta_thujaplicinol_ht15_bivar_model3)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model3, beta_thujaplicinol_ht15_bivar_model4)
REMLRT.asreml(beta_thujaplicinol_ht15_bivar_model3, beta_thujaplicinol_ht15_bivar_model5)

varcomps(beta_thujaplicinol_ht15_bivar_model1)
asreml::vpredict(beta_thujaplicinol_ht15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicinol_ht15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(beta_thujaplicinol_ht15_bivar_model3)
asreml::vpredict(beta_thujaplicinol_ht15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicinol_ht15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# bthujaplicinol - dbh15
beta_thujaplicinol_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicinol_dbh15_bivar_model1$G.param, 
R.param = beta_thujaplicinol_dbh15_bivar_model1$R.param)

beta_thujaplicinol_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                               residual = ~ id(units):us(trait),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)

beta_thujaplicinol_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)#, G.param = beta_thujaplicinol_dbh15_bivar_model3$G.param, 
R.param = beta_thujaplicinol_dbh15_bivar_model3$R.param)

beta_thujaplicinol_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)

beta_thujaplicinol_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                               fixed = cbind(beta_thujaplicinol, dbh15) ~ trait + trait:site,
                                               random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                               residual = ~ dsum(~ id(units):us(trait) | site),
                                               na.action = na.method(x = "include", y = "include"),
                                               maxit = 100, workspace = 80e+06)

REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model2, beta_thujaplicinol_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model3, beta_thujaplicinol_dbh15_bivar_model1)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model2, beta_thujaplicinol_dbh15_bivar_model3)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model3, beta_thujaplicinol_dbh15_bivar_model4)
REMLRT.asreml(beta_thujaplicinol_dbh15_bivar_model3, beta_thujaplicinol_dbh15_bivar_model5)

varcomps(beta_thujaplicinol_dbh15_bivar_model1)
asreml::vpredict(beta_thujaplicinol_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicinol_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(beta_thujaplicinol_dbh15_bivar_model3)
asreml::vpredict(beta_thujaplicinol_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(beta_thujaplicinol_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595

# ht15 - dbh15
ht15_dbh15_bivar_model1 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)#, G.param = ht15_dbh15_bivar_model1$G.param, 
R.param = ht15_dbh15_bivar_model1$R.param)

ht15_dbh15_bivar_model2 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:idh(site):REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                residual = ~ id(units):us(trait),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)

ht15_dbh15_bivar_model3 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:site:REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)#, G.param = ht15_dbh15_bivar_model3$G.param, 
R.param = ht15_dbh15_bivar_model3$R.param)

ht15_dbh15_bivar_model4 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:site:REP + trait:idh(site):REP:SET + us(trait):vm(tree, wrc.ainv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)

ht15_dbh15_bivar_model5 <- asreml(data = wrc.pheno.dat,
                                                fixed = cbind(ht15, dbh15) ~ trait + trait:site,
                                                random = ~ trait:idh(site):REP + trait:site:REP:SET + us(trait):vm(tree, wrc.ainv),
                                                residual = ~ dsum(~ id(units):us(trait) | site),
                                                na.action = na.method(x = "include", y = "include"),
                                                maxit = 100, workspace = 80e+06)

REMLRT.asreml(ht15_dbh15_bivar_model2, ht15_dbh15_bivar_model1)
REMLRT.asreml(ht15_dbh15_bivar_model3, ht15_dbh15_bivar_model1)
REMLRT.asreml(ht15_dbh15_bivar_model2, ht15_dbh15_bivar_model3)
REMLRT.asreml(ht15_dbh15_bivar_model3, ht15_dbh15_bivar_model4)
REMLRT.asreml(ht15_dbh15_bivar_model3, ht15_dbh15_bivar_model5)

varcomps(ht15_dbh15_bivar_model1)
asreml::vpredict(ht15_dbh15_bivar_model1, gc ~ V8 / sqrt(V7 * V9)) # 0.2187494 +/-0.1793586
asreml::vpredict(ht15_dbh15_bivar_model1, pc ~ (V8 + (V12 + V16 + V20)/3)/ sqrt((V7 +(V11 + V15 + V19)/3) * (V9 + (V13 + V17 + V21)/3))) # 0.08814571 +/-0.02860595

varcomps(ht15_dbh15_bivar_model3)
asreml::vpredict(ht15_dbh15_bivar_model3, gc ~ V4 / sqrt(V3 * V5)) # 0.2187494 +/-0.1793586
asreml::vpredict(ht15_dbh15_bivar_model3, pc ~ (V4 + (V8 + V12 + V16)/3)/ sqrt((V3 +(V7 + V11 + V15)/3) * (V5 + (V9 + V13 + V17)/3))) # 0.08814571 +/-0.02860595
