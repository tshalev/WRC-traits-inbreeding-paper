rm(list = ls())

library(asreml)
library(asremlPlus)
library(nadiv)
library(synbreed)
library(sommer)
library(genio)
library(ggsci)
library(tidyverse)

source("publication_theme.r")

library(extrafont)
library(extrafontdb)

wrc.pheno.dat <- read.csv("master_phenotype_dataset.csv", header = T)
S_lines_ibc <- read.table("corrected_snps_S_lines_1654_r201_mac3.ibc", header = T)

wrc.pheno.dat[1:14] <- lapply(wrc.pheno.dat[1:14], as.factor)
wrc.pheno.dat[15:74] <- lapply(wrc.pheno.dat[15:74], as.numeric)

S_lines_ibc[c(1,2,7,9,10,12,13)] <- lapply(S_lines_ibc[c(1,2,7,9,10,12,13)], as.factor)

# Pedigree
wrc.pheno.ped <- read.csv("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/cedar/FINAL_DATA_SET/filtered_normalised/PED_FINAL_training_only_3_cols.ped.csv", header = F, na.strings = "NA") 
head(wrc.pheno.ped)

wrc.ainv <- ainverse(wrc.pheno.ped)

# S line relationship matrix
S_lines_ind_names <- read.table("S_line_ind_names.txt")
S_lines_loci_names <- read.table("S_line_loci_names.txt")
S_lines_012 <- read_bed("corrected_snps_S_lines_binary.bed", names_loci = S_lines_loci_names$V1, names_ind = S_lines_ind_names$V1)

# Subtract to make -101
S_lines_012 <- S_lines_012 - 1
S_lines_012 <- t(S_lines_012)
S_lines_012[1:10,1:2]

# Make G matrix
G <- A.mat(S_lines_012)

# Make G near positive definite
require(Matrix)
RealizedPD <- nearPD(G, keepDiag = T)
G <- matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
G <- G + diag(0.01, nrow(G))
attr(G, "dimnames") = RealizedPD[[1]]@Dimnames
class(G) = "relationshipMatrix"

# Make G a relationship matrix
S_lines.giv <- write.relationshipMatrix(G, sorting = "ASReml", type = "ginv", file = NULL)
head(attr(S_lines.giv, "rowNames"))
names(S_lines.giv) <- c("row", "column", "coefficient")
head(S_lines.giv)
attr(S_lines.giv, "INVERSE") <- TRUE

### Plots (including Figure 3) ###
fhat_gen_lm <- lm(Fhat3 ~ generation, data = S_lines_ibc)
summary(fhat_gen_lm)

id_gen_plot <- ggplot(S_lines_ibc, aes(x = generation, y = Fhat3, colour = type, group = type)) +
  geom_point() +
  geom_smooth(se = F) +
  theme_Publication() +
  scale_colour_nejm() +
  ylab(expression(bolditalic("F"))) +
  xlab("Generation") +
  theme(legend.position = "none")
  

id_gen_plot

# ggsave("F_gen_plot.tiff", height = 8, width = 9, dpi = 300)

S_lines_ibc_no_S5 <- S_lines_ibc %>% 
  filter(!(generation == "S5"))

ht_gen_lm <- lm(HT_BLUP ~ generation, data = S_lines_ibc_no_S5)
summary(ht_gen_lm)

ht_gen_plot <- ggplot(S_lines_ibc_no_S5, aes(x = generation, y = HT_BLUP, colour = type, group = type)) +
  geom_point() +
  geom_abline(slope = 0.004406, intercept = -0.0003571, linetype = "solid", colour= "black", size = 1) +
  stat_smooth(se = F) +
  theme_Publication() +
  # scale_colour_nejm() +
  ylab("Height BV") +
  xlab("Generation") +
  scale_colour_manual(values = col, name = "Line Type", labels = c("Random", "Select")) +
  theme(legend.position = "none") 

ht_gen_plot

ggsave("ht_gen_plot_no_legend_reg_line.tiff", height = 8, width = 9, dpi = 300)

ht_bv_gen_boxplot <- ggplot(S_lines_ibc, aes(x = generation, y = HT_BLUP, fill = type)) +
  geom_boxplot() +
  theme_Publication() +
  scale_fill_nejm()

ht_bv_gen_boxplot

ht_bv_line_gen_boxplot <- ggplot(S_lines_ibc, aes(x = line, y = Fhat3, fill = type)) +
  geom_boxplot() +
  theme_Publication() +
  scale_fill_nejm()

ht_bv_line_gen_boxplot


### Models for training population - inbreeding depression ########

## height

ht15_mlm <- asreml(data = wrc.pheno.dat,
                   fixed = ht15 ~ site + Fhat3,
                   random = ~ idh(site):REP + idh(site):REP:SET,
                   residual = ~ dsum(~ id(units) | site),
                   na.action = na.method(x = "include", y = "include"),
                   maxit = 100, workspace = 80e+06)
summary(ht15_mlm)

wald.asreml(ht15_mlm)

summary(ht15_mlm, coef = T)$coef.fixed

## DBH

dbh15_mlm <- asreml(data = wrc.pheno.dat,
                   fixed = dbh15 ~ site + Fhat3,
                   random = ~ site:REP + site:REP:SET,
                   residual = ~ dsum(~ id(units) | site),
                   na.action = na.method(x = "include", y = "include"),
                   maxit = 100, workspace = 80e+06)
summary(dbh15_mlm)

wald.asreml(dbh15_mlm)

summary(dbh15_mlm, coef = T)$coef.fixed

## heartwood radius

heartwood_mlm <- asreml(data = wrc.pheno.dat,
                   fixed = heartwood ~ site + dbh15 + pith_age + Fhat3,
                   random = ~ idh(site):REP + idh(site):REP:SET,
                   residual = ~ dsum(~ id(units) | site),
                   na.action = na.method(x = "include", y = "include"),
                   maxit = 100, workspace = 80e+06)

summary(heartwood_mlm)

wald.asreml(heartwood_mlm)

summary(heartwood_mlm, coef = T)$coef.fixed

## sapwood radius
sapwood_mlm <- asreml(data = wrc.pheno.dat,
                        fixed = sapwood ~ site + dbh15 + pith_age + Fhat3,
                        random = ~ idh(site):REP + idh(site):REP:SET,
                        residual = ~ dsum(~ id(units) | site),
                        na.action = na.method(x = "include", y = "include"),
                        maxit = 100, workspace = 80e+06)

summary(sapwood_mlm)

wald.asreml(sapwood_mlm)

summary(sapwood_mlm, coef = T)$coef.fixed

## sapwood-heartwood ratio radius
sapwood_heartwood_mlm <- asreml(data = wrc.pheno.dat,
                      fixed = sapwood_heartwood_ratio ~ site + dbh15 + pith_age + Fhat3,
                      random = ~ idh(site):REP + idh(site):REP:SET,
                      residual = ~ dsum(~ id(units) | site),
                      na.action = na.method(x = "include", y = "include"),
                      maxit = 100, workspace = 80e+06)

summary(sapwood_heartwood_mlm)

wald.asreml(sapwood_heartwood_mlm)

summary(sapwood_heartwood_mlm, coef = T)$coef.fixed

## total foliar monoterpenes
fmono_mlm <- asreml(data = wrc.pheno.dat,
                      fixed = fmono ~ site + Fhat3,
                      random = ~ idh(site):REP + site:REP:SET,
                      residual = ~ dsum(~ id(units) | site),
                      na.action = na.method(x = "include", y = "include"),
                      maxit = 100, workspace = 80e+06)

summary(fmono_mlm)

wald.asreml(fmono_mlm)

summary(fmono_mlm, coef = T)$coef.fixed

## total thujaplicins
thujaplicins_mlm <- asreml(data = wrc.pheno.dat,
                      fixed = thujaplicins ~ site + dbh15 + pith_age + Fhat3,
                      random = ~ idh(site):REP + site:REP:SET,
                      residual = ~ dsum(~ id(units) | site),
                      na.action = na.method(x = "include", y = "include"),
                      maxit = 100, workspace = 80e+06)

summary(thujaplicins_mlm)

wald.asreml(thujaplicins_mlm)

summary(thujaplicins_mlm, coef = T)$coef.fixed


### HT BLUP models for S lines #################
S_lines_ibc_naomit <- S_lines_ibc %>% 
  filter(!(is.na(HT_BLUP)))

S_lines_select <- S_lines_ibc_naomit %>% 
  filter(type =="s")

S_lines_random <- S_lines_ibc_naomit %>% 
  filter(type =="r")

wrc_naomit <- wrc.pheno.dat %>% 
  filter(!(is.na(HT_BLUP)))

S_lines_HT_BLUP_F <- lm(HT_BLUP ~ Fhat3, data = S_lines_ibc_naomit)

mlm_training_HT_BLUP_F <- asreml(data = wrc.pheno.dat,
                     fixed = HT_BLUP ~ Fhat3,
                     residual = ~ dsum(~ id(units) | site),
                     na.action = na.method(x = "include", y = "include"),
                     maxit = 100, workspace = 80e+06)


anova(S_lines_HT_BLUP_F)

wald.asreml(mlm_training_HT_BLUP_F)
summary(mlm_training_HT_BLUP_F, coef = T)$coef.fixed

m <- lm(HT_BLUP ~ Fhat3, data = S_lines_ibc_naomit)

m2 <- lm(HT_BLUP ~ Fhat3 + type, data = S_lines_ibc_naomit)

select_model <- lm(HT_BLUP ~ , data = S_lines_select)
random_model <- lm(HT_BLUP ~ Fhat3 + generation, data = S_lines_random)

anova(m, m2)

m2  # get the estimated coefficients for the model
anova(m2)   # get the Anova table for the model, including the F test
anova(select_model)   # get the Anova table for the model, including the F test
car::Anova(select_model, type = "III")
summary(select_model)

anova(random_model)   # get the Anova table for the model, including the F test
car::Anova(random_model, type = "III")
summary(random_model)


summary(m) # get t-tests and some fit statistics for the model
ibc_naomit$yhat.m<-fitted(m)  # the estimated y values
ibc_naomit$resid.m<-resid(m)

# get diagnostic plots
par(mfrow=c(2,2),mai=c(0.6,0.6,0.6,0.6),cex=0.7)
plot(ibc_naomit$yhat.m,ibc_naomit$resid.m, main="Model 3, Residual Plot",
     xlab="yhat", ylab="residual")
plot(ibc_naomit$HT_BLUP,ibc_naomit$yhat.m, main="Model 3, Fitted line plot",
     ylab="yhat", xlab="volume")
qqnorm(ibc_naomit$resid.m, main="Model 3, Normality plot")
hist(ibc_naomit$resid.m, breaks =8 , density=10,col="green", border="black",
     main="Model 3, Error Distribution") 
par(mfrow=c(1,1),mai=c(1.0,1.0,1.0,1.0),cex=1.0)

# get the normality test results
shapiro.test(ibc_naomit$resid.m) # Shapiro-Wilk normality test

# Use Breusch-Pagan test for unequal variances
require(car)
ncvTest(m)

interaction.plot(ibc_naomit$Fhat3, ibc_naomit$type, ibc_naomit$yhat.m)


### Get HT BLUPs from pedigree ######
ht15_ABLUP_model <- asreml(data = wrc.pheno.dat,
                    fixed = ht15 ~ site,
                    random = ~ idh(site):REP + idh(site):REP:SET + vm(tree, wrc.ainv),
                    residual = ~ dsum(~ id(units) | site),
                    na.action = na.method(x = "include", y = "include"),
                    maxit = 100, workspace = 80e+06)



# BVs
ht15_solutions <- rbind(data.frame(summary(ht15_ABLUP_model, coef = T)$coef.fixed),
                             data.frame(summary(ht15_ABLUP_model, coef = T)$coef.random))

str(ht15_solutions)
head(ht15_solutions, 30)

ht15_solutions <- tibble::rownames_to_column(ht15_solutions, "model_term")

ht15_solutions_ped_tree <- ht15_solutions %>% 
  filter(str_detect(model_term, 'vm\\(tree, wrc.ainv\\)\\_[0-9]') == TRUE)

(mean_accuracy_ht15_univariate <- mean(sqrt(1-(ht15_solutions_ped_tree$std.error^2/varcomps(ht15_ABLUP_model)$component[5])))) # 0.6136018

ht15_BVs_df <- data.frame(ID = sub("vm(tree, wrc.ainv)_", "", ht15_solutions_ped_tree$model_term, fixed = T),
                               BV = (ht15_solutions_ped_tree$solution))

write.table(ht15_BVs_df, "ht15_ABLUP_BVs_ASREML4.txt", quote = F, row.names = F)


#### Models to get BVs for S lines ######################################

## total monoterpenes GBLUP model ####
total_mono_model <- asreml(data = S_lines_ibc,
                          fixed = total_mono ~ location + generation,
                          random = ~ vm(IID, S_lines.giv),
                          residual = ~ id(units),
                          na.action = na.method(x = "include", y = "include"),
                          maxit = 100, workspace = 80e+06)

summary(total_mono_model)
wald.asreml(total_mono_model)

# BVs
total_mono_solutions <- rbind(data.frame(summary.asreml(total_mono_model, coef = T)$coef.fixed),
                             data.frame(summary.asreml(total_mono_model, coef = T)$coef.random))

str(total_mono_solutions)
head(total_mono_solutions, 30)

total_mono_solutions <- tibble::rownames_to_column(total_mono_solutions, "model_term")

total_mono_solutions_ped_tree <- total_mono_solutions %>% 
  filter(str_detect(model_term, 'vm\\(IID, S_lines.giv\\)\\_[0-9]') == TRUE)

(mean_accuracy_total_mono_univariate <- mean(sqrt(1-(total_mono_solutions_ped_tree$std.error^2/varcomps(total_mono_model)$component[2])))) # 0.7125621

total_mono_BVs_df <- data.frame(ID = sub("vm(IID, S_lines.giv)_", "", total_mono_solutions_ped_tree$model_term, fixed = T),
                               BV = (total_mono_solutions_ped_tree$solution))

# write.table(total_mono_BVs_df, "total_mono_GBLUP_BVs_ASREML4.txt", quote = F, row.names = F)

total_mono_BVs_df <- total_mono_BVs_df %>% 
  arrange(ID)

S_lines_ibc <- S_lines_ibc %>% 
  arrange(IID) %>% 
  mutate(total_mono_BV = total_mono_BVs_df$BV)


mono_BV_lm <- lm(total_mono_BV ~ Fhat3, S_lines_ibc)
mono_BV_lm2 <- lm(total_mono_BV ~ location + line + Fhat3, S_lines_ibc)

mono_BV_lm
anova(mono_BV_lm, mono_BV_lm2)

anova(mono_BV_lm)
S_lines_ibc$yhat <- fitted(mono_BV_lm)
interaction.plot()


## Distribution of Fhat3
ggplot(wrc.pheno.dat, aes(x = Fhat3)) +
  geom_histogram(bins = 100) +
  xlab(expression(bolditalic("F"))) +
  ylab("Count") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

ggsave("Fhat3_dist.tiff", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = Fhat3, y = ht15)) +
  geom_point(alpha = 0.5) +
  xlab(expression(bolditalic("F"))) +
  ylab("Height (cm)") +
  theme_Publication()

ggsave("height_fhat_scatter.tiff", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = Fhat3, y = fmono)) +
  geom_point(alpha = 0.5) +
  xlab(expression(bolditalic("F"))) +
  ylab(expression(bold(paste("Total foliar monoterpenes (", mu, "g/g DW)")))) + 
  theme_Publication()

ggsave("fmono_fhat_scatter.tiff", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = Fhat3, y = thujaplicins)) +
  geom_point(alpha = 0.5) +
  xlab(expression(bolditalic("F"))) +
  ylab(expression(bold(paste("Total wood thujaplicins (", mu, "g/g CW)")))) + 
  theme_Publication()
  
ggsave("total_thuj_fhat_scatter.tiff", dpi = 300, width = 9, height = 8)
