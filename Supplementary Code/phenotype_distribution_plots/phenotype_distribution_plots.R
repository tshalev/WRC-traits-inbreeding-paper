rm(list=ls())

library(tidyverse)

source("publication_theme.r")

# Read phenotype data
wrc.pheno.dat <- read.table("master_phenotype_dataset.csv", na.strings = "NA", header = T, sep = ",")

col <- c("#BC3C29FF", "#0072B5FF", "#E18727FF",  "#7876B1FF", "#6F99ADFF", "#FFDC91FF", "#EE4C97FF", "#20854EFF")

wrc.pheno.dat[1:14] <- lapply(wrc.pheno.dat[1:14], as.factor)
wrc.pheno.dat[15:74] <- lapply(wrc.pheno.dat[15:74], as.numeric)


### Plot traits by sites ###############


# boxplot_df <- phen %>% 
#   select(site, ht15, dbh15, fath, fmono, Alpha.Thujaplicin, Sum.Thuj, Sum.PA.B, Sum.import.6)
# 
# melted_boxplot_df <- reshape2::melt(boxplot_df)

# ggplot(melted_boxplot_df, aes(x = site, y = value)) +
#   geom_violin() +
#   geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
#   facet_wrap(~ variable) +
#   theme_Publication()


ggplot(wrc.pheno.dat, aes(x = site, y = ht15)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab("Height (cm)") +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("ht15_site_distribution.tiff", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = dbh15)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab("DBH (mm)") +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("DBH15_site_distribution.tiff", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = fath)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab(expression(bold(paste(alpha, "-thujone (", mu, "g/g DW)")))) +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("a-thujone_site_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = fmono)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab(expression(bold(paste("Total foliar monoterpenes (", mu, "g/g DW)")))) +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("fmono_site_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = alpha_thujaplicin)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab(expression(bold(paste(alpha, "-thujaplicin (", mu, "g/g CW)")))) +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("a-thujaplicin_site_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = Sum.Thuj)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab(expression(bold(paste("Total thujaplicins (", mu, "g/g CW)")))) + 
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("total_thujaplicins_site_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = Sum.PA.B)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab("Total lignans (PAR/g CW)") +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("total_lignans_site_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = site, y = Sum.import.6)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.alpha = 0.8) +
  ylab("Total wood (PAR/g CW)") +
  scale_x_discrete(name = "Site", labels = c("Jordan River", "Powell River", "Port McNeill")) +
  theme_Publication()

# ggsave("total_wood_site_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = ht15, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab("Height (cm)") +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("ht15_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = dbh15, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab("DBH (mm)") +
  ylab("Density") +
  theme_Publication() +
  scale_fill_manual(name="Site", labels=c("Jordan River", "Powell River", "Port McNeill"), values = col) +
  theme(legend.position = c(0.9, 0.92), legend.background = element_rect(colour = "black"))

# ggsave("dbh15_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = fath, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste(alpha, "-thujone (", mu, "g/g DW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("a-thujone_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = fbth, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste(beta, "-thujone (", mu, "g/g DW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("fbth_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = fsab, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste("Sabinene (", mu, "g/g DW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("fsab_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = alpha_thujaplicin, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste(alpha, "-thujaplicin (", mu, "g/g CW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("a-thujaplicin_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = beta_thujaplicin, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste(beta, "-thujaplicin (", mu, "g/g CW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("b-thujaplicin_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = gamma_thujaplicin, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste(gamma, "-thujaplicin (", mu, "g/g CW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("g-thujaplicin_phen_distribution.svg", dpi = 300, width = 9, height = 8)

ggplot(wrc.pheno.dat, aes(x = beta_thujaplicinol, fill = site)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste(beta, "-thujaplicinol (", mu, "g/g CW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_nejm() +
  theme(legend.position = "none")

# ggsave("b-thujaplicinol_phen_distribution.svg", dpi = 300, width = 9, height = 8)


### Height/DBH summaries / significance tests #######
as.data.frame(wrc.pheno.dat %>% 
  group_by(site) %>% 
  summarise(mean = mean(ht15), median = median(ht15), sd = sd(ht15)))

as.data.frame(wrc.pheno.dat %>% 
  group_by(site) %>% 
  summarise(mean = mean(dbh15), median = median(dbh15), sd = sd(dbh15)))

site_1_ht_dbh <- wrc.pheno.dat %>% 
  select(tree, site, ht15, dbh15) %>% 
  filter(site == 1)

site_2_ht_dbh <- wrc.pheno.dat %>% 
  select(tree, site, ht15, dbh15) %>% 
  filter(site == 2)

site_3_ht_dbh <- wrc.pheno.dat %>% 
  select(tree, site, ht15, dbh15) %>% 
  filter(site == 3)

t.test(site_2_ht_dbh$ht15, site_1_ht_dbh$ht15, alternative = "greater")
t.test(site_2_ht_dbh$ht15, site_3_ht_dbh$ht15, alternative = "greater")

t.test(site_2_ht_dbh$dbh15, site_1_ht_dbh$dbh15, alternative = "greater")
t.test(site_2_ht_dbh$dbh15, site_3_ht_dbh$dbh15, alternative = "greater")

## heritability plot ####
heritabilities <- read.csv("heritabilities.csv")

heritabilities$Trait <- factor(heritabilities$Trait, levels = c("a-thujone", "b-thujone", "sabinene", "a-thujaplicin", 
                                                                "b-thujaplicin", "g-thujaplicin", "b-thujaplicinol", "Height", "DBH"))

ggplot(heritabilities, aes(x = Trait, y = G_h2)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = G_h2 - G_h2_SE, ymax = G_h2 + G_h2_SE), width = 0.1) +
  ylab(expression(italic(bolditalic(h)^2))) +
  scale_x_discrete("Trait",
                   labels = c(expression(paste(alpha, "-thujone")), expression(paste(beta, "-thujone")), "Sabinene", expression(paste(alpha, "-thujaplicin")),
                              expression(paste(beta, "-thujaplicin")), expression(paste(gamma, "-thujaplicin")), expression(paste(beta, "-thujaplicinol")), "Height", "DBH")) +
  theme_Publication()

# ggsave("h2_plot.tiff", dpi = 300, width = 14, height = 8)

## density plots ###
melted_foliar <- melt(wrc.pheno.dat, id.vars = "tree", measure.vars = c("fath", "fbth", "fsab"))

ggplot(melted_foliar, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste("Monoterpenes (", "\u03bc", "g/g DW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_manual(values = col, name = NULL, labels = c(expression(paste(alpha, "-thujone")), expression(paste(beta, "-thujone")), "Sabinene"))
    # theme(legend.position = "none")

# ggsave("foliar_terpene_density_plot.tiff", dpi = 300, width = 9, height = 8)

melted_wood <- melt(wrc.pheno.dat, id.vars = "tree", measure.vars = c("alpha_thujaplicin", "beta_thujaplicin", "gamma_thujaplicin", "beta_thujaplicinol"))

ggplot(melted_wood, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  xlab(expression(bold(paste("Monoterpenes (", "\u03bc", "g/g CW)")))) +
  ylab("Density") +
  theme_Publication() +
  scale_fill_manual(values = col, name = NULL, labels = c(expression(paste(alpha, "-thujaplicin")), expression(paste(beta, "-thujaplicin")), expression(paste(gamma, "-thujaplicin")), expression(paste(beta, "-thujaplicinol"))))
# theme(legend.position = "none")

# ggsave("wood_terpene_density_plot.tiff", dpi = 300, width = 9, height = 8)

