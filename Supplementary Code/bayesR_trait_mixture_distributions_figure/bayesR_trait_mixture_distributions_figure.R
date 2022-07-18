rm(list=ls())

library(extrafont)
library(wesanderson)
library(ghibli)
library(ggsci)
library(reshape2)
library(tidyverse)
source("../publication_theme.r")

font_import()
loadfonts(device = "win")
windowsFonts()

# Read in model summaries
athuj_adj_phen_model_output_chain_1 <- read.table("bayesR_athujone_adj_phenotype/bayesR_athujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
athuj_adj_phen_model_output_chain_2 <- read.table("bayesR_athujone_adj_phenotype/bayesR_athujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
athuj_adj_phen_model_output_chain_3 <- read.table("bayesR_athujone_adj_phenotype/bayesR_athujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
athuj_adj_phen_model_output_chain_4 <- read.table("bayesR_athujone_adj_phenotype/bayesR_athujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
athuj_adj_phen_model_output_chain_5 <- read.table("bayesR_athujone_adj_phenotype/bayesR_athujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

bthuj_adj_phen_model_output_chain_1 <- read.table("bayesR_bthujone_adj_phenotype/bayesR_bthujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
bthuj_adj_phen_model_output_chain_2 <- read.table("bayesR_bthujone_adj_phenotype/bayesR_bthujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
bthuj_adj_phen_model_output_chain_3 <- read.table("bayesR_bthujone_adj_phenotype/bayesR_bthujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
bthuj_adj_phen_model_output_chain_4 <- read.table("bayesR_bthujone_adj_phenotype/bayesR_bthujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
bthuj_adj_phen_model_output_chain_5 <- read.table("bayesR_bthujone_adj_phenotype/bayesR_bthujone_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

sabinene_adj_phen_model_output_chain_1 <- read.table("bayesR_sabinene_adj_phenotype/bayesR_sabinene_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
sabinene_adj_phen_model_output_chain_2 <- read.table("bayesR_sabinene_adj_phenotype/bayesR_sabinene_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
sabinene_adj_phen_model_output_chain_3 <- read.table("bayesR_sabinene_adj_phenotype/bayesR_sabinene_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
sabinene_adj_phen_model_output_chain_4 <- read.table("bayesR_sabinene_adj_phenotype/bayesR_sabinene_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
sabinene_adj_phen_model_output_chain_5 <- read.table("bayesR_sabinene_adj_phenotype/bayesR_sabinene_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

athujaplicin_adj_phen_model_output_chain_1 <- read.table("bayesR_athujaplicin_adj_phenotype/bayesR_athujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
athujaplicin_adj_phen_model_output_chain_2 <- read.table("bayesR_athujaplicin_adj_phenotype/bayesR_athujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
athujaplicin_adj_phen_model_output_chain_3 <- read.table("bayesR_athujaplicin_adj_phenotype/bayesR_athujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
athujaplicin_adj_phen_model_output_chain_4 <- read.table("bayesR_athujaplicin_adj_phenotype/bayesR_athujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
athujaplicin_adj_phen_model_output_chain_5 <- read.table("bayesR_athujaplicin_adj_phenotype/bayesR_athujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

bthujaplicin_adj_phen_model_output_chain_1 <- read.table("bayesR_bthujaplicin_adj_phenotype/bayesR_bthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
bthujaplicin_adj_phen_model_output_chain_2 <- read.table("bayesR_bthujaplicin_adj_phenotype/bayesR_bthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
bthujaplicin_adj_phen_model_output_chain_3 <- read.table("bayesR_bthujaplicin_adj_phenotype/bayesR_bthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
bthujaplicin_adj_phen_model_output_chain_4 <- read.table("bayesR_bthujaplicin_adj_phenotype/bayesR_bthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
bthujaplicin_adj_phen_model_output_chain_5 <- read.table("bayesR_bthujaplicin_adj_phenotype/bayesR_bthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

gthujaplicin_adj_phen_model_output_chain_1 <- read.table("bayesR_gthujaplicin_adj_phenotype/bayesR_gthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
gthujaplicin_adj_phen_model_output_chain_2 <- read.table("bayesR_gthujaplicin_adj_phenotype/bayesR_gthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
gthujaplicin_adj_phen_model_output_chain_3 <- read.table("bayesR_gthujaplicin_adj_phenotype/bayesR_gthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
gthujaplicin_adj_phen_model_output_chain_4 <- read.table("bayesR_gthujaplicin_adj_phenotype/bayesR_gthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
gthujaplicin_adj_phen_model_output_chain_5 <- read.table("bayesR_gthujaplicin_adj_phenotype/bayesR_gthujaplicin_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

bthujaplicinol_adj_phen_model_output_chain_1 <- read.table("bayesR_bthujaplicinol_adj_phenotype/bayesR_bthujaplicinol_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
bthujaplicinol_adj_phen_model_output_chain_2 <- read.table("bayesR_bthujaplicinol_adj_phenotype/bayesR_bthujaplicinol_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
bthujaplicinol_adj_phen_model_output_chain_3 <- read.table("bayesR_bthujaplicinol_adj_phenotype/bayesR_bthujaplicinol_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
bthujaplicinol_adj_phen_model_output_chain_4 <- read.table("bayesR_bthujaplicinol_adj_phenotype/bayesR_bthujaplicinol_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
bthujaplicinol_adj_phen_model_output_chain_5 <- read.table("bayesR_bthujaplicinol_adj_phenotype/bayesR_bthujaplicinol_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

ht15_adj_phen_model_output_chain_1 <- read.table("bayesR_ht15_noclb_adj_phenotype/bayesR_ht15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
ht15_adj_phen_model_output_chain_2 <- read.table("bayesR_ht15_noclb_adj_phenotype/bayesR_ht15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
ht15_adj_phen_model_output_chain_3 <- read.table("bayesR_ht15_noclb_adj_phenotype/bayesR_ht15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
ht15_adj_phen_model_output_chain_4 <- read.table("bayesR_ht15_noclb_adj_phenotype/bayesR_ht15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
ht15_adj_phen_model_output_chain_5 <- read.table("bayesR_ht15_noclb_adj_phenotype/bayesR_ht15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

dbh15_adj_phen_model_output_chain_1 <- read.table("bayesR_dbh15_noclb_adj_phenotype/bayesR_dbh15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_1.model")
dbh15_adj_phen_model_output_chain_2 <- read.table("bayesR_dbh15_noclb_adj_phenotype/bayesR_dbh15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_2.model")
dbh15_adj_phen_model_output_chain_3 <- read.table("bayesR_dbh15_noclb_adj_phenotype/bayesR_dbh15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_3.model")
dbh15_adj_phen_model_output_chain_4 <- read.table("bayesR_dbh15_noclb_adj_phenotype/bayesR_dbh15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_4.model")
dbh15_adj_phen_model_output_chain_5 <- read.table("bayesR_dbh15_noclb_adj_phenotype/bayesR_dbh15_adj_phen_300000it_100000burnin_10thin_4dist_permute_chain_5.model")

# Heritabilities
heritabilities <- read.csv("heritabilities.csv")

# Add mean column
athuj_all_chains <- data.frame(cbind(chain_1 = athuj_adj_phen_model_output_chain_1, chain_2 = athuj_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = athuj_adj_phen_model_output_chain_3$V2, chain_4 = athuj_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = athuj_adj_phen_model_output_chain_5$V2), row.names = 1)

athuj_all_chains <- athuj_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(athuj_all_chains)) %>% 
  column_to_rownames()

bthuj_all_chains <- data.frame(cbind(chain_1 = bthuj_adj_phen_model_output_chain_1, chain_2 = bthuj_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = bthuj_adj_phen_model_output_chain_3$V2, chain_4 = bthuj_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = bthuj_adj_phen_model_output_chain_5$V2), row.names = 1)

bthuj_all_chains <- bthuj_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(bthuj_all_chains)) %>% 
  column_to_rownames()

sabinene_all_chains <- data.frame(cbind(chain_1 = sabinene_adj_phen_model_output_chain_1, chain_2 = sabinene_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = sabinene_adj_phen_model_output_chain_3$V2, chain_4 = sabinene_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = sabinene_adj_phen_model_output_chain_5$V2), row.names = 1)

sabinene_all_chains <- sabinene_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(sabinene_all_chains)) %>% 
  column_to_rownames()

athujaplicin_all_chains <- data.frame(cbind(chain_1 = athujaplicin_adj_phen_model_output_chain_1, chain_2 = athujaplicin_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = athujaplicin_adj_phen_model_output_chain_3$V2, chain_4 = athujaplicin_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = athujaplicin_adj_phen_model_output_chain_5$V2), row.names = 1)

athujaplicin_all_chains <- athujaplicin_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(athujaplicin_all_chains)) %>% 
  column_to_rownames()

bthujaplicin_all_chains <- data.frame(cbind(chain_1 = bthujaplicin_adj_phen_model_output_chain_1, chain_2 = bthujaplicin_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = bthujaplicin_adj_phen_model_output_chain_3$V2, chain_4 = bthujaplicin_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = bthujaplicin_adj_phen_model_output_chain_5$V2), row.names = 1)

bthujaplicin_all_chains <- bthujaplicin_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(bthujaplicin_all_chains)) %>% 
  column_to_rownames()

gthujaplicin_all_chains <- data.frame(cbind(chain_1 = gthujaplicin_adj_phen_model_output_chain_1, chain_2 = gthujaplicin_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = gthujaplicin_adj_phen_model_output_chain_3$V2, chain_4 = gthujaplicin_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = gthujaplicin_adj_phen_model_output_chain_5$V2), row.names = 1)

gthujaplicin_all_chains <- gthujaplicin_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(gthujaplicin_all_chains)) %>% 
  column_to_rownames()

bthujaplicinol_all_chains <- data.frame(cbind(chain_1 = bthujaplicinol_adj_phen_model_output_chain_1, chain_2 = bthujaplicinol_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = bthujaplicinol_adj_phen_model_output_chain_3$V2, chain_4 = bthujaplicinol_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = bthujaplicinol_adj_phen_model_output_chain_5$V2), row.names = 1)

bthujaplicinol_all_chains <- bthujaplicinol_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(bthujaplicinol_all_chains)) %>% 
  column_to_rownames()

ht15_all_chains <- data.frame(cbind(chain_1 = ht15_adj_phen_model_output_chain_1, chain_2 = ht15_adj_phen_model_output_chain_2$V2, 
                                     chain_3 = ht15_adj_phen_model_output_chain_3$V2, chain_4 = ht15_adj_phen_model_output_chain_4$V2, 
                                     chain_5 = ht15_adj_phen_model_output_chain_5$V2), row.names = 1)

ht15_all_chains <- ht15_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(ht15_all_chains)) %>% 
  column_to_rownames()

dbh15_all_chains <- data.frame(cbind(chain_1 = dbh15_adj_phen_model_output_chain_1, chain_2 = dbh15_adj_phen_model_output_chain_2$V2, 
                                    chain_3 = dbh15_adj_phen_model_output_chain_3$V2, chain_4 = dbh15_adj_phen_model_output_chain_4$V2, 
                                    chain_5 = dbh15_adj_phen_model_output_chain_5$V2), row.names = 1)

dbh15_all_chains <- dbh15_all_chains %>%
  rownames_to_column() %>% 
  mutate(mean = rowMeans(dbh15_all_chains)) %>% 
  column_to_rownames()


model_outputs <- data.frame(cbind(athuj_all_chains$mean, bthuj_all_chains$mean, sabinene_all_chains$mean, 
                                  athujaplicin_all_chains$mean, bthujaplicin_all_chains$mean,
                                  gthujaplicin_all_chains$mean, bthujaplicinol_all_chains$mean, ht15_all_chains$mean, dbh15_all_chains$mean))

model_outputs <- data.frame(t(model_outputs))
model_outputs <- data.frame(cbind(trait = c("a-thujone (2,206 SNPs)", "b-thujone (2,704 SNPs)", "sabinene (3,602 SNPs)", 
                                            "a-thujaplicin (2,508 SNPs)", "b-thujaplicin (2,559 SNPs)", "g-thujaplicin (1,787 SNPs)",
                                            "b-thujaplicinol (2,026 SNPs)", "ht15 (3,016 SNPs)", "dbh15 (3,505 SNPs)"), 
                                  V3 = (model_outputs$Vk2/model_outputs$Va), 
                                  V2 = (model_outputs$Vk3/model_outputs$Va), 
                                  V1 = (model_outputs$Vk4/model_outputs$Va)))

Vk <- melt(model_outputs, id.vars = 1)
Vk$value <- as.numeric(Vk$value)

Vk$variable <- factor(Vk$variable, levels = rev(levels(Vk$variable)))
Vk$trait <- factor(Vk$trait, levels = c("a-thujone (2,206 SNPs)", "b-thujone (2,704 SNPs)", "sabinene (3,602 SNPs)", 
                                        "a-thujaplicin (2,508 SNPs)", "b-thujaplicin (2,559 SNPs)", "g-thujaplicin (1,787 SNPs)",
                                        "b-thujaplicinol (2,026 SNPs)", "ht15 (3,016 SNPs)", "dbh15 (3,505 SNPs)"))


Vk$n_snps <- rep(c(2629, 2863, 3505, 2860, 2589, 2024, 2322, 3521, 3458), 3)

Vk$heritabilities <- c(heritabilities[3,2], heritabilities[4,2],heritabilities[5,2], heritabilities[6,2], heritabilities[7,2], heritabilities[8,2],
                       heritabilities[9,2], heritabilities[1,2], heritabilities[2,2])

V1 <- Vk %>% 
  filter(variable == "V1")

V2 <- Vk %>% 
  filter(variable == "V2")

V3 <- Vk %>%
  filter(variable == "V3")

# Some correlations
cor(V1$value, V1$n_snps, method = "spearman")
cor(V2$value, V2$n_snps)
cor(V3$value, V3$n_snps)

cor(V1$heritabilities, V1$n_snps)

cor(V1$value, V1$heritabilities)
cor(V2$value, V2$heritabilities)
cor(V3$value, V3$heritabilities)

col <- pal_nejm()(8)
grey <- scale_colour_grey()(3)

# Figure 1

(mixture <- ggplot(Vk, aes(x = trait, y = value, fill = variable)) +
    geom_col(aes(fill = variable)) +
    scale_y_continuous("Variance Explained", breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
    scale_x_discrete("Trait", breaks = unique(Vk$trait),
                     labels = c(expression(paste(alpha, "-thujone")), expression(paste(beta, "-thujone")), "Sabinene", 
                                expression(paste(alpha, "-thujaplicin")), expression(paste(beta, "-thujaplicin")), expression(paste(gamma, "-thujaplicin")),
                                expression(paste(beta, "-thujaplicinol")), "Height", "DBH")) +
    theme_Publication() +
    # scale_fill_nejm() +
    # scale_colour_manual(values = c("#1B1919FF", "#1B1919FF", "#1B1919FF")) +
    scale_fill_grey(name = "Mixture\nDistribution",
                      breaks = c("V1", "V2", "V3"),
                      labels = c(expression(10^-2%*%sigma[g]^2), 
                                 expression(10^-3%*%sigma[g]^2),
                                 expression(10^-4%*%sigma[g]^2))) +
    guides(colour = "none")) +
    theme(axis.text.x = element_text(size = rel(0.75)))

ggsave("bayesR_trait_mixture_components_grey.tiff", dpi = 300, width = 14, height = 8)

### Correlations for # snps ####
bayesr_table <- read.csv("bayesr_table.csv")
cor(bayesr_table$putative_causal_SNPs, bayesr_table$Annotated_genes)
cor(bayesr_table$putative_causal_SNPs, bayesr_table$Large_effect_SNPs)
