rm(list=ls())

library(gtsummary)
library(reshape2)
library(ggsci)
library(tidyverse)
source("~/UBC/GSAT/PhD/WRC/r_scripts/publication_theme.r")

library(extrafont)
library(extrafontdb)

setwd("E:/Tal/Documents/UBC/GSAT/PhD/WRC/papers/genetic_architecture/supplementary_code/snp_annotation/")

# Read in files

# Variant Effect Prediction (VEP)
vep <- read.delim("vep_training", header = T, stringsAsFactors = T, na.strings = "-")

# Outlier SNPs from Shalev et al. (2022)
outlier_snps <- read.table("outlier_snps_all_snps.txt")

# Data tables from UpsetR plots for extracting overlapping SNPs
traits_all_snps <- read.csv("upsetr_traits_all_snps.csv")
traits_top_snps <- read.csv("upsetr_traits_top_snps.csv")

# SNPs with effect for each trait (from BayesR)
a_thujone_snps_some_effect <- read.table("athujone_snps_with_effect.txt", header = T)
b_thujone_snps_some_effect <- read.table("bthujone_snps_with_effect.txt", header = T)
sabinene_snps_some_effect <- read.table("sabinene_snps_with_effect.txt", header = T)
a_thujaplicin_snps_some_effect <- read.table("athujaplicin_snps_with_effect.txt", header = T)
b_thujaplicin_snps_some_effect <- read.table("bthujaplicin_snps_with_effect.txt", header = T)
b_thujaplicinol_snps_some_effect <- read.table("gthujaplicin_snps_with_effect.txt", header = T)
g_thujaplicin_snps_some_effect <- read.table("bthujaplicinol_snps_with_effect.txt", header = T)
ht15_snps_some_effect <- read.table("ht15_snps_with_effect.txt", header = T)
dbh15_snps_some_effect <- read.table("dbh15_snps_with_effect.txt", header = T)

a_thujone_snps_top_effect <- read.table("athujone_snps_with_top_effect.txt", header = T)
b_thujone_snps_top_effect <- read.table("bthujone_snps_with_top_effect.txt", header = T)
sabinene_snps_top_effect <- read.table("sabinene_snps_with_top_effect.txt", header = T)
a_thujaplicin_snps_top_effect <- read.table("athujaplicin_snps_with_top_effect.txt", header = T)
b_thujaplicin_snps_top_effect <- read.table("bthujaplicin_snps_with_top_effect.txt", header = T)
b_thujaplicinol_snps_top_effect <- read.table("gthujaplicin_snps_with_top_effect.txt", header = T)
g_thujaplicin_snps_top_effect <- read.table("bthujaplicinol_snps_with_top_effect.txt", header = T)
ht15_snps_top_effect <- read.table("ht15_snps_with_top_effect.txt", header = T)
dbh15_snps_top_effect <- read.table("dbh15_snps_with_top_effect.txt", header = T)

a_thujone_snps_small_effect <- read.table("athujone_snps_with_small_effect.txt", header = T)
b_thujone_snps_small_effect <- read.table("bthujone_snps_with_small_effect.txt", header = T)
sabinene_snps_small_effect <- read.table("sabinene_snps_with_small_effect.txt", header = T)
a_thujaplicin_snps_small_effect <- read.table("athujaplicin_snps_with_small_effect.txt", header = T)
b_thujaplicin_snps_small_effect <- read.table("bthujaplicin_snps_with_small_effect.txt", header = T)
b_thujaplicinol_snps_small_effect <- read.table("gthujaplicin_snps_with_small_effect.txt", header = T)
g_thujaplicin_snps_small_effect <- read.table("bthujaplicinol_snps_with_small_effect.txt", header = T)
ht15_snps_small_effect <- read.table("ht15_snps_with_small_effect.txt", header = T)
dbh15_snps_small_effect <- read.table("dbh15_snps_with_small_effect.txt", header = T)

a_thujone_snps_poly_effect <- read.table("athujone_snps_with_poly_effect.txt", header = T)
b_thujone_snps_poly_effect <- read.table("bthujone_snps_with_poly_effect.txt", header = T)
sabinene_snps_poly_effect <- read.table("sabinene_snps_with_poly_effect.txt", header = T)
a_thujaplicin_snps_poly_effect <- read.table("athujaplicin_snps_with_poly_effect.txt", header = T)
b_thujaplicin_snps_poly_effect <- read.table("bthujaplicin_snps_with_poly_effect.txt", header = T)
b_thujaplicinol_snps_poly_effect <- read.table("gthujaplicin_snps_with_poly_effect.txt", header = T)
g_thujaplicin_snps_poly_effect <- read.table("bthujaplicinol_snps_with_poly_effect.txt", header = T)
ht15_snps_poly_effect <- read.table("ht15_snps_with_poly_effect.txt", header = T)
dbh15_snps_poly_effect <- read.table("dbh15_snps_with_poly_effect.txt", header = T)

# top effect intersections
traits_top_a_b_thujone <- traits_top_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1)

traits_top_a_b_thujone_sabinene <- traits_top_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1 & sabinene == 1)

traits_top_a_b_thujone_not_sabinene <- traits_top_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1 & sabinene == 0)

traits_top_ht_dbh <- traits_top_snps %>% 
  filter(ht15 == 1 & dbh15 == 1)

traits_top_a_thujone_sabinene <- traits_top_snps %>% 
  filter(a.thujone == 1 & sabinene == 1)

traits_top_a_b_thujaplicin <- traits_top_snps %>% 
  filter(a.thujaplicin == 1 & b.thujaplicin == 1 & b.thujaplicinol == 0)

traits_top_a_g_thujaplicin_b_thujaplicinol <- traits_top_snps %>% 
  filter(a.thujaplicin == 1 & g.thujaplicin == 1 & b.thujaplicinol == 1)

traits_top_g_thujaplicin_b_thujaplicinol <- traits_top_snps %>% 
  filter(a.thujaplicin == 0 & g.thujaplicin == 1 & b.thujaplicinol == 1 & b.thujaplicin == 0)

# all SNPs unique to specific intersections
traits_a_b_thujone_unique <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_a_b_thujone_sabinene_unique <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1 & sabinene == 1 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_b_thujone_sabinene_unique <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 1 & sabinene == 1 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_ht_dbh_unique <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 1 & dbh15 == 1)

traits_a_thujone_sabinene_unique <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 0 & sabinene == 1 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_a_b_thujaplicin_unique <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 1 & b.thujaplicin == 1 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_a_g_thujaplicin_b_thujaplicinol_unique <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 1 & b.thujaplicin == 0 & g.thujaplicin == 1 & b.thujaplicinol == 1 & ht15 == 0 & dbh15 == 0)

traits_g_thujaplicin_b_thujaplicinol_unique <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 1 & b.thujaplicinol == 1 & ht15 == 0 & dbh15 == 0)

## all SNPs unique to each trait
# all SNPs
traits_a_thujone <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_b_thujone <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 1 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_sabinene <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 1 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_a_thujaplicin <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 1 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_b_thujaplicin <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 1 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_g_thujaplicin <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 1 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 0)

traits_b_thujaplicinol <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 1 & ht15 == 0 & dbh15 == 0)

traits_ht <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 1 & dbh15 == 0)

traits_dbh <- traits_all_snps %>% 
  filter(a.thujone == 0 & b.thujone == 0 & sabinene == 0 & a.thujaplicin == 0 & b.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & ht15 == 0 & dbh15 == 1)

# all SNPs non-unique
traits_a_b_thujone <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1)

traits_a_b_thujone_sabinene <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1 & sabinene == 1)

traits_b_thujone_sabinene <- traits_all_snps %>% 
  filter(b.thujone == 1 & sabinene == 1)

traits_a_b_thujone_not_sabinene <- traits_all_snps %>% 
  filter(a.thujone == 1 & b.thujone == 1 & sabinene == 0)

traits_ht_dbh <- traits_all_snps %>% 
  filter(ht15 == 1 & dbh15 == 1)

traits_a_thujone_sabinene <- traits_all_snps %>% 
  filter(a.thujone == 1 & sabinene == 1)

traits_a_b_thujaplicin <- traits_all_snps %>% 
  filter(a.thujaplicin == 1 & b.thujaplicin == 1)

traits_a_g_thujaplicin_b_thujaplicinol <- traits_all_snps %>% 
  filter(a.thujaplicin == 1 & g.thujaplicin == 1 & b.thujaplicinol == 1 & b.thujaplicin == 0)

traits_g_thujaplicin_b_thujaplicinol <- traits_all_snps %>% 
  filter(g.thujaplicin == 1 & b.thujaplicinol == 1)

traits_all_wood <- traits_all_snps %>% 
  filter(a.thujaplicin == 1 & b.thujaplicin == 1 & g.thujaplicin == 1 & b.thujaplicinol == 1)

traits_b_thujaplicin <- traits_all_snps %>% 
  filter(a.thujaplicin == 0 & g.thujaplicin == 0 & b.thujaplicinol == 0 & b.thujaplicin == 1)

### VEP for all SNPs /each trait
a_thujone_vep <- vep %>% 
  filter(SNP %in% a_thujone_snps_some_effect$SNP)

b_thujone_vep <- vep %>% 
  filter(SNP %in% b_thujone_snps_some_effect$SNP)

sabinene_vep <- vep %>% 
  filter(SNP %in% sabinene_snps_some_effect$SNP)

a_thujaplicin_vep <- vep %>% 
  filter(SNP %in% a_thujaplicin_snps_some_effect$SNP)

b_thujaplicin_vep <- vep %>% 
  filter(SNP %in% b_thujaplicin_snps_some_effect$SNP)

g_thujaplicin_vep <- vep %>% 
  filter(SNP %in% g_thujaplicin_snps_some_effect$SNP)

b_thujaplicinol_vep <- vep %>% 
  filter(SNP %in% b_thujaplicinol_snps_some_effect$SNP)

ht15_vep <- vep %>% 
  filter(SNP %in% ht15_snps_some_effect$SNP)

dbh15_vep <- vep %>% 
  filter(SNP %in% dbh15_snps_some_effect$SNP)

### VEP for top SNPs/each trait
a_thujone_top_vep <- vep %>% 
  filter(SNP %in% a_thujone_snps_top_effect$SNP)

b_thujone_top_vep <- vep %>% 
  filter(SNP %in% b_thujone_snps_top_effect$SNP)

sabinene_top_vep <- vep %>% 
  filter(SNP %in% sabinene_snps_top_effect$SNP)

a_thujaplicin_top_vep <- vep %>% 
  filter(SNP %in% a_thujaplicin_snps_top_effect$SNP)

b_thujaplicin_top_vep <- vep %>% 
  filter(SNP %in% b_thujaplicin_snps_top_effect$SNP)

g_thujaplicin_top_vep <- vep %>% 
  filter(SNP %in% g_thujaplicin_snps_top_effect$SNP)

b_thujaplicinol_top_vep <- vep %>% 
  filter(SNP %in% b_thujaplicinol_snps_top_effect$SNP)

ht15_top_vep <- vep %>% 
  filter(SNP %in% ht15_snps_top_effect$SNP)

dbh15_top_vep <- vep %>% 
  filter(SNP %in% dbh15_snps_top_effect$SNP)

### traits VEP for all SNPs unique to each trait
a_thujone_traits_vep <- vep %>% 
  filter(SNP %in%  traits_a_thujone$SNP)

b_thujone_traits_vep <- vep %>% 
  filter(SNP %in%  traits_b_thujone$SNP)

sabinene_traits_vep <- vep %>% 
  filter(SNP %in%  traits_sabinene$SNP)

a_thujaplicin_traits_vep <- vep %>% 
  filter(SNP %in%  traits_a_thujaplicin$SNP)

b_thujaplicin_traits_vep <- vep %>% 
  filter(SNP %in%  traits_b_thujaplicin$SNP)

g_thujaplicin_traits_vep <- vep %>% 
  filter(SNP %in%  traits_g_thujaplicin$SNP)

b_thujaplicinol_traits_vep <- vep %>% 
  filter(SNP %in%  traits_b_thujaplicinol$SNP)

ht15_traits_vep <- vep %>% 
  filter(SNP %in%  traits_ht$SNP)

dbh15_traits_vep <- vep %>% 
  filter(SNP %in%  traits_dbh$SNP)


### VEP for top interacting SNPs ####
a_b_thujone_top_vep <- vep %>% 
  filter(SNP %in% traits_top_a_b_thujone$SNP)

a_b_thujone_not_sabinene_top_vep <- vep %>% 
  filter(SNP %in% traits_top_a_b_thujone_not_sabinene$SNP)

a_b_thujone_sabinene_top_vep <- vep %>% 
  filter(SNP %in% traits_top_a_b_thujone_sabinene$SNP)

ht_dbh_top_vep <- vep %>% 
  filter(SNP %in% traits_top_ht_dbh$SNP)

a_thujone_sabinene_top_vep <- vep %>% 
  filter(SNP %in% traits_top_a_thujone_sabinene$SNP)

a_b_thujaplicin_top_vep <- vep %>% 
  filter(SNP %in% traits_top_a_b_thujaplicin$SNP)

a_g_thujaplicin_b_thujaplicinol_top_vep <- vep %>% 
  filter(SNP %in% traits_top_a_g_thujaplicin_b_thujaplicinol$SNP)

g_thujaplicin_b_thujaplicinol_top_vep <- vep %>% 
  filter(SNP %in% traits_top_g_thujaplicin_b_thujaplicinol$SNP)

### VEP for all interacting SNPs ####
a_b_thujone_vep <- vep %>% 
  filter(SNP %in% traits_a_b_thujone$SNP)

a_b_thujone_not_sabinene_vep <- vep %>% 
  filter(SNP %in% traits_a_b_thujone_not_sabinene$SNP)

a_b_thujone_sabinene_vep <- vep %>% 
  filter(SNP %in% traits_a_b_thujone_sabinene$SNP)

ht_dbh_vep <- vep %>% 
  filter(SNP %in% traits_ht_dbh$SNP)

a_thujone_sabinene_vep <- vep %>% 
  filter(SNP %in% traits_a_thujone_sabinene$SNP)

a_b_thujaplicin_vep <- vep %>% 
  filter(SNP %in% traits_a_b_thujaplicin$SNP)

a_g_thujaplicin_b_thujaplicinol_vep <- vep %>% 
  filter(SNP %in% traits_a_g_thujaplicin_b_thujaplicinol$SNP)

g_thujaplicin_b_thujaplicinol_vep <- vep %>% 
  filter(SNP %in% traits_g_thujaplicin_b_thujaplicinol$SNP)

all_wood_vep <- vep %>% 
  filter(SNP %in% traits_all_wood$SNP)


### JGI annotations from genome ##########################

gene_annotations <- read.delim("gene_functions.txt", header = F, stringsAsFactors = F, sep = "\t")
gff <- read.table("Tplicatav3.1.primaryTrs.gff3")

# GFF with simplified gene names
simplified_gff <- read.table("Tplicatav3.1.primaryTrs_simplified_genes.gff3")

gff <- gff %>% 
  mutate(V10 = simplified_gff$V9)

gff_genes <- gff %>% 
  filter(V3 == "gene")

all_snp_annotations <- gene_annotations %>%
  filter(V1 %in% vep$Gene) %>% 
  rename(Gene = V1)

all_snp_annotations <- merge(vep, all_snp_annotations, by = "Gene", all = T) 

# write.csv(all_snp_annotations, "all_snp_annotations.csv")

## Genes for each trait ####
a_thujone_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% a_thujone_vep$Gene) %>% 
  rename(Gene = V1)

a_thujone_snp_annotations <- merge(a_thujone_vep, a_thujone_snp_annotations, by = "Gene", all = T)

a_thujone_snp_annotations_only_genes <- a_thujone_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_thujone_snp_annotations_only_genes$Gene))
summary(a_thujone_snp_annotations$Consequence)

# write.csv(a_thujone_snp_annotations_only_genes, "a_thujone_snp_annotations_only_genes.csv")

b_thujone_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% b_thujone_vep$Gene) %>% 
  rename(Gene = V1)

b_thujone_snp_annotations <- merge(b_thujone_vep, b_thujone_snp_annotations, by = "Gene", all = T)

b_thujone_snp_annotations_only_genes <- b_thujone_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(b_thujone_snp_annotations_only_genes$Gene))
summary(b_thujone_snp_annotations$Consequence)

# write.csv(b_thujone_snp_annotations_only_genes, "b_thujone_snp_annotations_only_genes.csv")

sabinene_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% sabinene_vep$Gene) %>% 
  rename(Gene = V1)

sabinene_snp_annotations <- merge(sabinene_vep, sabinene_snp_annotations, by = "Gene", all = T)

sabinene_snp_annotations_only_genes <- sabinene_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(sabinene_snp_annotations_only_genes$Gene))
summary(sabinene_snp_annotations$Consequence)

# write.csv(sabinene_snp_annotations_only_genes, "sabinene_snp_annotations_only_genes.csv")

a_thujaplicin_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% a_thujaplicin_vep$Gene) %>% 
  rename(Gene = V1)

a_thujaplicin_snp_annotations <- merge(a_thujaplicin_vep, a_thujaplicin_snp_annotations, by = "Gene", all = T)

a_thujaplicin_snp_annotations_only_genes <- a_thujaplicin_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_thujaplicin_snp_annotations_only_genes$Gene))
summary(a_thujaplicin_snp_annotations$Consequence)

b_thujaplicin_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% b_thujaplicin_vep$Gene) %>% 
  rename(Gene = V1)

b_thujaplicin_snp_annotations <- merge(b_thujaplicin_vep, b_thujaplicin_snp_annotations, by = "Gene", all = T)

b_thujaplicin_snp_annotations_only_genes <- b_thujaplicin_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(b_thujaplicin_snp_annotations_only_genes$Gene))
summary(b_thujaplicin_snp_annotations$Consequence)

g_thujaplicin_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% g_thujaplicin_vep$Gene) %>% 
  rename(Gene = V1)

g_thujaplicin_snp_annotations <- merge(g_thujaplicin_vep, g_thujaplicin_snp_annotations, by = "Gene", all = T)

g_thujaplicin_snp_annotations_only_genes <- g_thujaplicin_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(g_thujaplicin_snp_annotations_only_genes$Gene))
summary(g_thujaplicin_snp_annotations$Consequence)

b_thujaplicinol_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% b_thujaplicinol_vep$Gene) %>% 
  rename(Gene = V1)

b_thujaplicinol_snp_annotations <- merge(b_thujaplicinol_vep, b_thujaplicinol_snp_annotations, by = "Gene", all = T)

b_thujaplicinol_snp_annotations_only_genes <- b_thujaplicinol_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(b_thujaplicinol_snp_annotations_only_genes$Gene))
summary(b_thujaplicinol_snp_annotations$Consequence)

ht15_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% ht15_vep$Gene) %>% 
  rename(Gene = V1)

ht15_snp_annotations <- merge(ht15_vep, ht15_snp_annotations, by = "Gene", all = T)

ht15_snp_annotations_only_genes <- ht15_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(ht15_snp_annotations_only_genes$Gene))
summary(ht15_snp_annotations$Consequence)

dbh15_snp_annotations <- gene_annotations %>% 
  filter(V1 %in% dbh15_vep$Gene) %>% 
  rename(Gene = V1)

dbh15_snp_annotations <- merge(dbh15_vep, dbh15_snp_annotations, by = "Gene", all = T)

dbh15_snp_annotations_only_genes <- dbh15_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(dbh15_snp_annotations_only_genes$Gene))
levels(factor(dbh15_snp_annotations_only_genes$SNP))

all_snp_vep <- rbind(a_thujone_vep, b_thujone_vep, sabinene_vep, a_thujaplicin_vep, b_thujaplicin_vep, g_thujaplicin_vep, b_thujaplicinol_vep, ht15_vep, dbh15_vep)
# write.csv(all_snp_vep, "all_snp_vep.csv")

## Genes for snps unique to each trait ####
a_thujone_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% a_thujone_traits_vep$Gene) %>% 
  rename(Gene = V1)

a_thujone_snp_traits_annotations <- merge(a_thujone_traits_vep, a_thujone_snp_traits_annotations, by = "Gene", all = T)

a_thujone_snp_traits_annotations_only_genes <- a_thujone_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_thujone_snp_traits_annotations_only_genes$Gene))
summary(a_thujone_snp_traits_annotations$Consequence)

# write.csv(a_thujone_snp_traits_annotations_only_genes, "a_thujone_unique_snps_annotations_only_genes.csv")

b_thujone_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% b_thujone_traits_vep$Gene) %>% 
  rename(Gene = V1)

b_thujone_snp_traits_annotations <- merge(b_thujone_traits_vep, b_thujone_snp_traits_annotations, by = "Gene", all = T)

b_thujone_snp_traits_annotations_only_genes <- b_thujone_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(b_thujone_snp_traits_annotations_only_genes$Gene))
summary(b_thujone_snp_traits_annotations$Consequence)

# write.csv(b_thujone_snp_traits_annotations_only_genes, "b_thujone_unique_snps_annotations_only_genes.csv")

sabinene_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% sabinene_traits_vep$Gene) %>% 
  rename(Gene = V1)

sabinene_snp_traits_annotations <- merge(sabinene_traits_vep, sabinene_snp_traits_annotations, by = "Gene", all = T)

sabinene_snp_traits_annotations_only_genes <- sabinene_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(sabinene_snp_traits_annotations_only_genes$Gene))
summary(sabinene_snp_traits_annotations$Consequence)

# write.csv(sabinene_snp_traits_annotations_only_genes, "sabinene_unique_snps_annotations_only_genes.csv")

a_thujaplicin_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% a_thujaplicin_traits_vep$Gene) %>% 
  rename(Gene = V1)

a_thujaplicin_snp_traits_annotations <- merge(a_thujaplicin_traits_vep, a_thujaplicin_snp_traits_annotations, by = "Gene", all = T)

a_thujaplicin_snp_traits_annotations_only_genes <- a_thujaplicin_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_thujaplicin_snp_traits_annotations_only_genes$Gene))
summary(a_thujaplicin_snp_traits_annotations$Consequence)

b_thujaplicin_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% b_thujaplicin_traits_vep$Gene) %>% 
  rename(Gene = V1)

b_thujaplicin_snp_traits_annotations <- merge(b_thujaplicin_traits_vep, b_thujaplicin_snp_traits_annotations, by = "Gene", all = T)

b_thujaplicin_snp_traits_annotations_only_genes <- b_thujaplicin_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(b_thujaplicin_snp_traits_annotations_only_genes$Gene))
summary(b_thujaplicin_snp_traits_annotations$Consequence)

g_thujaplicin_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% g_thujaplicin_traits_vep$Gene) %>% 
  rename(Gene = V1)

g_thujaplicin_snp_traits_annotations <- merge(g_thujaplicin_traits_vep, g_thujaplicin_snp_traits_annotations, by = "Gene", all = T)

g_thujaplicin_snp_traits_annotations_only_genes <- g_thujaplicin_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(g_thujaplicin_snp_traits_annotations_only_genes$Gene))
summary(g_thujaplicin_snp_traits_annotations$Consequence)

b_thujaplicinol_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% b_thujaplicinol_traits_vep$Gene) %>% 
  rename(Gene = V1)

b_thujaplicinol_snp_traits_annotations <- merge(b_thujaplicinol_traits_vep, b_thujaplicinol_snp_traits_annotations, by = "Gene", all = T)

b_thujaplicinol_snp_traits_annotations_only_genes <- b_thujaplicinol_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(b_thujaplicinol_snp_traits_annotations_only_genes$Gene))
summary(b_thujaplicinol_snp_traits_annotations$Consequence)

ht15_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% ht15_traits_vep$Gene) %>% 
  rename(Gene = V1)

ht15_snp_traits_annotations <- merge(ht15_traits_vep, ht15_snp_traits_annotations, by = "Gene", all = T)

ht15_snp_traits_annotations_only_genes <- ht15_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(ht15_snp_traits_annotations_only_genes$Gene))
levels(factor(ht15_snp_traits_annotations_only_genes$SNP))

dbh15_snp_traits_annotations <- gene_annotations %>% 
  filter(V1 %in% dbh15_traits_vep$Gene) %>% 
  rename(Gene = V1)

dbh15_snp_traits_annotations <- merge(dbh15_traits_vep, dbh15_snp_traits_annotations, by = "Gene", all = T)

dbh15_snp_traits_annotations_only_genes <- dbh15_snp_traits_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(dbh15_snp_traits_annotations_only_genes$Gene))
levels(factor(dbh15_snp_traits_annotations_only_genes$SNP))

## annotations for all interacting snps
## a-b-thujone snps ###
a_b_thujone_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujone_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujone_snp_annotations <- merge(a_b_thujone_vep, a_b_thujone_snp_annotations, by = "Gene", all = T)  

a_b_thujone_snp_annotations_only_genes <- a_b_thujone_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujone_snp_annotations_only_genes$Gene))
levels(factor(a_b_thujone_snp_annotations_only_genes$SNP))

## a-b-thujone-not sabinene snps ###
a_b_thujone_not_sabinene_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujone_not_sabinene_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujone_not_sabinene_snp_annotations <- merge(a_b_thujone_not_sabinene_vep, a_b_thujone_not_sabinene_snp_annotations, by = "Gene", all = T)  

a_b_thujone_not_sabinene_snp_annotations_only_genes <- a_b_thujone_not_sabinene_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujone_not_sabinene_snp_annotations_only_genes$Gene))
levels(factor(a_b_thujone_not_sabinene_snp_annotations_only_genes$SNP))
summary(a_b_thujone_sabinene_snp_annotations$Consequence)

## a-b-thujone-sabinene snps ###
a_b_thujone_sabinene_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujone_sabinene_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujone_sabinene_snp_annotations <- merge(a_b_thujone_sabinene_vep, a_b_thujone_sabinene_snp_annotations, by = "Gene", all = T)  

a_b_thujone_sabinene_snp_annotations_only_genes <- a_b_thujone_sabinene_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujone_sabinene_snp_annotations_only_genes$Gene))
levels(factor(a_b_thujone_sabinene_snp_annotations_only_genes$SNP))
summary(a_b_thujone_sabinene_snp_annotations$Consequence)

## ht -dbh snps
ht_dbh_snp_annotations <- gene_annotations %>%
  filter(V1 %in% ht_dbh_vep$Gene) %>% 
  rename(Gene = V1)

ht_dbh_snp_annotations <- merge(ht_dbh_vep, ht_dbh_snp_annotations, by = "Gene", all = T) 

ht_dbh_snp_annotations_only_genes <- ht_dbh_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(ht_dbh_snp_annotations_only_genes$Gene))
levels(factor(ht_dbh_snp_annotations_only_genes$SNP))

## a-thujone - sabinene snps
a_thujone_sabinene_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_thujone_sabinene_vep$Gene) %>% 
  rename(Gene = V1)

a_thujone_sabinene_snp_annotations <- merge(a_thujone_sabinene_vep, a_thujone_sabinene_snp_annotations, by = "Gene", all = T) 

a_thujone_sabinene_snp_annotations_only_genes <- a_thujone_sabinene_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_thujone_sabinene_snp_annotations_only_genes$Gene))
levels(factor(a_thujone_sabinene_snp_annotations$SNP))
summary(a_thujone_sabinene_snp_annotations$Consequence)

## a-b-thujaplicin snps
a_b_thujaplicin_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujaplicin_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujaplicin_snp_annotations <- merge(a_b_thujaplicin_vep, a_b_thujaplicin_snp_annotations, by = "Gene", all = T) 

a_b_thujaplicin_snp_annotations_only_genes <- a_b_thujaplicin_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujaplicin_snp_annotations_only_genes$Gene))
levels(factor(a_b_thujaplicin_snp_annotations_only_genes$SNP))
summary(a_b_thujaplicin_snp_annotations$Consequence)

## a-g-b-thujaplicinol snps
a_g_thujaplicin_b_thujaplicinol_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_g_thujaplicin_b_thujaplicinol_vep$Gene) %>% 
  rename(Gene = V1)

a_g_thujaplicin_b_thujaplicinol_snp_annotations <- merge(a_g_thujaplicin_b_thujaplicinol_vep, a_g_thujaplicin_b_thujaplicinol_snp_annotations, by = "Gene", all = T) 

a_g_thujaplicin_b_thujaplicinol_snp_annotations_only_genes <- a_g_thujaplicin_b_thujaplicinol_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_g_thujaplicin_b_thujaplicinol_snp_annotations_only_genes$Gene))
summary(a_g_thujaplicin_b_thujaplicinol_snp_annotations$Consequence)

## g-b-thujaplicinol snps
g_thujaplicin_b_thujaplicinol_snp_annotations <- gene_annotations %>%
  filter(V1 %in% g_thujaplicin_b_thujaplicinol_vep$Gene) %>% 
  rename(Gene = V1)

g_thujaplicin_b_thujaplicinol_snp_annotations <- merge(g_thujaplicin_b_thujaplicinol_vep, g_thujaplicin_b_thujaplicinol_snp_annotations, by = "Gene", all = T) 

g_thujaplicin_b_thujaplicinol_snp_annotations_only_genes <- g_thujaplicin_b_thujaplicinol_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(g_thujaplicin_b_thujaplicinol_snp_annotations_only_genes$Gene))
levels(factor(g_thujaplicin_b_thujaplicinol_snp_annotations_only_genes$SNP))

## annotations for top interacting snps
## a-b-thujone snps ###
a_b_thujone_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujone_top_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujone_top_snp_annotations <- merge(a_b_thujone_top_vep, a_b_thujone_top_snp_annotations, by = "Gene", all = T)  

a_b_thujone_top_snp_annotations_only_genes <- a_b_thujone_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujone_top_snp_annotations_only_genes$Gene))
levels(factor(a_b_thujone_top_snp_annotations_only_genes$SNP))
summary(a_b_thujone_top_snp_annotations$Consequence)

## a-b-thujone-sabinene snps ###
a_b_thujone_sabinene_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujone_sabinene_top_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujone_sabinene_top_snp_annotations <- merge(a_b_thujone_sabinene_top_vep, a_b_thujone_sabinene_top_snp_annotations, by = "Gene", all = T)  

a_b_thujone_sabinene_top_snp_annotations_only_genes <- a_b_thujone_sabinene_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujone_sabinene_top_snp_annotations_only_genes$Gene))
summary(a_b_thujone_sabinene_top_snp_annotations$Consequence)

## ht -dbh snps
ht_dbh_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% ht_dbh_top_vep$Gene) %>% 
  rename(Gene = V1)

ht_dbh_top_snp_annotations <- merge(ht_dbh_top_vep, ht_dbh_top_snp_annotations, by = "Gene", all = T) 

ht_dbh_top_snp_annotations_only_genes <- ht_dbh_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(ht_dbh_top_snp_annotations_only_genes$Gene))
summary(ht_dbh_top_snp_annotations$Consequence)

## a-thujone - sabinene snps
a_thujone_sabinene_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_thujone_sabinene_top_vep$Gene) %>% 
  rename(Gene = V1)

a_thujone_sabinene_top_snp_annotations <- merge(a_thujone_sabinene_top_vep, a_thujone_sabinene_top_snp_annotations, by = "Gene", all = T) 

a_thujone_sabinene_top_snp_annotations_only_genes <- a_thujone_sabinene_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_thujone_sabinene_top_snp_annotations_only_genes$Gene))
summary(a_thujone_sabinene_top_snp_annotations$Consequence)

## a-b-thujaplicin snps
a_b_thujaplicin_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_b_thujaplicin_top_vep$Gene) %>% 
  rename(Gene = V1)

a_b_thujaplicin_top_snp_annotations <- merge(a_b_thujaplicin_top_vep, a_b_thujaplicin_top_snp_annotations, by = "Gene", all = T) 

a_b_thujaplicin_top_snp_annotations_only_genes <- a_b_thujaplicin_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_b_thujaplicin_top_snp_annotations_only_genes$Gene))
summary(a_b_thujaplicin_top_snp_annotations$Consequence)

## a-g-b-thujaplicinol snps
a_g_thujaplicin_b_thujaplicinol_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% a_g_thujaplicin_b_thujaplicinol_top_vep$Gene) %>% 
  rename(Gene = V1)

a_g_thujaplicin_b_thujaplicinol_top_snp_annotations <- merge(a_g_thujaplicin_b_thujaplicinol_top_vep, a_g_thujaplicin_b_thujaplicinol_top_snp_annotations, by = "Gene", all = T) 

a_g_thujaplicin_b_thujaplicinol_top_snp_annotations_only_genes <- a_g_thujaplicin_b_thujaplicinol_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(a_g_thujaplicin_b_thujaplicinol_top_snp_annotations_only_genes$Gene))
summary(a_g_thujaplicin_b_thujaplicinol_top_snp_annotations$Consequence)

## g-b-thujaplicinol snps
g_thujaplicin_b_thujaplicinol_top_snp_annotations <- gene_annotations %>%
  filter(V1 %in% g_thujaplicin_b_thujaplicinol_top_vep$Gene) %>% 
  rename(Gene = V1)

g_thujaplicin_b_thujaplicinol_top_snp_annotations <- merge(g_thujaplicin_b_thujaplicinol_top_vep, g_thujaplicin_b_thujaplicinol_top_snp_annotations, by = "Gene", all = T) 

g_thujaplicin_b_thujaplicinol_top_snp_annotations_only_genes <- g_thujaplicin_b_thujaplicinol_top_snp_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant"))

levels(factor(g_thujaplicin_b_thujaplicinol_top_snp_annotations_only_genes$Gene))
summary(g_thujaplicin_b_thujaplicinol_top_snp_annotations$Consequence)

### Outlier SNPs ####
# All outlier SNP annotations
outlier_annotations <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_genes <- outlier_annotations %>% 
  filter(!(Consequence == "downstream_gene_variant" | Consequence == "upstream_gene_variant" | Consequence == "intergenic_variant")) %>% 
  droplevels()

table(outlier_genes$Gene)

unique(outlier_genes$Gene)

### Outlier SNPs with some effect #
all_outlier_snps <- outlier_snps %>% 
  filter(V1 %in% a_thujone_snps_some_effect$SNP | V1 %in% b_thujone_snps_some_effect$SNP | V1 %in% sabinene_snps_some_effect$SNP | 
           V1 %in% a_thujaplicin_snps_some_effect$SNP | V1 %in% b_thujaplicin_snps_some_effect$SNP | V1 %in% g_thujaplicin_snps_some_effect$SNP | 
           V1 %in% b_thujaplicinol_snps_some_effect$SNP | V1 %in% ht15_snps_some_effect$SNP | V1 %in% dbh15_snps_some_effect$SNP)

all_outlier_snps$V1

# Foliar SNPs
outlier_athujone_snps <- a_thujone_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_bthujone_snps <- b_thujone_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_sabinene_snps <- sabinene_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_foliar_snps <- as.vector(c(outlier_athujone_snps$SNP,outlier_bthujone_snps$SNP, outlier_sabinene_snps$SNP))
unique(outlier_foliar_snps)

# Wood SNPs
outlier_athujaplicin_snps <- a_thujaplicin_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicin_snps <-  b_thujaplicin_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_gthujaplicin_snps <- g_thujaplicin_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicinol_snps <-  b_thujaplicinol_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_wood_snps <- as.vector(c(outlier_athujaplicin_snps$SNP,outlier_bthujaplicin_snps$SNP, outlier_gthujaplicin_snps$SNP, outlier_bthujaplicinol_snps$SNP))
unique(outlier_wood_snps)

# Growth SNPs
outlier_ht15_snps <- ht15_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_dbh15_snps <- dbh15_snps_some_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_growth_snps <- as.vector(c(outlier_ht15_snps$SNP,outlier_dbh15_snps$SNP))
unique(outlier_growth_snps)

### Outlier SNPs with top effect #

top_outlier_snps <- outlier_snps %>% 
  filter(V1 %in% a_thujone_snps_top_effect$SNP | V1 %in% b_thujone_snps_top_effect$SNP | V1 %in% sabinene_snps_top_effect$SNP | 
           V1 %in% a_thujaplicin_snps_top_effect$SNP | V1 %in% b_thujaplicin_snps_top_effect$SNP | V1 %in% g_thujaplicin_snps_top_effect$SNP | 
           V1 %in% b_thujaplicinol_snps_top_effect$SNP | V1 %in% ht15_snps_top_effect$SNP | V1 %in% dbh15_snps_top_effect$SNP)

top_outlier_snps$V1

# Foliar SNPs
outlier_athujone_top_snps <- a_thujone_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_bthujone_top_snps <- b_thujone_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_sabinene_top_snps <- sabinene_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_foliar_top_snps <- as.vector(c(outlier_athujone_top_snps$SNP,outlier_bthujone_top_snps$SNP, outlier_sabinene_top_snps$SNP))
unique(outlier_foliar_top_snps)

# Wood SNPs
outlier_athujaplicin_top_snps <- a_thujaplicin_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicin_top_snps <-  b_thujaplicin_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_gthujaplicin_top_snps <- g_thujaplicin_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicinol_top_snps <-  b_thujaplicinol_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_wood_top_snps <- as.vector(c(outlier_athujaplicin_top_snps$SNP,outlier_bthujaplicin_top_snps$SNP, outlier_gthujaplicin_top_snps$SNP, outlier_bthujaplicinol_top_snps$SNP))
unique(outlier_wood_top_snps)

# Growth SNPs
outlier_ht15_top_snps <- ht15_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_dbh15_top_snps <- dbh15_snps_top_effect %>% 
  filter(SNP %in% outlier_snps$V1)

outlier_growth_top_snps <- as.vector(c(outlier_ht15_top_snps$SNP,outlier_dbh15_top_snps$SNP))
unique(outlier_growth_top_snps)

# Outlier SNPs with small effect

small_outlier_snps <- outlier_snps %>% 
  filter(V1 %in% a_thujone_snps_small_effect$SNP | V1 %in% b_thujone_snps_small_effect$SNP | V1 %in% sabinene_snps_small_effect$SNP | 
           V1 %in% a_thujaplicin_snps_small_effect$SNP | V1 %in% b_thujaplicin_snps_small_effect$SNP | V1 %in% g_thujaplicin_snps_small_effect$SNP | 
           V1 %in% b_thujaplicinol_snps_small_effect$SNP | V1 %in% ht15_snps_small_effect$SNP | V1 %in% dbh15_snps_small_effect$SNP)

small_outlier_snps$V1

# Foliar SNPs
outlier_athujone_small_snps <- a_thujone_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_bthujone_small_snps <- b_thujone_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_sabinene_small_snps <- sabinene_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_foliar_small_snps <- as.vector(c(outlier_athujone_small_snps$SNP,outlier_bthujone_small_snps$SNP, outlier_sabinene_small_snps$SNP))
unique(outlier_foliar_small_snps)

# wood SNPs
outlier_athujaplicin_small_snps <- a_thujaplicin_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicin_small_snps <-  b_thujaplicin_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_gthujaplicin_small_snps <- g_thujaplicin_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicinol_small_snps <-  b_thujaplicinol_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_wood_small_snps <- as.vector(c(outlier_athujaplicin_small_snps$SNP,outlier_bthujaplicin_small_snps$SNP, outlier_gthujaplicin_small_snps$SNP, outlier_bthujaplicinol_small_snps$SNP))
unique(outlier_wood_small_snps)

# growth SNPs
outlier_ht15_small_snps <- ht15_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_dbh15_small_snps <- dbh15_snps_small_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_growth_small_snps <- as.vector(c(outlier_ht15_small_snps$SNP,outlier_dbh15_small_snps$SNP))
unique(outlier_growth_small_snps)

# Outlier SNPs with poly effect

poly_outlier_snps <- outlier_snps %>% 
  filter(V1 %in% a_thujone_snps_poly_effect$SNP | V1 %in% b_thujone_snps_poly_effect$SNP | V1 %in% sabinene_snps_poly_effect$SNP | 
           V1 %in% a_thujaplicin_snps_poly_effect$SNP | V1 %in% b_thujaplicin_snps_poly_effect$SNP | V1 %in% g_thujaplicin_snps_poly_effect$SNP | 
           V1 %in% b_thujaplicinol_snps_poly_effect$SNP | V1 %in% ht15_snps_poly_effect$SNP | V1 %in% dbh15_snps_poly_effect$SNP)

poly_outlier_snps$V1

# Foliar SNPs
outlier_athujone_poly_snps <- a_thujone_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_bthujone_poly_snps <- b_thujone_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_sabinene_poly_snps <- sabinene_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_foliar_poly_snps <- as.vector(c(outlier_athujone_poly_snps$SNP,outlier_bthujone_poly_snps$SNP, outlier_sabinene_poly_snps$SNP))
unique(outlier_foliar_poly_snps)

# wood SNPs
outlier_athujaplicin_poly_snps <- a_thujaplicin_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicin_poly_snps <-  b_thujaplicin_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_gthujaplicin_poly_snps <- g_thujaplicin_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_bthujaplicinol_poly_snps <-  b_thujaplicinol_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_wood_poly_snps <- as.vector(c(outlier_athujaplicin_poly_snps$SNP,outlier_bthujaplicin_poly_snps$SNP, outlier_gthujaplicin_poly_snps$SNP, outlier_bthujaplicinol_poly_snps$SNP))
unique(outlier_wood_poly_snps)

# Growth SNPs
outlier_ht15_poly_snps <- ht15_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_dbh15_poly_snps <- dbh15_snps_poly_effect %>%
  filter(SNP %in% outlier_snps$V1)

outlier_growth_poly_snps <- as.vector(c(outlier_ht15_poly_snps$SNP,outlier_dbh15_poly_snps$SNP))
unique(outlier_growth_poly_snps)

### Outlier annotation with some effect ##
outlier_trait_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & (SNP %in% a_thujone_snps_some_effect$SNP | SNP %in% b_thujone_snps_some_effect$SNP | SNP %in% sabinene_snps_some_effect$SNP | 
           SNP %in% a_thujaplicin_snps_some_effect$SNP | SNP %in% b_thujaplicin_snps_some_effect$SNP | SNP %in% g_thujaplicin_snps_some_effect$SNP | 
           SNP %in% b_thujaplicinol_snps_some_effect$SNP | SNP %in% ht15_snps_some_effect$SNP | SNP %in% dbh15_snps_some_effect$SNP))

table(droplevels(outlier_trait_annotation$Gene))

# Foliar terpene traits

outlier_trait_foliar_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & (SNP %in% a_thujone_snps_some_effect$SNP | SNP %in% b_thujone_snps_some_effect$SNP | SNP %in% sabinene_snps_some_effect$SNP))
                                     
table(droplevels(outlier_trait_foliar_annotation$Gene))
table(droplevels(outlier_trait_foliar_annotation$SNP))
unique(droplevels(outlier_trait_foliar_annotation$Gene))

outlier_athujone_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% a_thujone_snps_some_effect$SNP)

outlier_bthujone_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% b_thujone_snps_some_effect$SNP)

outlier_sabinene_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% sabinene_snps_some_effect$SNP)

table(droplevels(outlier_athujone_annotation$Gene))
table(droplevels(outlier_bthujone_annotation$Gene))
table(droplevels(outlier_sabinene_annotation$Gene))

table(droplevels(outlier_athujone_annotation$SNP))
table(droplevels(outlier_bthujone_annotation$SNP))
table(droplevels(outlier_sabinene_annotation$SNP))

# Wood terpene traits
outlier_trait_wood_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & (SNP %in% a_thujaplicin_snps_some_effect$SNP | SNP %in% b_thujaplicin_snps_some_effect$SNP | SNP %in% g_thujaplicin_snps_some_effect$SNP |
                                       SNP %in% b_thujaplicinol_snps_some_effect$SNP))

table(droplevels(outlier_trait_wood_annotation$Gene))
table(droplevels(outlier_trait_wood_annotation$SNP))
unique(droplevels(outlier_trait_wood_annotation$Gene))

outlier_athujaplicin_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% a_thujaplicin_snps_some_effect$SNP)

outlier_bthujaplicin_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% b_thujaplicin_snps_some_effect$SNP)

outlier_gthujaplicin_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% g_thujaplicin_snps_some_effect$SNP)

outlier_bthujaplicinol_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% b_thujaplicinol_snps_some_effect$SNP)

# Growth traits
outlier_trait_growth_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & (SNP %in% ht15_snps_some_effect$SNP | SNP %in% dbh15_snps_some_effect$SNP))

table(droplevels(outlier_trait_growth_annotation$Gene))
table(droplevels(outlier_trait_growth_annotation$SNP))
unique(droplevels(outlier_trait_growth_annotation$Gene))

outlier_ht15_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% ht15_snps_some_effect$SNP)

outlier_dbh15_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% dbh15_snps_some_effect$SNP)


# Outlier annotation with top effect
outlier_athujone_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% a_thujone_snps_top_effect$SNP)

outlier_bthujone_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% b_thujone_snps_top_effect$SNP)

outlier_sabinene_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% sabinene_snps_top_effect$SNP)

outlier_athujaplicin_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% a_thujaplicin_snps_top_effect$SNP)

outlier_bthujaplicin_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% b_thujaplicin_snps_top_effect$SNP)

outlier_bthujaplicinol_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% b_thujaplicinol_snps_top_effect$SNP)

outlier_gthujaplicin_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% g_thujaplicin_snps_top_effect$SNP)

outlier_ht15_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% ht15_snps_top_effect$SNP)

outlier_dbh15_top_annotation <- all_snp_annotations %>% 
  filter(SNP %in% outlier_snps$V1 & SNP %in% dbh15_snps_top_effect$SNP)

outlier_snps %>% filter(V1 %in% outlier_athujaplicin_snps$SNP | V1 %in% outlier_bthujaplicin_snps$SNP | V1 %in% outlier_gthujaplicin_snps$SNP | V1 %in% outlier_bthujaplicinol_snps$SNP)

outlier_snps %>% filter(V1 %in% outlier_athujone_snps$SNP | V1 %in% outlier_bthujone_snps$SNP| V1 %in% outlier_sabinene_snps$SNP)

dbh_ht_outlier_annotations <- all_snp_annotations %>% filter(SNP %in% outlier_snps$V1) %>% 
  filter(SNP %in% outlier_ht15_snps$SNP | SNP %in% outlier_dbh15_snps$SNP)

