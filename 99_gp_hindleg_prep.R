# get subset of data for survival gwas to test JWAS 
library(lme4)
library(tidyverse)
library(broom.mixed)
library(snpStats)
library(data.table)
library(furrr)
library(here)
library(glue)
source("gg_themes.R")

# fitness data
load("../sheep_ID/data/survival_mods_data.RData")
fitness <- fitness_data %>% 
        filter_at(vars(hindleg, froh_all, birth_year, sheep_year), ~ !is.na(.)) %>% 
        # only adult sheep
        filter(age > 1) %>% 
        mutate(age_cent = age - mean(age, na.rm = TRUE),
               age_cent2 = age_cent^2,
               age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               # times 10 to estimate a 10% percent increase
               froh_all10 = froh_all * 10,
               froh_all10_cent = froh_all10 - mean(froh_all10, na.rm = TRUE),
               lamb = ifelse(age == 0, 1, 0),
               lamb_cent = lamb - mean(lamb, na.rm = TRUE),
               lamb = as.factor(lamb)) %>% 
        as.data.frame() 

load("../sheep_ID/data/sheep_ped.RData")
pedigree <- sheep_ped %>% 
                mutate(dam = replace_na(MOTHER, 0),
                       sire = replace_na(FATHER, 0)) %>% 
                select(ID, sire, dam) %>% 
        write_delim("data/processed/pedigree.txt", delim = ",", col_names = FALSE)

# roh
# roh data
#file_path <- "../sheep_ID/output/ROH/roh.hom"
#roh_lengths <- fread(file_path)

# snps
# plink name
sheep_plink_name <- "../sheep_ID/data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map

# top hit leg length "s23172.1"
snps_map_sub <- snps_map %>% 
        filter(chromosome == 16) %>% 
        filter((position < 69726554+2500000) & (position > 69726554-2500000))

# join additive and roh data to survival for gwas
leg_gwas <- fitness %>% 
        dplyr::select(id,hindleg, sex, sheep_year, age_std) %>% 
        as_tibble() %>% 
        mutate(sex = ifelse(sex == "M", 0, 1)) %>% 
        rename(ID = id) %>% 
        mutate(ID = as.character(ID))

library(imputeTS)
# additive genotypes
genos <- as_tibble(as(full_sample$genotypes[, snps_map_sub$snp], Class = "numeric"),
                   rownames = "id") %>% 
        na_mean() %>% 
        filter(id %in% leg_gwas$ID) %>% 
        #right_join(leg_gwas["ID"], by = c("id" = "ID")) %>% 
        rename(ID = id)

# genotypes
write_delim(genos, "data/processed/jwas_genos_leg.txt", delim = ",")
write_delim(leg_gwas, "data/processed/jwas_phenos_leg.txt", delim = ",")
