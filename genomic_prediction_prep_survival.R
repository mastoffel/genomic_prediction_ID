# get subset of data for survival gwas to test JWAS 
library(lme4)
library(tidyverse)
library(broom.mixed)
library(snpStats)
library(data.table)
library(furrr)
library(here)
library(glue)

# fitness data
load("../sheep_ID/data/survival_mods_data.RData")
fitness <- fitness_data %>% 
        filter_at(vars(survival, froh_all, birth_year, sheep_year,
                        sex, twin, age), ~ !is.na(.)) %>% 
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
# pcs
# GRM PCs from plink
pcs <- read_delim("../sheep_ID/data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

# roh
# roh data
file_path <- "../sheep_ID/output/ROH/roh.hom"
roh_lengths <- fread(file_path)

# snps
# plink name
sheep_plink_name <- "../sheep_ID/data/sheep_geno_imputed_oar_filt"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_map <- full_sample$map

# gwas top snps
gwas_res <- read_rds("../sheep_ID/output/gwas_res_oar_roh_sep.rds")
#add map positions
gwas_out <- gwas_res %>% 
        rename(snp.name = snp) %>% 
        left_join(snps_map) %>% 
        rename(snp = snp.name)

# get region around top snp
snps_map_sub <- gwas_out %>% 
        filter(chromosome == 14) %>% 
        filter(state == "roh_2") %>% 
        filter((position < 45000000) & (position > 44000000))

# additive genotypes
genos <- as_tibble(as(full_sample$genotypes[, snps_map_sub$snp], Class = "numeric"),
                   rownames = "id") %>% 
        na_mean() %>% 
        rename(ID = id) %>% 
        na_mean() %>% 
        filter(ID %in% fitness$id)
   

# check whether snp is in roh for given individual
setkey(roh_lengths, IID)
roh_id_per_snp <- function(i) {
        position <- as.numeric(snps_map_sub[i, "position"])
        chromosome <- as.numeric(snps_map_sub[i, "chromosome"])
        roh_lengths[, roh := as.numeric((CHR == chromosome) & (POS1 <= position) & (POS2 >= position))]
        roh_id <- roh_lengths[,  .(roh = max(roh)), by = c("IID")]$roh
}
roh_ind <- map(1:nrow(snps_map_sub), roh_id_per_snp)
roh_df <- as.data.frame(do.call(cbind, roh_ind))
names(roh_df) <- paste0("roh_", snps_map_sub$snp)
roh_df$id <- as.character(unique(roh_lengths$IID))
roh_df <- roh_df %>% filter(id %in% fitness$id)

# which chromosome
chrs <- unique(snps_map_sub$chromosome)
froh_no_chr <- paste0("froh_no_chr", chrs)

# join additive and roh data to survival for gwas
annual_survival_gwas <- fitness %>% 
        dplyr::select(id, survival, sex, twin, lamb, birth_year, sheep_year, age_std, age_std2, {{ froh_no_chr }}) %>%
        as_tibble() %>% 
        rename(ID = id)

snp_names <- snps_map_sub$snp

roh_0_df <- map_dfc(snp_names, function(i) as.numeric((genos[[i]] == 0) & (roh_df[[paste0("roh_", i)]] == 1))) %>% 
                setNames(paste0("roh_0_", snp_names)) %>% 
                add_column(ID = genos$ID, .before = 1) 
roh_2_df <- map_dfc(snp_names, function(i) as.numeric((genos[[i]] == 2) & (roh_df[[paste0("roh_", i)]] == 1))) %>% 
                setNames(paste0("roh_2_", snp_names)) %>% 
                add_column(ID = genos$ID, .before = 1 )

# write phenotypes
write_delim(annual_survival_gwas %>% mutate(survival = ifelse(survival == 0, 1, 2)), 
            "data/processed/jwas_phenos_surv.txt", delim = ",")

# genotypes
write_delim(genos, "data/processed/jwas_genos_surv.txt", delim = ",")
write_delim(roh_0_df, "data/processed/jwas_roh_0_surv.txt", delim = ",")
write_delim(roh_2_df, "data/processed/jwas_roh_2_surv.txt", delim = ",")


