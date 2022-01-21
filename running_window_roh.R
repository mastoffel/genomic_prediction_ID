# initial analysis: running window ROH on chromosome 14
library(tidyverse)
library(data.table)
library(glue)
library(lme4)
library(furrr)
# ROH > 1Mb
roh <- fread("output/ROH/roh.hom") %>% 
        filter(CHR == 14) %>% 
        rename(ID = IID) %>% 
        mutate(ID = as.character(ID)) 

# Chr lengths
chr_data <- read_delim("data/raw/chromosome_info_oar31.txt", delim = "\t") %>% 
        rename(size_BP = Length,
               CHR = Part) %>% 
        mutate(size_KB = size_BP / 1000)

autosomal_genome_size <- chr_data %>% 
        .[2:27, ] %>% 
        summarise(sum_KB = sum(size_KB)) %>% 
        as.numeric()

chr14_length <- chr_data %>% filter(CHR == "Chromosome 14") %>% .[["size_BP"]]

# get ROH genotypes for running windows along chromosome 14
# make windows
window_maker <- function(chr_size, window_size, step_size=window_size/2){
        start_bp <- seq(0, chr_size-window_size, step_size)
        end_bp <- seq(window_size, chr_size, step_size)
        return(tibble(start_bp, end_bp))
}
windows <- window_maker(chr14_length, 5e5, 0.1 * 1e6)

# get roh genotypes for each window
roh_geno <- function(start_bp, end_bp){
        start_mb <- start_bp/1e6
        end_mb <- end_bp/1e6
        roh %>% 
                mutate(roh_overlap = as.numeric(roh$POS1 <= start_bp & roh$POS2 >= end_bp)) %>% 
                group_by(ID) %>% 
                summarise("roh_{start_mb}_{end_mb}" := sum(roh_overlap))
}
# genos to table
roh_genos <- pmap(windows, roh_geno) %>% 
                reduce(left_join)

# PCs
pcs <- read_delim("../sheep_ID/data/ann_surv_pca.txt", " ", col_names = TRUE) %>% 
        mutate(id = as.character(id))

# fitness data
load("../sheep_ID/data/survival_mods_data.RData")
fitness <- fitness_data %>% 
        filter_at(vars(survival, froh_all, birth_year, sheep_year,
                       sex, twin, age), ~ !is.na(.)) %>% 
        mutate(age_std = as.numeric(scale(age)),
               age_std2 = age_std^2,
               lamb = ifelse(age == 0, 1, 0)) %>% 
        dplyr::select(id, survival, sex, twin, lamb, birth_year, sheep_year, age_std, age_std2, froh_no_chr14) %>%
        as_tibble() %>% 
        rename(ID = id) %>% 
        mutate(survival = ifelse(survival == 0, 1, 2)) %>% 
        mutate(ID = as.character(ID)) %>% 
        left_join(pcs, by = c("ID" = "id"))

# combine
gwas_df <- fitness %>% 
                #left_join(pcs) %>% 
                left_join(roh_genos) %>% 
                imputeTS::na_mean()
   
as_jwas <- gwas_df %>% 
                select(ID:pc7)

as_genos <- gwas_df %>% 
                select(one_of(names(roh_genos))) %>% 
                group_by(ID) %>% 
                slice(1) %>% 
                ungroup()

#names(as_genos)[2:ncol(as_genos)] <- str_replace_all(names(as_genos)[2:ncol(as_genos)], '\\.', '')
#as_genos <- as_genos[map_lgl(as_genos, function(x) any(x != 0))]
# save for jwas
write_delim(as_jwas, "data/processed/jwas_phenos_roh_windows.txt", delim = ",")
write_delim(as_genos, "data/processed/jwas_genos_roh_windows.txt", delim = "," )

# models
nlopt <- function(par, fn, lower, upper, control) {
        .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper, 
                                  opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
                                              maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
        list(par = res$solution,
             fval = res$objective,
             conv = if (res$status > 0) 0 else res$status,
             message = res$message
        )
}

run_gwas <- function(window, df) {
        # for mean froh without focal chr
        #chr <- as.numeric(snps_map_sub[snps_map_sub$snp.name == snp, "chromosome"])
        #froh_no_chr <- paste0("froh_no_chr", chr)
        print(window)
        df <- cbind(df[, 1:18], df[names(df) == window])
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age_std2 + froh_no_chr14", 
                                         " + ",
                                         #"pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         #"pc1 + pc2 + pc3 + pc4 +",
                                         window, "+ (1|sheep_year) + (1|id)")) # (1|birth_year) + 
        #snp, "+ ", paste0("roh_", snp), " + (1|sheep_year) + (1|id)"))
        mod <- glmer(formula = formula_snp,
                     data = df, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}

# running window names
windows <- names(roh_genos)[2:ncol(roh_genos)]

#plan(multisession, workers = 4)
gwas <- map(windows, run_gwas, gwas_df)

saveRDS(gwas, file = "output/gwas.rds")

all_windows <- gwas %>% 
        bind_rows(.id = "window") %>% 
        filter(str_detect(term, "^roh")) 

plot(all_windows$estimate)
all_windows %>% 
        arrange(p.value)
