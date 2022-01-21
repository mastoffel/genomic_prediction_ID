# idea: call haplotypes in running windows and fit additive/dominance effect
# for each haplotype

library(data.table)
library(tidyverse)
#library(collapse)
library(glue)
library(here)
library(furrr)
library(data.table)
library(lme4)
# get individuals with genotypes
inds <- fread(here("data", "plink", "sheep.fam")) %>% 
        select(V2) %>% 
        rename(id = V2) %>% 
        .$id

chr <- 14

# get haplotypes, remove everything but genotypes to make it a matrix
haps_raw <- fread(here("data", "processed", glue("sheep_phased_50K_chr_{chr}.hap.gz"))) %>% 
        select(-V1, -V2, -V3, -V4, -V5) %>% 
        as.matrix()

ind_names <- system(glue("zgrep '^#CHROM*' ../haplotype_homozygosity/output/phased/sheep_phased_50K_chr_{chr}.vcf.gz"),
                    intern = TRUE) %>% 
        str_split("\t") %>%
        unlist() %>% 
        .[-(1:9)] %>% 
        str_split("_") %>% 
        map_chr(2)

# double each name and add _a _b for haplotype 1 / haplotype 2
ind_names_hap <- rep(ind_names, each = 2) %>% 
        paste0(., c("_a", "_b"))

# add to matrix
colnames(haps_raw) <- ind_names_hap

start_snp <- 1
hap_length <- 30

# reshape 
call_haps <- function(start_snp, hap_length, haps_raw) {
        
        if ((start_snp + hap_length-1) > nrow(haps_raw)) {
                stop_snp <- nrow(haps_raw)
        } else {
                stop_snp <- start_snp+hap_length-1
        }
        haps <- apply(haps_raw[start_snp:stop_snp, ,drop=F], 
                      paste, collapse = "", MARGIN = 2) %>% 
                enframe(name = "id_hap", value = "hap") %>% 
                mutate(id = str_remove(id_hap, "_[a-z]")) %>%
                select(-id_hap) %>% 
                mutate(hap_pos = rep(c("hap_a", "hap_b"), nrow(.)/2), .before = 1) %>% 
                mutate(id = as.numeric(id)) %>% 
                pivot_wider(names_from = hap_pos, values_from = hap)
        return(haps)
}

haps_across_chr <- function(start_snp, hap_length, haps_raw) {
        
        haps <- call_haps(start_snp, hap_length, haps_raw)
        # list haplotypes
        haps_tab <- table(as.matrix(haps[, c("hap_a", "hap_b")]))
        # haplotypes with high frequencies hf > 1%
        haps_hf <- haps_tab[haps_tab/sum(haps_tab) > 0.05]
        # make data.frame with all haps
        haps_df <- tibble(id = paste0("hap_",start_snp,"_",start_snp+hap_length, "_", 
                                      1:length(haps_hf)),
                          hap = names(haps_hf),
                          hap_freq = haps_hf)
        
        # get additive / dominance genotypes for top haps
        hap_to_geno <- function(hap_id, haps_df, haps) {
                focal_hap <- haps_df %>% filter(id == hap_id) %>% .[["hap"]]
                haps %>% 
                        mutate(hap_a = hap_a == focal_hap,
                               hap_b = hap_b == focal_hap) %>% 
                        mutate("{hap_id}" := hap_a + hap_b) %>% 
                        select(-hap_a, -hap_b)
        }
        
        map(haps_df$id, hap_to_geno, haps_df, haps) %>% 
                reduce(left_join) %>% 
                mutate(id = as.character(id))
        
}

# non overlapping windows
hap_length <- 20
start_snps <- seq(1, nrow(haps_raw), hap_length)
# additive haps (0,1,2)
hap_genos <- map(start_snps, haps_across_chr, hap_length, haps_raw) %>% 
                reduce(left_join) 

# make haps with 0/1 for not homozygous/homozygous
hap_genos_hom <- hap_genos
hap_genos_hom[hap_genos_hom == 1] <- 0
hap_genos_hom[hap_genos_hom == 2] <- 1

# change names of additive haplotypes
names(hap_genos)[2:ncol(hap_genos)] <- paste0("add_",names(hap_genos)[2:ncol(hap_genos)])

write_delim(hap_genos %>% rename(ID = id), "data/processed/jwas_haps_add_windows.txt", delim = "," )
write_delim(hap_genos_hom %>% rename(ID = id), "data/processed/jwas_haps_windows.txt", delim = "," )

# modeling
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
        #rename(ID = id) %>% 
       # mutate(survival = ifelse(survival == 0, 1, 2)) %>% 
        mutate(id = as.character(id)) %>% 
        left_join(pcs, by = c("id"))


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

all_haps <- names(hap_genos_hom)[2:ncol(hap_genos_hom)]

run_mods <- function(hap){
        df <- fitness %>% 
                left_join(hap_genos_hom[c("id", hap)]) %>% 
                left_join(hap_genos[c("id", paste0("add_", hap))])
        
        formula_snp <- as.formula(paste0("survival ~ 1 + sex + twin + age_std + age_std2 + froh_no_chr14", 
                                         " + ", 
                                         #"pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + ",
                                         hap, "+", paste0("add_", hap), "+ (1|sheep_year) + (1|id)"))
        
        mod <- glmer(formula = formula_snp,
                     data = df, family = "binomial",
                     control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
        out <- broom.mixed::tidy(mod)
        out
}
