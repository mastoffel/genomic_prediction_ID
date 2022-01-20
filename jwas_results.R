# genomic prediction results
library(tidyverse)

marker_effs <- read_delim("results12/marker_effects_genotypes1.txt")
marker_effs %>% 
        arrange(desc(Model_Frequency))

ebvs <- read_delim("results2/EBV_hindleg.txt")

mcmc <- read_delim("results/MCMC_samples_heritability.txt", ",")

lines(mcmc$s23172.1)

h <- read_delim("results2/MCMC_samples_heritability.txt", ",")
marker_effs %>% 
        filter(Marker_ID == "s23172.1")

plot(marker_effs$Estimate)

marker_effs %>% 
        arrange(Estimate)
snps_map_sub %>% 
        left_join(marker_effs, by = c("snp.name" = "Marker_ID")) %>% 
        arrange(Model_Frequency)

marker_effs %>% 
        arrange(desc(Model_Frequency))

snps_map %>% 
        select(snp.name, chromosome, position) %>% 
        filter(chromosome == 16) %>% 
        rename(markerID = snp.name) %>% 
        filter(markerID %in% marker_effs$Marker_ID) %>% 
        as_tibble() %>% 
        write_delim("data/processed/map.csv", delim = ",")

snps_map %>% 
        filter(snp.name == "s23172.1")
snps_map_sub <- snps_map %>% 
                filter(chromosome == 16)
which(snps_map_sub$snp.name == "s23172.1")
