using JWAS, CSV, DataFrames, StatsBase, Statistics
using Queryverse
# load phenotypes
phenotypes = CSV.read("data/processed/jwas_phenos_surv.txt",DataFrame,delim = ',',header=true)
# add permanent environment effect
phenotypes[!,:pe]=phenotypes[!,:ID]
show(phenotypes, allcols=true)

# load genotypes
genotypes1  = get_genotypes("data/processed/jwas_genos_surv.txt", method = "BayesA", separator = ',', header=true, quality_control=false,
                            estimateScale = true, estimatePi = false)
genotypes2  = get_genotypes("data/processed/jwas_roh_0_surv.txt", method = "BayesC", separator = ',', header=true, quality_control=false,
                            estimateScale = true)       
genotypes3  = get_genotypes("data/processed/jwas_roh_2_surv.txt", method = "BayesC", separator = ',', header=true, quality_control=false,
                            estimateScale = true)               # pedigree
pedigree   = get_pedigree("data/processed/pedigree.txt",separator=",",header=false)
pedigree

model_equation = "survival = intercept + sex + froh_no_chr14 + age_std + age_std2 + lamb + twin + ID + pe + sheep_year + birth_year + genotypes1 + genotypes2 + genotypes3"
model=build_model(model_equation)

# set covariates (quantitative vars)
set_covariate(model, "age_std", "age_std2", "froh_no_chr14")
set_random(model, "ID", pedigree)
set_random(model, "sheep_year")
set_random(model, "pe")
set_random(model, "birth_year")
# run model

out=runMCMC(model, phenotypes, burnin = 10000, pedigree = true, missing_phenotypes = false, chain_length=200000,
            categorical_trait=true) # chain_length=100000


results    = innerjoin(out["EBV_survival"], phenotypes, on = :ID) 
results
accuracy  = cor(results[!,:EBV],results[!,:survival])

keys(out)

# gwas
marker_effects_file="results/MCMC_samples_marker_effects_genotypes_hindleg.txt"
out = GWAS(marker_effects_file,header=true)

