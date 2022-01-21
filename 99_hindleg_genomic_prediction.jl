using JWAS, CSV, DataFrames, StatsBase, Statistics
using Queryverse
# load phenotypes
phenotypes = CSV.read("data/processed/jwas_phenos_leg.txt",DataFrame,delim = ',',header=true)
# add permanent environment effect
phenotypes[!,:pe]=phenotypes[!,:ID]
show(phenotypes, allcols=true)

# load genotypes
genotypes  = get_genotypes("data/processed/jwas_genos_leg.txt", method = "BayesC", separator = ',', header=true, quality_control=false,
                            estimateScale = true)
                            # pedigree
pedigree   = get_pedigree("data/processed/pedigree.txt",separator=",",header=false)
pedigree

model_equation = "hindleg = intercept + sex + age_std + ID + pe + sheep_year + genotypes"
model=build_model(model_equation)

# set covariates (quantitative vars)
set_covariate(model, "age_std")
set_random(model, "ID", pedigree)
set_random(model, "sheep_year")
set_random(model, "pe")

# run model
out=runMCMC(model, phenotypes, burnin = 10000, pedigree = true, missing_phenotypes = false, chain_length=100000) # chain_length=100000
out
results    = innerjoin(out["EBV_hindleg"], phenotypes, on = :ID) 
results
accuracy  = cor(results[!,:EBV],results[!,:hindleg])

keys(out)

# gwas
marker_effects_file="results/MCMC_samples_marker_effects_genotypes_hindleg.txt"
out = GWAS(marker_effects_file,header=true)

map_file="data/processed/map.csv"
marker_effects_file="results3/MCMC_samples_marker_effects_genotypes_hindleg.txt"
out=GWAS(model,map_file,marker_effects_file,header=true,window_size="1 Mb")
out

#CSV.write("output/gwas.csv", out)
CSV.write("gwas.csv", out[1:5, 1:5])
save("output/gwas.csv", out)
writetable("output/gwas.csv", out)
out
CSV.write("out.csv", out)
DataFrame(rand(5, 2))