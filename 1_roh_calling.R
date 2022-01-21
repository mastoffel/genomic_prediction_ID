# call ROH > 1Mb
system(paste0("/usr/local/bin/plink --bfile data/processed/sheep_geno_imputed_oar_filt ",
              "--sheep --out output/ROH/roh ",
              "--homozyg --homozyg-window-snp 25 --homozyg-snp 25 --homozyg-kb 1000 ",
              "--homozyg-gap 500 --homozyg-density 100 --homozyg-window-missing 2 ",
              "--homozyg-het 2 ",
              "--homozyg-window-het 2"))
