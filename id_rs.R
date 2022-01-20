# test for inbreeding depression in breeding success
library(tidyverse)
library(here)
library(lme4)
library(broom.mixed)
library(performance)
library(sjPlot)

# fitness and froh data
load(here("data", "processed", "fitness_roh.RData"))
fitness <- fitness_data %>% 
                mutate(age_std = scale(age),
                       age2_std = age_std^2) 

fitness_f <- fitness %>% 
        filter(sex == "F") %>% 
        mutate(rs = ifelse(offspring_born > 0, 1, 0)) 

fitness_m <- fitness %>% 
        filter(sex == "M") %>% 
        mutate(rs = offspring_born,
               olre = as.factor(1:nrow(.)))
        
fit_f <- glmer(rs ~ froh_all + twin + age_std + age2_std + (1|sheep_year) + (1|birth_year) + (1|id),
             data = fitness_f, family = binomial)
fit_f %>% 
        tidy(conf.int = TRUE) %>% 
        mutate_if(is.numeric, round, 5)
plot_model(fit_f, type = "pred", terms = "froh_all [all]", show.data = TRUE)
