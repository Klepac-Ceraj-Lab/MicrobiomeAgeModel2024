#####
# Code for Figure 1 - Community-level descriptive analysis
#####

library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
library(viridisLite)

loaded_data <- read.csv("/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/results/2024AgeModelManuscript/combined_inputs.csv")
# loaded_data <- read.csv("/home/guilherme/Repos/Leap/results/2024_06_26_abundance_allvars_prescree/combined_inputs.csv")
clock_md <- loaded_data[, c("ageMonths", "site", "datasource")]
clock_profiles <- loaded_data[, 11:(ncol(loaded_data)-1)]
clock_profiles <- clock_profiles / 100.0

D_clock <- vegdist(clock_profiles)
set.seed(1)

fit_adonis_clock_onlyds <- adonis2(D_clock ~ datasource, data = clock_md)
print(fit_adonis_clock_onlyds)

fit_adonis_clock_onlyage <- adonis2(D_clock ~ ageMonths, data = clock_md)
print(fit_adonis_clock_onlyage)

fit_adonis_clock_dsage <- adonis2(D_clock ~ datasource + ageMonths, data = clock_md)
print(fit_adonis_clock_dsage)

fit_adonis_clock_siteage <- adonis2(D_clock ~ site + ageMonths, data = clock_md)
print(fit_adonis_clock_siteage)

fit_adonis_clock_ageds <- adonis2(D_clock ~ ageMonths + datasource, data = clock_md)
print(fit_adonis_clock_ageds)
