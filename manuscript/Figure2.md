# Figure 2 - Model benchmark, Importance analysis

## Pre-configuration

### Loading Packages 
```julia
using Chain
using XLSX
using DataFrames
using MultivariateStats
using Microbiome
using Distances
using Diversity
using Random
using JLD2
using KernelDensity
using Statistics
using CairoMakie
using Leap
using CategoricalArrays
using GLM
using StatsBase
using StableRNGs
using Polynomials
```

### Configurable parameters
```julia
master_colors = Dict(
    "ECHO" => "purple",
    "ECHO-RESONANCE" => "purple",
    "1kDLEAP-GERMINA" => "blue",
    "1kDLEAP-CORK" => "orange",
    "1kDLEAP-COMBINE" => "orange",
    "1kDLEAP-KHULA" => "red",
    "1kDLEAP-M4EFAD" => "darkgreen",
    "DIABIMMUNE" => "lightblue",
    "CMD" => "lightblue"
)

experiment_name = "2024AgeModelManuscript"
outdir = joinpath(pwd(), "results", experiment_name)
figdir = joinpath(outdir, "figures")
deepdivemonodir, deepdivecolordir = ( joinpath(figdir, "species_monocolor_scatterplots"), joinpath(figdir, "species_colored_scatterplots") )
isdir(outdir) ? @warn("Directory $(outdir) already exists! This notebook will overwrite files already there.") : ( mkpath(outdir), mkpath(figdir), mkpath(deepdivemonodir), mkpath(deepdivecolordir) )
presence_absence = false # This argument will control whether the model will be based on abundances or binary presence/absence
```

## Loading data

### Loading taxonomic profiles from all the cohorts
This line will evoke the auxiliary notebook that contains the code to load data from all cohorts
```julia
include("notebooks/allcohorts_data_loading_nofeed.jl")
```

### Filtering and preparing the data
Computing samples and features before and after filters:

```julia
prevalence_threshold = 0.05

println("Table of combined inputs has $(nrow(combined_inputs)) samples and $(sum(map(sum, eachcol(combined_inputs[:, 11:end])) .> 0.0) - 2) detected taxa before prevalence filtering")
# -1 because one of the features is Shannon

println("Of those samples:\n\t$(sum(combined_inputs.richness .< 5 )) have 4 or less taxa detected;\n\t$(sum(combined_inputs.richness .< 4 )) have 3 or less taxa detected;\n\t$(sum(combined_inputs.richness .< 3 )) have 2 or less taxa detected.\n\t$(sum(combined_inputs.richness .< 2 )) have 1 or less taxa detected; \n\t$(sum(combined_inputs.richness .< 1 )) have 0 taxa detected; " )

println("Of the $(sum(map(sum, eachcol(combined_inputs[:, 11:end])) .> 0.0)) features, only $(ncol(filter_prevalence(combined_inputs, prevalence_threshold)[:, 11:end]) - 1) pass a universal prevalence of $(prevalence_threshold) filter")

println("After prevalence filtering, there are $(sum(map(sum, eachrow(filter_prevalence(combined_inputs, prevalence_threshold)[:, 11:end-2])) .== 0.0 ) ) samples that end up with no abundance on the remaining taxa")

filtered_inputs = filter_prevalence(combined_inputs, prevalence_threshold)

subset!(filtered_inputs, :richness => x -> x .>= 3) # Minimum smaple richness should be more than 1% of the final number of predictors (~150, posthoc), rounded to the ceiling (so, 2). Hence, richness has to be >= 3.
select!(filtered_inputs, Not(:richness))

CSV.write("manuscript/final_manuscript_inputs2.csv", filtered_inputs)
```

##  RF Model training

### Actual function for model training
The following block of code will train the model on the combination of cohorts, performing crossvalidation and grid-search hyperparameter optimization. Training can take several hours if the hyperparameter grid is large. It is advised to train once and store the result on a `JLD2` object so it can be accessed with `JLD2.load`-like methods for downstream analysis and plotting. Hence, the block of code should be run only once per data update.
```julia
# regression_Age_FullCV = probe_regression_randomforest(
#     "regression_Age_FullCV",
#     filtered_inputs,
#     identity,
#     collect(11:ncol(combined_inputs)),
#     :ageMonths;
#     split_strat = "subject",
#     custom_input_group = nothing,
#     unique_col = :sample,
#     n_folds = 5,
#     n_replicas = 50,
#     n_rngs = 5,
#     tuning_space = (; #PRODUCTION
#         maxnodes_range = [ -1 ],
#         nodesize_range = [ 5, 7 ],
#         min_samples_split = [ 2 ],
#         sampsize_range = [ 0.7, 0.8 ],
#         mtry_range = [ -1, 0, 10 ],
#         ntrees_range = [ 100, 200 ]
#     )    
# )
# @show sort(report_regression_merits(regression_Age_FullCV), :Val_RMSE_mean)
# JLD2.@save joinpath(outdir, "AgeModel_FullCV_Results.jld") regression_Age_FullCV
```

After the code is run at leat once, the results can then be loaded with:
```julia
JLD2.@load joinpath(outdir, "AgeModel_FullCV_Results.jld") regression_Age_FullCV
```
Assuming that `outdir` points to the same place where the model was stored when trained. If not, the argumetn can be customized accordingly. Remember that, per documentation, programatically-built filenames require explicit variable names on the macro to work. Loading variables into current scope requires literal file names.

# Creating Master Figure 2
```julia
figure2_master = Figure(; size = (1200, 850))

A_Subfig = GridLayout(figure2_master[1,1], alignmode=Inside()) 
B_Subfig = GridLayout(figure2_master[1,2], alignmode=Inside())
CDEFG_Subfig  = GridLayout(figure2_master[2,1:2], alignmode=Inside())
```

## Scatterplot
```julia
age_bins = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
hp_idx = 15 # checked with `sort(report_regression_merits(regression_Age_FullCV), :Val_RMSE_mean)`

axA = Axis(
    A_Subfig[1, 1];
    xlabel = "Sample collection age (months)",
    xticks = (2:1:18),
    ylabel = "Microbial age (months)",
    yticks = (2:1:18),
    aspect = AxisAspect(1.0)
    # alignmode = Outside(),
    # title = "Predictions for Test/Validation data"
)
hidedecorations!(axA, label = false, ticklabels = false, ticks = false, minorgrid = true, minorticks = true)
xlims!(axA, [1.99, 18.1]); ylims!(axA, [1.99, 18.1])
@chain regression_Age_FullCV begin
    predictions_to_plot(filtered_inputs, age_bins, "val"; hp = hp_idx)
    scatter!(axA, _[:, :ageMonths],  _[:, :test_prediction], color = [ (ccol, 0.6) for ccol in _[:, "datacolor"] ], marker = :circle)
end
ablines!(axA, 0, 1; linestyle = :dash, linewidth=2, color = :gray )
@chain regression_Age_FullCV begin
    report_regression_merits()
    sort(:Val_RMSE_mean)
    annotations!(
        axA,
        [
            "RMSE: " * string(round(_[:, :Val_RMSE_mean][1]; digits = 4)) * " (months)",
            "r: " * string(round(_[:, :Val_Cor_mean][1]; digits = 4))
        ],
        [Point(15.0, 3.5), Point(15.0, 2.5)];
        fontsize = 14,
        align = (:center, :bottom)
    )
end

Legend(
    A_Subfig[2, 1],
    [
        MarkerElement(marker = :circle, color = :lightblue, markersize = 14),
        MarkerElement(marker = :circle, color = :purple, markersize = 14),
        MarkerElement(marker = :circle, color = :blue, markersize = 14),
        # MarkerElement(marker = :circle, color = :darkgreen, markersize = 14),
        MarkerElement(marker = :circle, color = :red, markersize = 14),
        MarkerElement(marker = :circle, color = :orange, markersize = 14),
        MarkerElement(marker = :circle, color = :darkgreen, markersize = 14)
    ],
    [
        "CMD",
        "ECHO-Resonance",
        "1kDLEAP-Brainrise",
        # "1kDLEAP-KhulaMW",
        "1kDLEAP-KhulaSA",
        "1kDLEAP-Combine",
        "1kDLEAP-M4EFaD"
    ],
    orientation = :vertical,
    nbanks = 3,
    labelsize = 12,
    tellheight = true,
    tellwidth = true
)
```

## Figure 2, Panel B - Importance plots
```julia
importances_table = hpimportances(regression_Age_FullCV, hp_idx)[1:31,:]
importances_table.correl = [ cor(filtered_inputs[:, ccol], filtered_inputs.ageMonths) for ccol in importances_table.variable ]
importances_table.impsign = importances_table.weightedImportance .* sign.(importances_table.correl)

onlyspecies_importances = subset(importances_table, :variable => x -> x .!= "shannon_index")

println("Number of species positively correlated with age: $(sum( onlyspecies_importances.correl .> 0.0 ))")
println("Mean positive correlation: $(mean( onlyspecies_importances.correl[ onlyspecies_importances.correl .> 0.0 ])) +- $(Statistics.std( onlyspecies_importances.correl[ onlyspecies_importances.correl .> 0.0 ]))")
println("Mean negative correlation: $(mean( onlyspecies_importances.correl[ onlyspecies_importances.correl .< 0.0 ])) +- $(Statistics.std( onlyspecies_importances.correl[ onlyspecies_importances.correl .< 0.0 ]))")

## bugs present in every cohort?
cohorts = unique(filtered_inputs.datasource)
bugs = onlyspecies_importances.variable
check_prevalences_mat = Matrix{Float64}(undef, length(bugs), length(cohorts))

for (j, cohort) in enumerate(cohorts)
    for (i, bug) in enumerate(bugs)
        check_prevalences_mat[i,j] = sum( subset(filtered_inputs, :datasource => x -> x .== cohort)[:, bug] .> 0.0 )
    end
end

@show findall(check_prevalences_mat .== 0.0)

for bug in ["Faecalibacterium_prausnitzii", "Anaerostipes_hadrus", "Flavonifractor_plautii", "Eubacterium_rectale", "Bifidobacterium_longum", "Bifidobacterium_breve", "Ruminococcus_gnavus"]
    println("$(bug): $(cor(filtered_inputs[:, bug], filtered_inputs.ageMonths))")
end

## If doing signs, watch comment on next line!
sort!(importances_table, :impsign; rev = true)

axB = Axis(
    B_Subfig[1, 1];
    xlabel = "sign(R) x proportional importance",
    yticks = (reverse(collect(1:31)), importances_table.variable[1:31]),
    ylabel = "Predictor"
)

barplot!(
    axB,
    reverse(collect(1:31)),
    importances_table.impsign[1:31],
    color = [ ( (el > 0) ? "blue" : "red" ) for el in importances_table.impsign[1:31] ],
    direction=:x
)
```

## Add labels
```julia
Label(A_Subfig[1, 1, TopLeft()], "A", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(B_Subfig[1, 1, TopLeft()], "B", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(CDEFG_Subfig[1, 1, TopLeft()], "C", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(CDEFG_Subfig[1, 2, TopLeft()], "D", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(CDEFG_Subfig[1, 3, TopLeft()], "E", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(CDEFG_Subfig[1, 4, TopLeft()], "F", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(CDEFG_Subfig[1, 5, TopLeft()], "G", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
```

## Fix layout
```julia
colgap!(figure2_master.layout, 5)
rowgap!(figure2_master.layout, 5)
colgap!(A_Subfig, 5)
rowgap!(A_Subfig, 5)
colgap!(B_Subfig, 5)
rowgap!(B_Subfig, 5)
colgap!(CDEFG_Subfig, 5)
rowgap!(CDEFG_Subfig, 5)

# colsize!(figure1_master.layout, 2, Relative(0.4))
rowsize!(figure2_master.layout, 2, Relative(0.3))
# rowsize!(AB_Subfig, 1, Relative(0.50))
# rowsize!(AB_Subfig, 2, Relative(0.15))
# rowsize!(AB_Subfig, 3, Relative(0.35))
```

# Export Figure 2

```julia
save(joinpath(outdir, "figures", "Figure2.png"), figure2_master)
figure2_master
```