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
using MLJ
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
include("/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/notebooks/allcohorts_data_loading_nofeed.jl")
combined_inputs.richness = map(x -> sum(x .> 0.0), eachrow(Matrix(combined_inputs[:, 11:ncol(combined_inputs)-1])))
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
filtered_inputs.richness = map(x -> sum(x .> 0.0), eachrow(Matrix(filtered_inputs[:, 11:ncol(filtered_inputs)-2])))

subset!(filtered_inputs, :richness => x -> x .>= 1) # Minimum sample richness should be at least 1.
select!(filtered_inputs, Not(:richness))

CSV.write("manuscript/final_manuscript_inputs.csv", filtered_inputs)
```

##  RF Model training

### Actual function for model training
The following block of code will train the model on the combination of cohorts, performing crossvalidation and grid-search hyperparameter optimization. Training can take several hours if the hyperparameter grid is large. It is advised to train once and store the result on a `JLD2` object so it can be accessed with `JLD2.load`-like methods for downstream analysis and plotting. Hence, the block of code should be run only once per data update.
```julia
regression_Age_LeaveEchoOut = probe_regression_randomforest(
    "regression_Age_LeaveEchoOut",
    subset(filtered_inputs, :datasource => x -> x .!== "ECHO-RESONANCE"),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .== "ECHO-RESONANCE"),
    ext_firstinputcol = 11,
    ext_uniquecol = :sample,
    unique_col = :sample,
    n_folds = 5,
    n_replicas = 100,
    n_rngs = 5,
    tuning_space = (; #PRODUCTION
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        min_samples_split = [ 2 ],
        sampsize_range = [ 0.8 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 200 ]
    )
)

@show sort(report_regression_merits(regression_Age_LeaveEchoOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveEchoOut.jld") regression_Age_LeaveEchoOut

regression_Age_LeaveCMDOut = probe_regression_randomforest(
    "regression_Age_LeaveCMDOut",
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["CMD-OTHER", "CMD-DIABIMMUNE"])),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .∈ Ref(["CMD-OTHER", "CMD-DIABIMMUNE"])),
    ext_firstinputcol = 11,
    ext_uniquecol = :sample,
    unique_col = :sample,
    n_folds = 5,
    n_replicas = 100,
    n_rngs = 5,
    tuning_space = (; #PRODUCTION
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        min_samples_split = [ 2 ],
        sampsize_range = [ 0.8 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 200 ]
    )
)

@show sort(report_regression_merits(regression_Age_LeaveCMDOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveCMDOut.jld") regression_Age_LeaveCMDOut

regression_Age_LeaveKhulaOut = probe_regression_randomforest(
    "regression_Age_LeaveKhulaOut",
    subset(filtered_inputs, :datasource => x -> x .!== "1kDLEAP-KHULA"),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .== "1kDLEAP-KHULA"),
    ext_firstinputcol = 11,
    ext_uniquecol = :sample,
    unique_col = :sample,
    n_folds = 5,
    n_replicas = 100,
    n_rngs = 5,
    tuning_space = (; #PRODUCTION
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        min_samples_split = [ 2 ],
        sampsize_range = [ 0.8 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 200 ]
    )
)

@show sort(report_regression_merits(regression_Age_LeaveKhulaOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveKhulaOut.jld") regression_Age_LeaveKhulaOut

regression_Age_LeaveGerminaOut = probe_regression_randomforest(
    "regression_Age_LeaveGerminaOut",
    subset(filtered_inputs, :datasource => x -> x .!== "1kDLEAP-GERMINA"),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .== "1kDLEAP-GERMINA"),
    ext_firstinputcol = 11,
    ext_uniquecol = :sample,
    unique_col = :sample,
    n_folds = 5,
    n_replicas = 100,
    n_rngs = 5,
    tuning_space = (; #PRODUCTION
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        min_samples_split = [ 2 ],
        sampsize_range = [ 0.8 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 200 ]
    )
)

@show sort(report_regression_merits(regression_Age_LeaveGerminaOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveGerminaOut.jld") regression_Age_LeaveGerminaOut

regression_Age_LeaveCombineOut = probe_regression_randomforest(
    "regression_Age_LeaveCombineOut",
    subset(filtered_inputs, :datasource => x -> x .!== "1kDLEAP-COMBINE"),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .== "1kDLEAP-COMBINE"),
    ext_firstinputcol = 11,
    ext_uniquecol = :sample,
    unique_col = :sample,
    n_folds = 5,
    n_replicas = 100,
    n_rngs = 5,
    tuning_space = (; #PRODUCTION
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        min_samples_split = [ 2 ],
        sampsize_range = [ 0.8 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 200 ]
    )
)

@show sort(report_regression_merits(regression_Age_LeaveCombineOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveCombineOut.jld") regression_Age_LeaveCombineOut

regression_Age_LeaveM4EFADOut = probe_regression_randomforest(
    "regression_Age_LeaveM4EFADOut",
    subset(filtered_inputs, :datasource => x -> x .!== "1kDLEAP-M4EFAD"),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .== "1kDLEAP-M4EFAD"),
    ext_firstinputcol = 11,
    ext_uniquecol = :sample,
    unique_col = :sample,
    n_folds = 5,
    n_replicas = 100,
    n_rngs = 5,
    tuning_space = (; #PRODUCTION
        maxnodes_range = [ -1 ],
        nodesize_range = [ 5 ],
        min_samples_split = [ 2 ],
        sampsize_range = [ 0.8 ],
        mtry_range = [ -1 ],
        ntrees_range = [ 200 ]
    )
)

@show sort(report_regression_merits(regression_Age_LeaveM4EFADOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveM4EFADOut.jld") regression_Age_LeaveM4EFADOut
```

After the code is run at least once, the results can then be loaded with:
```julia
JLD2.@load joinpath(outdir, "regression_Age_LeaveEchoOut.jld") regression_Age_LeaveEchoOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveCMDOut.jld") regression_Age_LeaveCMDOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveKhulaOut.jld") regression_Age_LeaveKhulaOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveGerminaOut.jld") regression_Age_LeaveGerminaOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveCombineOut.jld") regression_Age_LeaveCombineOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveM4EFADOut.jld") regression_Age_LeaveM4EFADOut
```

Reporting the joint figures of merit (Supplementary Table 1):
```julia
joint_merits_df = vcat(
    report_regression_merits(regression_Age_LeaveEchoOut),
    report_regression_merits(regression_Age_LeaveGerminaOut),
    report_regression_merits(regression_Age_LeaveKhulaOut),
    report_regression_merits(regression_Age_LeaveCombineOut),
    report_regression_merits(regression_Age_LeaveM4EFADOut),
    report_regression_merits(regression_Age_LeaveCMDOut)
)
println("RMSE of LOOCV: $(round(mean(joint_merits_df.Test_RMSE_mean); digits = 2)) +- $(round(conf_interval(joint_merits_df.Test_RMSE_mean); digits = 2))")
```

# Creating Supplementary Figure 2
```julia
supp_figure3_master = Figure(; size = (600, 600))

age_bins = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
results_vector = [ 
    regression_Age_LeaveEchoOut,
    regression_Age_LeaveGerminaOut,
    regression_Age_LeaveCombineOut,
    regression_Age_LeaveKhulaOut,
    regression_Age_LeaveM4EFADOut,
    regression_Age_LeaveCMDOut
]

cohort_pertinences = [ ["ECHO-RESONANCE"], ["1kDLEAP-GERMINA"], ["1kDLEAP-COMBINE"],["1kDLEAP-KHULA"], ["1kDLEAP-M4EFAD"], ["CMD-OTHER", "CMD-DIABIMMUNE"] ]

lococv_mat = Matrix{Float64}(undef, 6, 6)

for i in eachindex(results_vector) 
    for j in eachindex(results_vector)
        if i == j
            this_rmse = report_regression_merits(results_vector[i]).Test_RMSE_mean[1]
        else
            @show this_samples = @chain predictions_to_plot(results_vector[i], results_vector[i].original_df, age_bins) begin
                subset(:datasource => x -> x .∈ Ref(cohort_pertinences[j]))                        
            end
            this_rmse = MLJ.rmse(this_samples.test_prediction, this_samples.ageMonths)
        end
        lococv_mat[i,j] = this_rmse
    end
end

## Actual heatmap
ax = Axis(
    supp_figure3_master[1,1],
    xlabel = "Source Left Out",
    ylabel = "Source Metric",
    title = "Leave-One-Source-Out",
    xticks = (1:6, ["ECHO", "Germina", "Combine", "Khula", "M4EFaD", "CMD"]),
    yticks = (1:6, ["ECHO", "Germina", "Combine", "Khula", "M4EFaD", "CMD"]),
    xticklabelrotation = pi/4,
    yticklabelrotation = pi/4,
    yticklabelsize=20,
    xticklabelsize=20,
    titlesize = 24,
    aspect = 1
)
hm = heatmap!(ax, lococv_mat, colormap = cgrad(:magma, rev = false))

Colorbar(supp_figure3_master[1,2], hm, label = "RMSE")

save(joinpath(outdir, "figures", "FigureS3.png"), supp_figure3_master)
save(joinpath(outdir, "figures", "FigureS3.eps"), supp_figure3_master)
save(joinpath(outdir, "figures", "FigureS3.svg"), supp_figure3_master)
supp_figure3_master
```