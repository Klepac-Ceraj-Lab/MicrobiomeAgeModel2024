# Figure 2 - Model benchmark, Importance analysis

## Pre-configuration

### Loading Packages 
```julia
using Chain
using XLSX
using DataFrames
using MultivariateStats
using Microbiome
using Clustering
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
    "1kDLEAP-BRAINRISE" => "blue",
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
```julia
JLD2.@load joinpath(outdir, "AgeModel_FullCV_Results.jld") regression_Age_FullCV
taxonomic_profiles = regression_Age_FullCV.original_df

bins = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
```

## Finding the important predictors
```julia
# @show sort(report_regression_merits(regression_Age_FullCV), :Val_RMSE_mean) # To check the nest hyperparameter index
hp_idx = 15

importances_table = @chain Leap.hpimportances(regression_Age_FullCV, hp_idx) begin
    subset(:variable => x -> x .!= "richness")
    subset(:variable => x -> x .!= "shannon_index")
end

important_bugs = importances_table[1:30, :variable]
```

## Idea 1. Abundance and Prevalence Heatmap - CMD, ECHO and Khula samples, top most important variables
```julia
interval_bounds = collect(zip(bins[1:end-1], bins[2:end]))

full_abundance_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))
full_prevalence_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))

cmd_abundance_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))
cmd_prevalence_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))

khula_abundance_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))
khula_prevalence_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))

echo_abundance_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))
echo_prevalence_matrix = zeros(Float64, length(important_bugs), length(interval_bounds))

for (bnds_idx, bnds) in enumerate(interval_bounds)
    for (bug_idx, bug_name) in enumerate(important_bugs)

        bugvec = subset(taxonomic_profiles, :ageMonths => x -> ( (bnds[1] .< x) .& ( x .<= bnds[2]) ))
        bugvec = select(bugvec, bug_name)[:,1]

        if isempty(bugvec[bugvec .!= 0.0])
            full_abundance_matrix[bug_idx, bnds_idx] = 0.0
            full_prevalence_matrix[bug_idx, bnds_idx] = 0.0
            continue
        end
        
        full_abundance_matrix[bug_idx, bnds_idx] = mean(bugvec[bugvec .!= 0.0])
        full_prevalence_matrix[bug_idx, bnds_idx] = mean(bugvec .!= 0.0)

    end
end

for (bnds_idx, bnds) in enumerate(interval_bounds)
    for (bug_idx, bug_name) in enumerate(important_bugs)

        ## EST + RUS + FIN
        # cmd_bugvec = subset(taxonomic_profiles, :datasource => x -> x .== "DIABIMMUNE")
        cmd_bugvec = subset(taxonomic_profiles, :site => x -> x .∈ Ref(["EST", "RUS", "FIN", "SWE"]))
        cmd_bugvec = subset(cmd_bugvec, :ageMonths => x -> ( (bnds[1] .< x) .& ( x .<= bnds[2]) ))
        cmd_bugvec = select(cmd_bugvec, bug_name)[:,1]

        if isempty(cmd_bugvec[cmd_bugvec .!= 0.0])
            cmd_abundance_matrix[bug_idx, bnds_idx] = 0.0
            cmd_prevalence_matrix[bug_idx, bnds_idx] = 0.0
        else
            cmd_abundance_matrix[bug_idx, bnds_idx] = mean(cmd_bugvec[cmd_bugvec .!= 0.0])
            cmd_prevalence_matrix[bug_idx, bnds_idx] = mean(cmd_bugvec .!= 0.0)
        end
        

        ## KHULA
        khula_bugvec = subset(taxonomic_profiles, :datasource => x -> x .== "1kDLEAP-KHULA")
        khula_bugvec = subset(khula_bugvec, :ageMonths => x -> ( (bnds[1] .< x) .& ( x .<= bnds[2]) ))
        khula_bugvec = select(khula_bugvec, bug_name)[:,1]

        if isempty(khula_bugvec[khula_bugvec .!= 0.0])
            khula_abundance_matrix[bug_idx, bnds_idx] = 0.0
            khula_prevalence_matrix[bug_idx, bnds_idx] = 0.0
        else
            khula_abundance_matrix[bug_idx, bnds_idx] = mean(khula_bugvec[khula_bugvec .!= 0.0])
            khula_prevalence_matrix[bug_idx, bnds_idx] = mean(khula_bugvec .!= 0.0)    
        end
        

        ## ECHO
        # echo_bugvec = subset(taxonomic_profiles, :datasource => x -> x .== "ECHO")
        echo_bugvec = subset(taxonomic_profiles, :site => x -> x .∈ Ref(["USA"]))
        echo_bugvec = subset(echo_bugvec, :ageMonths => x -> ( (bnds[1] .< x) .& ( x .<= bnds[2]) ))
        echo_bugvec = select(echo_bugvec, bug_name)[:,1]

        if isempty(echo_bugvec[echo_bugvec .!= 0.0])
            echo_abundance_matrix[bug_idx, bnds_idx] = 0.0
            echo_prevalence_matrix[bug_idx, bnds_idx] = 0.0
        else
            echo_abundance_matrix[bug_idx, bnds_idx] = mean(echo_bugvec[echo_bugvec .!= 0.0])
            echo_prevalence_matrix[bug_idx, bnds_idx] = mean(echo_bugvec .!= 0.0)
        end
        
    end
end

## Computing vectors and correlations of matrices
khula_prevalence_vector = vec(khula_prevalence_matrix)
cmd_prevalence_vector = vec(cmd_prevalence_matrix)
echo_prevalence_vector = vec(echo_prevalence_matrix)

@show cor(echo_prevalence_vector, cmd_prevalence_vector)
@show cor(echo_prevalence_vector, khula_prevalence_vector)
@show cor(cmd_prevalence_vector, khula_prevalence_vector)

## Remember the idea of the boxplots# table (bugs x 3) then calculate all the vector correlations and plot boxplots, see which ones are actually close or not.

## Computing geomeans

full_geomeans_matrix = zeros(Float64, size(full_abundance_matrix))
for i in 1:size(full_abundance_matrix,1) for j in 1:size(full_abundance_matrix,2)
    full_geomeans_matrix[i,j] = sqrt( (full_abundance_matrix[i,j]/100) * full_prevalence_matrix[i,j] )
end end

cmd_geomeans_matrix = zeros(Float64, size(cmd_abundance_matrix))
for i in 1:size(cmd_abundance_matrix,1) for j in 1:size(cmd_abundance_matrix,2)
    cmd_geomeans_matrix[i,j] = sqrt( (cmd_abundance_matrix[i,j]/100) * cmd_prevalence_matrix[i,j] )
end end

khula_geomeans_matrix = zeros(Float64, size(khula_abundance_matrix))
for i in 1:size(khula_abundance_matrix,1) for j in 1:size(khula_abundance_matrix,2)
    khula_geomeans_matrix[i,j] = sqrt( (khula_abundance_matrix[i,j]/100) * khula_prevalence_matrix[i,j] )
end end

echo_geomeans_matrix = zeros(Float64, size(echo_abundance_matrix))
for i in 1:size(echo_abundance_matrix,1) for j in 1:size(echo_abundance_matrix,2)
    echo_geomeans_matrix[i,j] = sqrt( (echo_abundance_matrix[i,j]/100) * echo_prevalence_matrix[i,j] )
end end

## Scaling Abundances
abvec = vcat(vec(khula_abundance_matrix), vec(cmd_abundance_matrix), vec(echo_abundance_matrix))
minab = minimum(abvec[abvec .!= 0.0])

full_abundance_matrix = log2.(full_abundance_matrix .+ minab/2)
cmd_abundance_matrix = log2.(cmd_abundance_matrix .+ minab/2)
khula_abundance_matrix = log2.(khula_abundance_matrix .+ minab/2)
echo_abundance_matrix = log2.(echo_abundance_matrix .+ minab/2)

## Perform HCA on the bug dimension and store the order
dist_taxa_abundances = pairwise(Euclidean(), full_abundance_matrix; dims=1)
dist_taxa_prevalences = pairwise(Euclidean(), full_prevalence_matrix; dims=1)
dist_taxa_geomeans = pairwise(Euclidean(), full_geomeans_matrix; dims=1)

# hcl_taxa_abundances = hclust(dist_taxa_abundances; linkage=:complete, branchorder=:optimal)
# hcl_taxa_prevalences = hclust(dist_taxa_prevalences; linkage=:complete, branchorder=:optimal)
# hcl_taxa_geomeans = hclust(dist_taxa_geomeans; linkage=:complete, branchorder=:optimal)

hcl_taxa_abundances = hclust(dist_taxa_abundances; linkage=:complete, branchorder=:optimal)
hcl_taxa_prevalences = hclust(dist_taxa_prevalences; linkage=:ward, branchorder=:optimal)
hcl_taxa_geomeans = hclust(dist_taxa_geomeans; linkage=:complete, branchorder=:optimal)

abundance_order = hcl_taxa_abundances.order
prevalence_order = hcl_taxa_prevalences.order
geomeans_order = hcl_taxa_geomeans.order
```

# Creating Master Figure 3
```julia
figure3_master = Figure(; size = (1600, 900))
```

## Plotting all the heatmaps

```julia
# processed_bugnames = replace.(important_bugs[prevalence_order], "_" => " ")
# split_bugnames = map(x -> split(x, " "), processed_bugnames)
# split_bugnames[17][2] = "sp"
# pushfirst!(split_bugnames[17], " ")
# pushfirst!(split_bugnames[28], " ")
# rejoined_bugnames = map( x -> (x[1][1] * ". " * join(x[2:end], " ")), split_bugnames)

rejoined_bugnames = replace.(important_bugs[prevalence_order], "_" => " ")

axA = Axis(
    figure3_master[1,1],
    xlabel = "Age bin (months)",
    title = "Samples from Baltic\nEST+FIN+RUS+SWE/DIABIMMUNE",
    xticks = (eachindex(bins), string.(bins)),
    yticks = (eachindex(important_bugs), rejoined_bugnames),
    yticklabelfont="TeX Gyre Heros Makie Italic",
    yticklabelsize=24,
    xticklabelsize=24,
    ylabelsize=24,
    xlabelsize=24,
    titlesize = 24,
    yreversed=true
)

# hmA = heatmap!(axA, cmd_prevalence_matrix[prevalence_order, :]', colormap = cgrad(:lapaz, rev = true))
hmA = heatmap!(axA, cmd_prevalence_matrix[prevalence_order, :]', colormap = cgrad(:lapaz, rev = false))

axB = Axis(
    figure3_master[1,2],
    xlabel = "Age bin (months)",
    ylabel = "species important for the age model",
    title = "Samples from North America\nUSA/ECHO+CMD",
    xticks = (eachindex(bins), string.(bins)),
    yticks = (eachindex(important_bugs), rejoined_bugnames),
    yticklabelfont="TeX Gyre Heros Makie Italic",
    yticklabelsize=24,
    xticklabelsize=24,
    ylabelsize=24,
    xlabelsize=24,
    titlesize = 24,
    yreversed=true
)
hideydecorations!(axB)

# hmB = heatmap!(axB, echo_prevalence_matrix[prevalence_order, :]', colormap = cgrad(:lapaz, rev = true))
hmB = heatmap!(axB, echo_prevalence_matrix[prevalence_order, :]', colormap = cgrad(:lapaz, rev = false))

axC = Axis(
    figure3_master[1,3],
    xlabel = "Age bin (months)",
    ylabel = "species important for the age model",
    title = "Samples from South Africa\nZAF/KHULA ",
    xticks = (eachindex(bins), string.(bins)),
    yticks = (eachindex(important_bugs), replace.(important_bugs[prevalence_order], "_" => " ")),
    yticklabelfont="TeX Gyre Heros Makie Italic",
    yticklabelsize=24,
    xticklabelsize=24,
    ylabelsize=24,
    xlabelsize=24,
    titlesize = 24,
    yreversed=true
)
hideydecorations!(axC)

# hmC = heatmap!(axC, khula_prevalence_matrix[prevalence_order, :]', colormap = cgrad(:lapaz, rev = true))
hmC = heatmap!(axC, khula_prevalence_matrix[prevalence_order, :]', colormap = cgrad(:lapaz, rev = false))

Colorbar(figure3_master[1,4], hmA, label = "Average species prevalence on each age range on the 0-1 scale", labelsize = 24)
```

## Add labels
```julia
Label(figure3_master[1, 1, TopLeft()], "A", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(figure3_master[1, 2, TopLeft()], "B", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(figure3_master[1, 3, TopLeft()], "C", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
```

## Fix layout
```julia
colsize!(figure3_master.layout, 1, Relative(0.33))
colsize!(figure3_master.layout, 2, Relative(0.33))
colsize!(figure3_master.layout, 3, Relative(0.33))
```

# Export Figure 2

```julia
save(joinpath(outdir, "figures", "Figure3.png"), figure3_master)
figure3_master
```