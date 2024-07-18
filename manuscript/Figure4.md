# Figure 4 - Functional analysis heatmaps and top changing functions

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

## Helper functions
```julia
extremes(v::AbstractVector, n::Integer) = vcat(v[1:n], v[end-n:end])

function myfeaturefunc(s::String)
    s = replace(s, r"\|g__\w+\."=>"|")
    genefunction(s)
end

function myload(::ECProfiles; timepoint_metadata = load(Metadata()))
    comm = Leap.read_arrow(inputfiles("ecs.arrow"); featurefunc = myfeaturefunc)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end
```

## Loading data
```julia
JLD2.@load joinpath(outdir, "AgeModel_FullCV_Results.jld") regression_Age_FullCV
taxonomic_profiles = regression_Age_FullCV.original_df
subset!(taxonomic_profiles, :subject_id => x -> x .!= "khula-191-63214792")

functional_mdata = @chain taxonomic_profiles begin
    khula_bugvec = subset(:datasource => x -> x .== "1kDLEAP-KHULA")
    select([:subject_id, :sample, :ageMonths, :visit])
end 

functional_profiles = myload(ECProfiles(); timepoint_metadata = functional_mdata) # this can take a bit
relativeabundance!(functional_profiles)
ecs_mdata = DataFrame(get(functional_profiles))
samples_inrename = ecs_mdata.sample
samples_expected = functional_mdata.sample
@assert (all(samples_expected .∈ Ref(samples_inrename)) & all(samples_inrename .∈ Ref(samples_expected)))
setdiff(samples_expected, samples_inrename)

t1_samples = subset(ecs_mdata, :visit => x -> x .== "3mo")
t2_samples = subset(ecs_mdata, :visit => x -> x .== "6mo")
t3_samples = subset(ecs_mdata, :visit => x -> x .== "12mo")

longitudinal_samples = innerjoin(t1_samples, t3_samples, on = :subject_id, makeunique=true)
subset!(longitudinal_samples, :subject_id => x -> x .!= "khula-191-63214792")
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

## Selecting samples and functions for functional analysis
```julia
filtered_profiles = filter(feat-> hastaxon(feat), functional_profiles)
filtered_functions = union( name.(features(filtered_profiles)) )
selected_functions_widedfs = [ comm2wide(filter(feat-> name((feat)) == this_function, filtered_profiles)) for this_function in filtered_functions ]

# scores = Vector{Int64}(undef, length(filtered_functions))
scores = Vector{Float64}(undef, length(filtered_functions))
diffs = Vector{Float64}(undef, length(filtered_functions))

for i in eachindex(selected_functions_widedfs)
    
    selected_functions_widedfs[i].sum = map(x -> sum(x)*1e6, eachrow(Matrix(selected_functions_widedfs[i][:, 5:end])))
    
    select!(selected_functions_widedfs[i], [:sample, :subject_id, :ageMonths, :visit, :sum])
    youngsamps = subset(selected_functions_widedfs[i], :visit => x -> x .== "3mo")
    # youngsamps = subset(selected_functions_widedfs[i], :ageMonths => x -> x .<= 4.0)

    select!(youngsamps, [:subject_id, :sum])
    oldsamps = subset(selected_functions_widedfs[i], :visit => x -> x .== "12mo")
    # oldsamps = subset(selected_functions_widedfs[i], :ageMonths => x -> x .>= 9.0)

    select!(oldsamps, [:subject_id, :sum])
    @show allsamps = innerjoin(youngsamps, oldsamps, on = :subject_id; makeunique = true)

    scores[i] = sum( ((allsamps.sum_1 .> 10) .| (allsamps.sum .> 10) ) .& (allsamps.sum_1 .> allsamps.sum)) - sum( ((allsamps.sum_1 .> 10) .| (allsamps.sum .> 10) ) .& (allsamps.sum .> allsamps.sum_1)) # This is really good!
    # scores[i] = sum( ((allsamps.sum_1 .> 100) .| (allsamps.sum .> 100) ) .& (allsamps.sum_1 .> allsamps.sum)) - sum( ((allsamps.sum_1 .> 100) .| (allsamps.sum .> 100) ) .& (allsamps.sum .> allsamps.sum_1)) # This is really good!
    diffs[i] = mean(allsamps.sum_1) - mean(allsamps.sum)
end

func_idxes = extremes(sortperm(scores), 20)
selected_functions = sort(filtered_functions[func_idxes])[1:end-1]

#####
# Computing each taxa's contribution to each genefunction on each age range
#####
youngsamplemat = zeros(length(important_bugs), length(selected_functions))
oldsamplemat = zeros(length(important_bugs), length(selected_functions))

for i in eachindex(important_bugs)
    for j in eachindex(selected_functions)

        try            
            this_filtered_profiles = filter(feat-> (hastaxon(feat) && (name(feat) == selected_functions[j]) && (name(taxon(feat)) == important_bugs[i])), functional_profiles)

            this_widetable = comm2wide(this_filtered_profiles)
            this_widetable.sum = map(x -> sum(x)*1e7, eachrow(Matrix(this_widetable[:, 5:end])))
            select!(this_widetable, [:sample, :subject_id, :ageMonths, :visit, :sum])


            youngsamps = subset(this_widetable, :ageMonths => x -> x .<= 4.0)
            select!(youngsamps, [:subject_id, :sum])
            
            oldsamps = subset(this_widetable, :ageMonths => x -> x .>= 9.0)
            select!(oldsamps, [:subject_id, :sum])
        
            allsamps = innerjoin(youngsamps, oldsamps, on = :subject_id; makeunique = true)

            youngsamplemat[i,j] = maximum([log10(mean(allsamps.sum .+ 1e-3)), 0.0])
            oldsamplemat[i,j] = maximum([log10(mean(allsamps.sum_1 .+ 1e-3)), 0.0])

        catch
            continue
        end

    end
end

## Computing the optimal order for the functions and taxa via HCA

## Perform the actual HCA and store the order
dist_taxa = pairwise(Euclidean(), youngsamplemat + oldsamplemat; dims=1)
dist_functions = pairwise(Euclidean(), youngsamplemat + oldsamplemat; dims=2)
hcl_taxa = hclust(dist_taxa; linkage=:single, branchorder=:optimal)
hcl_functions = hclust(dist_functions; linkage=:single, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_function_order = hcl_functions.order

ordered_youngsamplemat = youngsamplemat[hclust_taxa_order, hclust_function_order]
ordered_oldsamplemat = oldsamplemat[hclust_taxa_order, hclust_function_order]

ordered_functions = selected_functions[hclust_function_order]
ordered_taxa = important_bugs[hclust_taxa_order]

## For paper results paragraph

absolute_differences_mat = ordered_oldsamplemat .- ordered_youngsamplemat

fourbugclust = absolute_differences_mat[27:30, 1:20]
@show(mean(fourbugclust[fourbugclust .!= 0.0]))
@show(conf_interval(fourbugclust[fourbugclust .!= 0.0]))

bigbugclust = absolute_differences_mat[19:26,21:40]
@show(mean(bigbugclust[bigbugclust .!= 0.0]))
@show(conf_interval(bigbugclust[bigbugclust .!= 0.0]))

#####
# Building Figure
#####

importances_table.correl = [ cor(taxonomic_profiles[:, ccol], taxonomic_profiles.ageMonths) for ccol in importances_table.variable ]

function reformat_taxa(ttaxa::String, imptab::DataFrame)

    sn = sign(subset(imptab, :variable => x -> x .== ttaxa).correl[1])
    ttaxa = replace(ttaxa, "_" => " ")
    if sn == +1.0
        return(rich(ttaxa, "   ", rich("▶"; color = :blue)))
    else
        return(rich(ttaxa, "   ", rich("◀"; color = :red)))
    end
end

function reformat_ecs(eecs::String, colset::Dict; pattern::Regex = r"(\d+\.\d+\.\d+\.?\d*):\s*(.*)")

    humann_ec_matches = eachmatch(pattern, eecs)

    ec_num = string(first(humann_ec_matches).captures[1])
    ec_fun = string(first(humann_ec_matches).captures[2])

    rich(rich("█\t"; font = :bold, color = colset[ec_num[1]]), ec_fun, rich(" [$(ec_num)]"; font = :bold, color = colset[ec_num[1]]))
    # Other examples of unicode rectangles: ▮ (version 1); █ ▉ ▊ ▋ ▌ ▍ ▎ ▏

end

ec_colors = Dict(
    '1' => "darkgreen",
    '2' => "chartreuse3",
    '3' => "orangered",
    '4' => "orange",
    '5' => "maroon",
    '6' => "midnightblue"
)

subset_to_plot = vcat(collect(1:8), collect(19:30))
fig = Figure(; size = (1500, 1400))

ax1 = Axis(
    fig[1, 1],
    xticks = (1:length(ordered_taxa[subset_to_plot]), [reformat_taxa(s, importances_table) for s in ordered_taxa[subset_to_plot] ]),
    yticks = (1:length(ordered_functions), [ reformat_ecs(el, ec_colors) for el in ordered_functions ]),
    xticklabelrotation = pi/2,
    xticklabelsize = 18,
    yticklabelsize = 18,
    yreversed = true,
    xticklabelfont = "TeX Gyre Heros Makie Italic",
    title = "Samples < 4 months old", titlesize = 20)
hideydecorations!(ax1)
ax2 = Axis(
    fig[1, 2],
    xticks = (1:length(ordered_taxa[subset_to_plot]), [reformat_taxa(s, importances_table) for s in ordered_taxa[subset_to_plot] ]),
    yticks = (1:length(ordered_functions), [ reformat_ecs(el, ec_colors) for el in ordered_functions ]),
    xticklabelrotation = pi/2,
    xticklabelsize = 18,
    yticklabelsize = 18,
    yaxisposition = :right,
    xticklabelfont = "TeX Gyre Heros Makie Italic",
    yreversed = true,
    title = "Samples > 9 months old", titlesize = 20)

## Plot with subsets to make the figure slightly smaller and simpler
hm = heatmap!(ax1, ordered_youngsamplemat[subset_to_plot, :], colormap = :magma, colorrange = (0.0, 3.0))
hm = heatmap!(ax2, ordered_oldsamplemat[subset_to_plot, :], colormap = :magma, colorrange = (0.0, 3.0))

Legend(
    fig[2, 1],
    [
        MarkerElement(marker = :circle, color = ec_colors['1'], markersize = 14),
        MarkerElement(marker = :circle, color = ec_colors['2'], markersize = 14),
        MarkerElement(marker = :circle, color = ec_colors['3'], markersize = 14),
        MarkerElement(marker = :circle, color = ec_colors['4'], markersize = 14),
        MarkerElement(marker = :circle, color = ec_colors['5'], markersize = 14),
        MarkerElement(marker = :circle, color = ec_colors['6'], markersize = 14)
        ],
    [
        "EC 1. Oxidoreductases",
        "EC 2. Transferases",
        "EC 3. Hydrolases",
        "EC 4. Lyases",
        "EC 5. Isomerases",
        "EC 6. Ligases"
    ],
    orientation = :vertical,
    nbanks = 2,
    tellheight = true,
    tellwidth = false
)

Colorbar(fig[2,2], hm, label = "log10(CPM)", vertical = false)

fig

save(joinpath(outdir, "figures", "functionalHeatmap.png"), fig)
save(joinpath(outdir, "figures", "functionalHeatmap.eps"), fig)
```