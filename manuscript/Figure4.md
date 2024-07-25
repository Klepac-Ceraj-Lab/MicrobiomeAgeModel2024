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
using Colors
```

### Configurable parameters
```julia
experiment_name = "2024AgeModelManuscript"
outdir = joinpath(pwd(), "results", experiment_name)
figdir = joinpath(outdir, "figures")
deepdivemonodir, deepdivecolordir = ( joinpath(figdir, "species_monocolor_scatterplots"), joinpath(figdir, "species_colored_scatterplots") )
```

## Helper functions
```julia
extremes(v::AbstractVector, n::Integer) = vcat(v[1:n], v[(end-n+1):end])

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

functional_mdata = @chain taxonomic_profiles begin
    khula_bugvec = subset(:datasource => x -> x .== "1kDLEAP-KHULA")
    select([:subject_id, :sample, :ageMonths, :visit])
end 

functional_profiles = myload(ECProfiles(); timepoint_metadata = functional_mdata) # this can take a bit
filtered_functional_profiles = relativeabundance(filter(f -> name(f) == "UNMAPPED" || hastaxon(f), functional_profiles))
ecs_mdata = DataFrame(get(filtered_functional_profiles))
# samples_inrename = ecs_mdata.sample
# samples_expected = functional_mdata.sample
# @assert (all(samples_expected .∈ Ref(samples_inrename)) & all(samples_inrename .∈ Ref(samples_expected)))
# setdiff(samples_expected, samples_inrename)
# setdiff(samples_inrename, samples_expected)

t1_samples = subset(ecs_mdata, :visit => x -> x .== "3mo")
t2_samples = subset(ecs_mdata, :visit => x -> x .== "6mo")
t3_samples = subset(ecs_mdata, :visit => x -> x .== "12mo")
longitudinal_samples = innerjoin(t1_samples, t3_samples, on = :subject_id, makeunique=true)
```

## Finding the important predictors
```julia
# @show sort(report_regression_merits(regression_Age_FullCV), :Val_RMSE_mean) # To check the nest hyperparameter index
hp_idx = 15

importances_table = @chain Leap.hpimportances(regression_Age_FullCV, hp_idx) begin
    subset(:variable => x -> x .!= "richness")
    subset(:variable => x -> x .!= "Shannon_index")
end

important_bugs = importances_table[1:30, :variable]
```

## Selecting samples and functions for functional analysis

Important note: the sum of all taxon-assigned abundances will result on the abundance without the taxon.
```julia
filtered_functions = union( name.(features(filtered_functional_profiles)) )
selected_functions_widedfs = [ comm2wide(filter(feat-> name((feat)) == this_function, filtered_functional_profiles)) for this_function in filtered_functions ]

scores = Vector{Float64}(undef, length(filtered_functions))
diff_MEANs = Vector{Float64}(undef, length(filtered_functions))
foldchanges = Vector{Float64}(undef, length(filtered_functions))
diff_STDs = Vector{Float64}(undef, length(filtered_functions))
diff_CIs = Vector{Float64}(undef, length(filtered_functions))

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

    scores[i] = round( (sum( ((allsamps.sum_1 .> 10) .| (allsamps.sum .> 10) ) .* sign.(allsamps.sum_1 .- allsamps.sum)) / nrow(allsamps)); digits = 2)

    diff_MEANs[i] = mean(allsamps.sum_1) - mean(allsamps.sum)
    foldchanges[i] = mean(log10.(allsamps.sum_1 .+ 1e-3)) - mean(log10.(allsamps.sum .+ 1e-3))
    diff_STDs[i] = Statistics.std((allsamps.sum_1) .- mean(allsamps.sum))
    diff_CIs[i] = ( 1.96 * Statistics.std((allsamps.sum_1) .- mean(allsamps.sum)) ) / sqrt(length(allsamps.sum_1))
end
```

We aimed to plot the 1.5% extremal ECs, or:
```julia
n_to_collect = ceil(Int64, (length(filtered_functions)*0.015)/2.0)
func_idxes = extremes(sortperm(scores), n_to_collect)
@show selected_functions = sort(filtered_functions[func_idxes])[1:end]

func_stats_df = DataFrame(
    :function_name => filtered_functions,
    :formatted_diff => [ "$(round(i; digits = 2)) ± $(round(j; digits = 2))" for (i,j) in zip(diff_MEANs,diff_CIs) ],
    :difference => diff_MEANs,
    :fold_change => foldchanges,
    :score => scores,
    :difference_CI => diff_CIs,
    :difference_SD => diff_STDs,    
)
sort!(func_stats_df, :score)
sort!(func_stats_df, :difference)
sort!(func_stats_df, :fold_change)

```

## Computing each taxa's contribution to each genefunction on each age range
```julia
youngsamplemat = zeros(length(important_bugs), length(selected_functions))
oldsamplemat = zeros(length(important_bugs), length(selected_functions))

for i in eachindex(important_bugs)
    for j in eachindex(selected_functions)

        try            
            this_filtered_profiles = filter(feat-> (hastaxon(feat) && (name(feat) == selected_functions[j]) && (name(taxon(feat)) == important_bugs[i])), functional_profiles)

            this_widetable = comm2wide(this_filtered_profiles)
            this_widetable.sum = map(x -> sum(x)*1e6, eachrow(Matrix(this_widetable[:, 5:end])))
            select!(this_widetable, [:sample, :subject_id, :ageMonths, :visit, :sum])


            youngsamps = subset(this_widetable, :visit => x -> x .== "3mo")
            # youngsamps = subset(this_widetable, :ageMonths => x -> x .<= 4.0)
            select!(youngsamps, [:subject_id, :sum])
            
            oldsamps = subset(this_widetable, :visit => x -> x .== "12mo")
            # oldsamps = subset(this_widetable, :ageMonths => x -> x .>= 9.0)
            select!(oldsamps, [:subject_id, :sum])
        
            allsamps = innerjoin(youngsamps, oldsamps, on = :subject_id; makeunique = true)

            # Fold change will be mean(log2(old)) - mean(log2(young))
            youngsamplemat[i,j] = mean(log2.(allsamps.sum .+ 1e-3))
            oldsamplemat[i,j] = mean(log2.(allsamps.sum_1 .+ 1e-3))

        catch
            continue
        end

    end
end

## Computing the optimal order for the functions and taxa via HCA

## Perform the actual HCA and store the order
dist_taxa = pairwise(Euclidean(), youngsamplemat + oldsamplemat; dims=1)
# dist_taxa = pairwise(Euclidean(), youngsamplemat - oldsamplemat; dims=1)
dist_functions = pairwise(Euclidean(), youngsamplemat + oldsamplemat; dims=2)
# dist_functions = pairwise(Euclidean(), youngsamplemat - oldsamplemat; dims=2)
hcl_taxa = hclust(dist_taxa; linkage=:single, branchorder=:optimal)
hcl_functions = hclust(dist_functions; linkage=:single, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_function_order = hcl_functions.order

ordered_youngsamplemat = youngsamplemat[hclust_taxa_order, hclust_function_order]
ordered_oldsamplemat = oldsamplemat[hclust_taxa_order, hclust_function_order]

ordered_diffmat = ordered_oldsamplemat - ordered_youngsamplemat

# a = Vector{Tuple{String, String, Float64}}()
a = Vector{NamedTuple{(:species, :ec, :diff), Tuple{String, String, Float64}}}()

for i in eachindex(hclust_taxa_order)
    for j in eachindex(hclust_function_order)
        push!(a, (species = important_bugs[hclust_taxa_order][i], ec =  selected_functions[hclust_function_order][j], diff = ordered_diffmat[i,j] ))
    end
end

a = DataFrame(a)
sort!(a, :diff)

ordered_functions = selected_functions[hclust_function_order]
ordered_taxa = important_bugs[hclust_taxa_order]

## For paper results paragraph

absolute_differences_mat = abs.(ordered_oldsamplemat .- ordered_youngsamplemat)

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

# subset_to_plot = vcat(collect(1:8), collect(19:30))
subset_to_plot = collect(1:30)
```

# Creating Master Figure 4
```julia
figure4_master = Figure(; size = (1500, 1400))
```

## Figure 2, Panels A and B - Functional Heatmaps
```julia
axA = Axis(
    figure4_master[1, 1],
    xticks = (1:length(ordered_taxa[subset_to_plot]), [reformat_taxa(s, importances_table) for s in ordered_taxa[subset_to_plot] ]),
    yticks = (1:length(ordered_functions), [ reformat_ecs(el, ec_colors) for el in ordered_functions ]),
    xticklabelrotation = pi/2,
    xticklabelsize = 18,
    yticklabelsize = 18,
    yreversed = true,
    xticklabelfont = "TeX Gyre Heros Makie Italic",
    title = "3 months timepoint", titlesize = 20)
hideydecorations!(axA)
axB = Axis(
    figure4_master[1, 2],
    xticks = (1:length(ordered_taxa[subset_to_plot]), [reformat_taxa(s, importances_table) for s in ordered_taxa[subset_to_plot] ]),
    yticks = (1:length(ordered_functions), [ reformat_ecs(el, ec_colors) for el in ordered_functions ]),
    xticklabelrotation = pi/2,
    xticklabelsize = 18,
    yticklabelsize = 18,
    yaxisposition = :right,
    xticklabelfont = "TeX Gyre Heros Makie Italic",
    yreversed = true,
    title = "12 months timepoint", titlesize = 20)

## Plot with subsets to make the figure slightly smaller and simpler
hm = heatmap!(axA, ordered_youngsamplemat[subset_to_plot, :], colormap = :magma, colorrange = (0.0, 3.0))
hm = heatmap!(axB, ordered_oldsamplemat[subset_to_plot, :], colormap = :magma, colorrange = (0.0, 3.0))

Legend(
    figure4_master[2, 1],
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

Colorbar(figure4_master[2,2], hm, label = "log10(CPM)", vertical = false)

save(joinpath(outdir, "figures", "Figure4.png"), figure4_master)
save(joinpath(outdir, "figures", "Figure4.eps"), figure4_master)
save(joinpath(outdir, "figures", "Figure4.svg"), figure4_master)
figure4_master
```