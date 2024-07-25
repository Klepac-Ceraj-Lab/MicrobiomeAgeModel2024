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
CSV.write("ecs_metadata.csv", ecs_mdata)
# samples_inrename = ecs_mdata.sample
# samples_expected = functional_mdata.sample
# @assert (all(samples_expected .∈ Ref(samples_inrename)) & all(samples_inrename .∈ Ref(samples_expected)))
# setdiff(samples_expected, samples_inrename)
# setdiff(samples_inrename, samples_expected)

t1_samples = subset(ecs_mdata, :visit => x -> x .== "3mo")
println("For t1 samples, M = $(round(mean(t1_samples.ageMonths); digits = 2)), SD = $(round(std(t1_samples.ageMonths); digits = 2)), range = $(extrema(t1_samples.ageMonths))")
t2_samples = subset(ecs_mdata, :visit => x -> x .== "6mo")
println("For t2 samples, M = $(round(mean(t2_samples.ageMonths); digits = 2)), SD = $(round(std(t2_samples.ageMonths); digits = 2)), range = $(extrema(t2_samples.ageMonths))")
t3_samples = subset(ecs_mdata, :visit => x -> x .== "12mo")
println("For t3 samples, M = $(round(mean(t3_samples.ageMonths); digits = 2)), SD = $(round(std(t3_samples.ageMonths); digits = 2)), range = $(extrema(t3_samples.ageMonths))")

longitudinal_samples = innerjoin(t1_samples, t3_samples, on = :subject_id, makeunique=true)
```

## Finding the important predictors
```julia
# @show sort(report_regression_merits(regression_Age_FullCV), :Val_RMSE_mean) # To check the nest hyperparameter index
hp_idx = 15

importances_table = hpimportances(regression_Age_FullCV, hp_idx)
importances_table.cumsum = cumsum(importances_table.weightedImportance)
nfeat_toplot = findfirst(importances_table.cumsum .> 0.7)
importances_table = importances_table[1:nfeat_toplot,:]

importances_table.correl = [ cor(taxonomic_profiles[:, ccol], taxonomic_profiles.ageMonths) for ccol in importances_table.variable ]
importances_table.impsign = importances_table.weightedImportance .* sign.(importances_table.correl)

onlyspecies_importances = subset(importances_table, :variable => x -> x .!= "Shannon_index")

important_bugs = onlyspecies_importances.variable
```

## Selecting samples and functions for functional analysis

Important note: the sum of all taxon-assigned abundances will result on the abundance without the taxon.
```julia
# filtered_functions = union( first.(split.(name.(features(filtered_functional_profiles)), '|')))
filtered_functions = union( name.(features(filtered_functional_profiles)))

func_stats_df = map(enumerate(filtered_functions)) do (i,f)
    i % 10 == 0 && @info i
    f == "UNGROUPED" && return (; function_name = f, score = NaN, difference = NaN, fold_change = NaN, difference_SD = NaN, difference_CI = NaN)
    f == "UNMAPPED" && return (; function_name = f, score = NaN, difference = NaN, fold_change = NaN, difference_SD = NaN, difference_CI = NaN)

    fnum = split(first(split(f, ':')), '.')
    srch = Regex("^" * join(fnum, raw"\.") * ":")

    df = comm2wide(filtered_functional_profiles[srch, :])
    df.sum = map(x -> sum(x)*1e6, eachrow(Matrix(df[:, 5:end])))

    select!(df, [:sample, :subject_id, :ageMonths, :visit, :sum])
    youngsamps = subset(df, :visit => x -> x .== "3mo")
    # youngsamps = subset(df, :ageMonths => x -> x .<= 4.0)
    rename!(youngsamps, :sum => :sum_young )
    select!(youngsamps, [:subject_id, :sum_young])

    oldsamps = subset(df, :visit => x -> x .== "12mo")
    # oldsamps = subset(df, :ageMonths => x -> x .>= 9.0)
    rename!(oldsamps, :sum => :sum_old )
    select!(oldsamps, [:subject_id, :sum_old])

    allsamps = innerjoin(youngsamps, oldsamps, on = :subject_id; makeunique = true)

    score = round( (sum( ((allsamps.sum_old .> 10) .| (allsamps.sum_young .> 10) ) .* sign.(allsamps.sum_old .- allsamps.sum_young)) / nrow(allsamps)); digits = 2)

    difference = mean(allsamps.sum_old) - mean(allsamps.sum_young)
    fold_change = mean(log10.(allsamps.sum_old .+ 1e-3)) - mean(log10.(allsamps.sum_young .+ 1e-3))
    difference_SD = Statistics.std((allsamps.sum_old) .- mean(allsamps.sum_young))
    difference_CI = ( 1.96 * Statistics.std((allsamps.sum_old) .- mean(allsamps.sum_young)) ) / sqrt(length(allsamps.sum_old))
    return (; function_name = f, score, difference, fold_change, difference_SD, difference_CI)
end |> DataFrame
```

We aimed to plot the 1.5% extremal ECs, or:
```julia
subset!(func_stats_df, :score => ByRow(x -> !isnan(x)))
sort!(func_stats_df, :score)
# sort!(func_stats_df, :difference)
# sort!(func_stats_df, :fold_change)

n_to_collect = ceil(Int64, (nrow(func_stats_df)*0.015)/2.0)
selected_functions = extremes(func_stats_df.function_name, n_to_collect)
```

## Computing each taxa's contribution to each genefunction on each age range
```julia
youngsamplemat = zeros(length(important_bugs), length(selected_functions))
oldsamplemat = zeros(length(important_bugs), length(selected_functions))

for i in eachindex(important_bugs)
    for j in eachindex(selected_functions)

        try            
            this_filtered_profiles = filter(feat-> (hastaxon(feat) && (name(feat) == selected_functions[j]) && (name(taxon(feat)) == important_bugs[i])), filtered_functional_profiles)

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
            # youngsamplemat[i,j] = mean(log10.(allsamps.sum .+ 1e-3))
            # oldsamplemat[i,j] = mean(log10.(allsamps.sum_1 .+ 1e-3))
            # youngsamplemat[i,j] = log10(mean(allsamps.sum) + 1e-3)
            # oldsamplemat[i,j] = log10(mean(allsamps.sum_1) + 1e-3)
            youngsamplemat[i,j] = mean(allsamps.sum)
            oldsamplemat[i,j] = mean(allsamps.sum_1)
        catch
            youngsamplemat[i,j] = 0.0
            oldsamplemat[i,j] = 0.0
            continue
        end

    end
end
```

## Computing the optimal order for the functions and taxa via HCA
```julia
## Perform the actual HCA and store the order
dist_taxa = pairwise(Euclidean(), youngsamplemat - oldsamplemat; dims=1)
dist_functions = pairwise(Euclidean(), youngsamplemat - oldsamplemat; dims=2)
hcl_taxa = hclust(dist_taxa; linkage=:ward, branchorder=:optimal)
hcl_functions = hclust(dist_functions; linkage=:ward, branchorder=:optimal)
hclust_taxa_order = hcl_taxa.order
hclust_function_order = hcl_functions.order

ordered_youngsamplemat = log10.(10e-3 .+ youngsamplemat[hclust_taxa_order, hclust_function_order])
ordered_oldsamplemat = log10.(10e-3 .+ oldsamplemat[hclust_taxa_order, hclust_function_order])

ordered_diffmat = ordered_oldsamplemat - ordered_youngsamplemat

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
```

## Some prints for paper results paragraph
```julia
vatanen2018_ECs = [
    "3.6.1.1: Inorganic diphosphatase",
    "2.7.7.56: tRNA nucleotidyltransferase",
    "2.6.1.42: Branched-chain-amino-acid transaminase",
    "3.1.22.4: Crossover junction endodeoxyribonuclease",
    "1.6.1.2: NAD(P)(+) transhydrogenase (Re/Si-specific)",
    "1.1.1.44: Phosphogluconate dehydrogenase (NADP(+)-dependent, decarboxylating)",
    "5.4.2.11: Phosphoglycerate mutase (2,3-diphosphoglycerate-dependent)",
    "6.3.4.18: 5-(carboxyamino)imidazole ribonucleotide synthase",
    "6.3.1.2: Glutamate--ammonia ligase",
    "2.3.1.117: 2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-succinyltransferase",
    "5.3.1.6: Ribose-5-phosphate isomerase",
    "4.2.1.1: Carbonate dehydratase",
    "2.7.7.72: CCA tRNA nucleotidyltransferase",
    "2.5.1.74: 1,4-dihydroxy-2-naphthoate polyprenyltransferase",
    "2.3.1.54: Formate C-acetyltransferase",
    "4.4.1.8: Cystathionine beta-lyase",
    "2.7.1.15: Ribokinase",
    "1.11.1.15: Peroxiredoxin",
    "6.3.4.14: Biotin carboxylase",
    "1.1.1.27: L-lactate dehydrogenase",
    "6.1.1.18: Glutamine--tRNA ligase",
    "2.8.4.4: [Ribosomal protein S12] (aspartate(89)-C(3))-methylthiotransferase",
    "2.2.1.1: Transketolase",
    "2.7.1.144: Tagatose-6-phosphate kinase",
    "2.3.1.179: Beta-ketoacyl-[acyl-carrier-protein] synthase II",
    "5.4.2.12: Phosphoglycerate mutase (2,3-diphosphoglycerate-independent)",
    "2.7.1.26: Riboflavin kinase",
    "2.7.1.11: 6-phosphofructokinase",
    "4.1.1.49: Phosphoenolpyruvate carboxykinase (ATP)",
    "2.6.1.83: LL-diaminopimelate aminotransferase",
    "2.3.1.51: 1-acylglycerol-3-phosphate O-acyltransferase",
    "1.7.99.1: Hydroxylamine reductase",
    "2.4.1.21: Starch synthase",
    "2.5.1.3: Thiamine-phosphate diphosphorylase",
    "4.1.99.17: Phosphomethylpyrimidine synthase",
    "2.7.1.90: Diphosphate--fructose-6-phosphate 1-phosphotransferase",
    "3.4.24.78: GPR endopeptidase",
    "2.5.1.49: O-acetylhomoserine aminocarboxypropyltransferase",
    "1.4.1.14: Glutamate synthase (NADH)",
    "4.1.1.3: Oxaloacetate decarboxylase"
]

absolute_differences_mat = abs.(ordered_oldsamplemat .- ordered_youngsamplemat)
@show func_stats_df

println("Number of functions from Vatanen2018 on our list: $(sum(ordered_functions .∈ Ref(vatanen2018_ECs))) or $(round(100*sum(ordered_functions .∈ Ref(vatanen2018_ECs))/length(vatanen2018_ECs); digits = 2))")

youngclustbugs = [ "Bifidobacterium_longum", "Bifidobacterium_breve", "Escherichia_coli", "Ruminococcus_gnavus" ]
oldclustbugs = ["Dorea_longicatena", "Blautia_obeum", "Blautia_wexlerae", "Anaerostipes_hadrus", "Faecalibacterium_prausnitzii", "Prevotella_copri"]

youngclust_trim = ordered_diffmat[ordered_taxa .∈ Ref(youngclustbugs), 1:n_to_collect]
oldclust_trim = ordered_diffmat[ordered_taxa .∈ Ref(oldclustbugs), (n_to_collect+1):(2*n_to_collect)]

println("For the cluster of species decreasing in abundance, mean fold change of EC abundance is $(mean(youngclust_trim[youngclust_trim .!= 0.0])) +- $(Leap.conf_interval(youngclust_trim[youngclust_trim .!= 0.0]))")
println("For the cluster of species increasing in abundance, mean fold change of EC abundance is $(mean(oldclust_trim[oldclust_trim .!= 0.0])) +- $(Leap.conf_interval(oldclust_trim[oldclust_trim .!= 0.0]))")
#####
# Building Figure
#####

function reformat_taxa(ttaxa::String, imptab::DataFrame)

    sn = sign(subset(imptab, :variable => x -> x .== ttaxa).correl[1])
    ttaxa = replace(ttaxa, "_" => " ")
    if sn == +1.0
        return(rich(ttaxa, "   ", rich("▶ "; color = :blue)))
    else
        return(rich(ttaxa, "   ", rich("◀ "; color = :red)))
    end
end

function reformat_ecs(eecs::String, colset::Dict; pattern::Regex = r"(\d+\.\d+\.\d+\.?\d*):\s*(.*)")

    humann_ec_matches = eachmatch(pattern, eecs)

    ec_num = string(first(humann_ec_matches).captures[1])
    ec_fun = string(first(humann_ec_matches).captures[2])

    rich(rich("█\t"; font = :bold, color = colset[ec_num[1]]), ec_fun, rich(" [$(ec_num)]"; font = :bold, color = :black))
    # rich(rich("█\t"; font = :bold, color = colset[ec_num[1]]), ec_fun, rich(" [$(ec_num)]"; font = :bold, color = colset[ec_num[1]]))
    # Other examples of unicode rectangles: ▮ (version 1); █ ▉ ▊ ▋ ▌ ▍ ▎ ▏

end

ec_colors = Dict(
    '1' => distinguishable_colors(30)[7],
    '2' => distinguishable_colors(30)[8],
    '3' => distinguishable_colors(30)[10],
    '4' => distinguishable_colors(30)[12],
    '5' => distinguishable_colors(30)[21],
    '6' => distinguishable_colors(30)[24],
    '7' => distinguishable_colors(30)[24]
)
```

## Selecting subset of Taxa to plot

```julia
manual_taxa_to_plot = [
 "Bifidobacterium_longum",
 "Bifidobacterium_breve",
 "Escherichia_coli",
 "Ruminococcus_gnavus",
 "Dorea_longicatena",
 "Roseburia_inulinivorans",
 "Eubacterium_rectale",
 "Ruminococcus_bromii",
 "Haemophilus_parainfluenzae",
 "Bacteroides_fragilis",
#  "Enterococcus_faecalis",
#  "Hungatella_hathewayi",
#  "Flavonifractor_plautii",
#  "Eggerthella_lenta",
#  "Intestinibacter_bartlettii",
#  "Firmicutes_bacterium_CAG_41",
#  "Roseburia_intestinalis",
#  "Clostridium_symbiosum",
#  "Clostridium_sp_AM22_11AC",
#  "Eubacterium_eligens",
#  "Agathobaculum_butyriciproducens",
#  "Fusicatenibacter_saccharivorans",
#  "Roseburia_faecis",
#  "Eubacterium_hallii",
 "Clostridium_innocuum",
 "Dorea_formicigenerans",
 "Erysipelatoclostridium_ramosum",
 "Clostridium_neonatale",
 "Streptococcus_thermophilus",
 "Blautia_obeum",
 "Blautia_wexlerae",
 "Anaerostipes_hadrus",
 "Faecalibacterium_prausnitzii",
 "Prevotella_copri"
]

# subset_to_plot = collect(1:nfeat_toplot-1) # For debugging purposes, plot all species
subset_taxa_plot = ordered_taxa .∈ Ref(manual_taxa_to_plot)
subset_function_plot = 1:length(ordered_functions)-1
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