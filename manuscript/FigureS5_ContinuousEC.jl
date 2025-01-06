# Figure 4 - Functional analysis heatmaps and top changing functions

## Pre-configuration

### Loading Packages
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
using CSV
using DataToolkit
using MixedModels
using MultipleTesting
using MicrobiomeAgeModel2024
Random.seed!(0)

### Configurable parameters and notebook set-up
outdir, figdir, deepdivemonodir, deepdivecolordir = setup_outdir(; experiment_name = "2024AgeModelRevisions")
presence_absence = false # This argument controls whether the analysis will be based on continous relative abundances or binary presence/absence of species.

#### UNCOMMENT ONLY ONE OF THE FOLLOWING 3 LINES TO PICK A SOURCE FOR THE ANALYSIS DATA
DataToolkit.loadcollection!("./Data_Local.toml")    ## Uncomment this line to use local files located on the "data" subfolder and the Local relative filesystem references
# DataToolkit.loadcollection!("./Data_AWS.toml")      ## Uncomment this line to use the datasets made available on the public AWS bucket
# DataToolkit.loadcollection!("./Data_Dryad.toml")    ## Uncomment this line to use the datasets published to Data Dryad (DOI: 10.5061/dryad.dbrv15f9z)

## Helper functions
extremes(v::AbstractVector, n::Integer) = vcat(v[1:n], v[(end-n+1):end])

function myfeaturefunc(s::String)
    s = replace(s, r"\|g__\w+\."=>"|")
    genefunction(s)
end

function myload(::ECProfiles, mypath; timepoint_metadata = load(Metadata()))
    comm = Leap.read_arrow(mypath; featurefunc = myfeaturefunc)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end

## Loading data
regression_Age_FullCV = d"cv_results"["regression_Age_FullCV"]
taxonomic_profiles = regression_Age_FullCV.original_df

bins = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

## Loading data

functional_mdata = d"ecs_metadata"
functional_profiles = myload(ECProfiles(), open(dataset("ec_profiles"), DataToolkit.FilePath).path; timepoint_metadata = functional_mdata) # this can take a bit
filtered_functional_profiles = relativeabundance(filter(f -> name(f) == "UNMAPPED" || hastaxon(f), functional_profiles))
ecs_mdata = DataFrame(get(filtered_functional_profiles))
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

longitudinal_samples = innerjoin(t1_samples, t2_samples, on = :subject_id, makeunique=true)
longitudinal_samples = innerjoin(longitudinal_samples, t3_samples, on = :subject_id, makeunique=true)
long_longitudinal_samples = vcat(
    select(longitudinal_samples, [:sample, :subject_id, :ageMonths, :visit]),
    rename(select(longitudinal_samples, [:sample_1, :subject_id, :ageMonths_1, :visit_1]), :sample_1 => :sample, :ageMonths_1 => :ageMonths, :visit_1 => :visit),
    rename(select(longitudinal_samples, [:sample_2, :subject_id, :ageMonths_2, :visit_2]), :sample_2 => :sample, :ageMonths_2 => :ageMonths, :visit_2 => :visit),
)

## Selecting samples and functions for functional analysis
# filtered_functions = union( first.(split.(name.(features(filtered_functional_profiles)), '|')))
filtered_functions = union( name.(features(filtered_functional_profiles)))

func_stats_df = map(enumerate(filtered_functions)) do (i,f)
    i % 10 == 0 && @info i
    f == "UNGROUPED" && return (; function_name = f, score = NaN, pvalue = NaN)
    f == "UNMAPPED" && return (; function_name = f, score = NaN, pvalue = NaN)

    fnum = split(first(split(f, ':')), '.')
    srch = Regex("^" * join(fnum, raw"\.") * ":")

    df = comm2wide(filtered_functional_profiles[srch, :])
    df.sum = map(x -> sum(x)*1e6, eachrow(Matrix(df[:, 5:end])))
    # print(df.sum)
    std(df.sum) == 0.0 && return (; function_name = f, score = NaN, pvalue = NaN)

    select!(df, [:sample, :sum])
    allsamps = innerjoin(long_longitudinal_samples, df, on = :sample)

    glmform = @formula( sum ~ ageMonths + (1 | subject_id) )

    println(f)
    fitted_line = LinearMixedModel(glmform, allsamps)
    try
        fit!(fitted_line)
    catch
        return (; function_name = f, score = NaN, pvalue = NaN)
    end

    fitresult = DataFrames.Tables.rowtable(coeftable(fitted_line))[2]

    score = fitresult.var"Coef."
    pvalue = fitresult.var"Pr(>|z|)"
    
    return (; function_name = f, score = score, pvalue = pvalue)
end |> DataFrame

subset!(func_stats_df, :score => ByRow(x -> !isnan(x)))

sort!(func_stats_df, :pvalue)

func_stats_df.qvalue = adjust(func_stats_df.pvalue, BenjaminiHochberg())
labelpval = function(pval)
    if pval < 0.001
        return("***")
    elseif pval < 0.01
        return("**")
    elseif pval < 0.05
        return("**")
    elseif pval < 0.1
        return("·")
    else
        return(" ")
    end
end
func_stats_df.label = map(labelpval, func_stats_df.pvalue)

sort!(func_stats_df, :score)

plot_df = func_stats_df[extremes(1:nrow(func_stats_df), 23), :]

## Some prints for paper results paragraph
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

figure4_ecs = [
    "2.7.1.15: Ribokinase",
    "1.17.4.1: Ribonucleoside-diphosphate reductase",
    "1.1.1.3: Homoserine dehydrogenase",
    "1.1.1.1: Alcohol dehydrogenase",
    "6.2.1.5: Succinate--CoA ligase (ADP-forming)",
    "3.6.4.13: RNA helicase",
    "2.2.1.2: Transaldolase",
    "2.7.1.35: Pyridoxal kinase",
    "1.2.1.12: Glyceraldehyde-3-phosphate dehydrogenase (phosphorylating)",
    "1.3.5.2: Dihydroorotate dehydrogenase (quinone)",
    "1.1.1.49: Glucose-6-phosphate dehydrogenase (NADP(+))",
    "4.6.1.1: Adenylate cyclase",
    "3.4.11.2: Membrane alanyl aminopeptidase",
    "2.3.1.117: 2,3,4,5-tetrahydropyridine-2,6-dicarboxylate N-succinyltransferase",
    "3.1.1.5: Lysophospholipase",
    "2.7.7.56: tRNA nucleotidyltransferase",
    "1.6.99.3: NADH dehydrogenase",
    "2.1.1.192: 23S rRNA (adenine(2503)-C(2))-methyltransferase",
    "4.2.3.1: Threonine synthase",
    "2.1.1.14: 5-methyltetrahydropteroyltriglutamate--homocysteine S-methyltransferase",
    "5.4.2.1: Transferred entry: 5.4.2.11 and 5.4.2.12",
    "1.4.1.14: Glutamate synthase (NADH)",
    "1.3.8.1: Short-chain acyl-CoA dehydrogenase",
    "4.2.1.55: 3-hydroxybutyryl-CoA dehydratase",
    "2.1.1.148: Thymidylate synthase (FAD)",
    "1.3.99.22: Coproporphyrinogen dehydrogenase",
    "2.8.4.4: [Ribosomal protein S12] (aspartate(89)-C(3))-methylthiotransferase",
    "2.4.99.17: S-adenosylmethionine:tRNA ribosyltransferase-isomerase",
    "1.5.1.7: Saccharopine dehydrogenase (NAD(+), L-lysine-forming)",
    "3.1.26.11: Ribonuclease Z",
    "5.4.2.12: Phosphoglycerate mutase (2,3-diphosphoglycerate-independent)",
    "3.4.11.4: Tripeptide aminopeptidase",
    "2.7.7.85: Diadenylate cyclase",
    "2.1.1.191: 23S rRNA (cytosine(1962)-C(5))-methyltransferase",
    "1.1.1.40: Malate dehydrogenase (oxaloacetate-decarboxylating) (NADP(+))",
    "4.1.1.49: Phosphoenolpyruvate carboxykinase (ATP)",
    "6.1.1.18: Glutamine--tRNA ligase",
    "2.6.1.83: LL-diaminopimelate aminotransferase",
    "5.1.3.13: dTDP-4-dehydrorhamnose 3,5-epimerase",
]

plot_df.previousmention = [ ( ( el ∈ union(vatanen2018_ECs, figure4_ecs) ) ? '★' : ' ' ) for el in plot_df.function_name ]

println("$(sum(plot_df.previousmention .== '★')) ($(mean(plot_df.previousmention .== '★')*100)%) of the previously-mentioned functions are present here as well!")

#####
# Building Figure
#####

function reformat_ecs(eecs::String, colset::Dict, previousmention::Char; pattern::Regex = r"(\d+\.\d+\.\d+\.?\d*):\s*(.*)")

    humann_ec_matches = eachmatch(pattern, eecs)

    ec_num = string(first(humann_ec_matches).captures[1])
    ec_fun = string(first(humann_ec_matches).captures[2])

    rich(rich("█\t"; font = :bold, color = colset[ec_num[1]]), ec_fun, rich(" [$(ec_num)]"; font = :bold, color = :black), rich(string(previousmention); font = :bold, color = :black) )
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

# Creating Master Figure 4
figureS5_master = Figure(; size = (1300, 1400))
A_subfig = GridLayout(figureS5_master[1,1], alignmode = Outside())
A_legend = GridLayout(figureS5_master[2,1], alignmode = Outside())
# BCDE_subfig = GridLayout(figureS5_master[3,1], alignmode = Outside())
# A_subfig = GridLayout(figureS5_master[1,1], alignmode = Inside())
# A_legend = GridLayout(figureS5_master[2,1], alignmode = Inside())
BCDE_subfig = GridLayout(figureS5_master[3,1], alignmode = Inside())
## Figure 2, Panel A - Functions from continuous analysis
axA = Axis(
    A_subfig[1, :],
    yticks = (1:nrow(plot_df), [ reformat_ecs(plot_df[i, "function_name"], ec_colors, plot_df[i, "previousmention"]) for i in 1:nrow(plot_df) ]),
    xticklabelsize = 18,
    yticklabelsize = 18,
    xlabel = "Effect size and direction (coefficient)",
    yaxisposition = :right,
    )

tightlimits!(axA, Top())
tightlimits!(axA, Bottom())
hidexdecorations!(axA, label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true)
hideydecorations!(axA, label = true, ticklabels = false, ticks = true, grid = true, minorgrid = true, minorticks = true)

barplot!(
    axA,
    collect(1:nrow(plot_df)),
    plot_df.score,
    color = [ ( (el > 0) ? "blue" : "red" ) for el in plot_df.score ],
    direction=:x
)

# annotations!(
#     axA,
#     plot_df.label,
#     [ Point(plot_df[i, "score"] .+ 2.0 .* sign.(plot_df[i, "score"]), i .- 0.5) for i in 1:nrow(plot_df)],
#     color = :gray10,
#     align = (:center, :bottom),
#     fontsize = 16
# )

vlines!(axA, [ 0.0 ]; color = :black)

lgd1 = Legend(
    A_legend[1,1],
    [   PolyElement(color = :blue, markersize = 20),
        PolyElement(color = :red, markersize = 20),
    ],
    [
        "Positive correlation with age",
        "Negative correlation with age",
    ],
    "Direction",
    orientation = :vertical,
    tellheight = true,
    tellwidth = true,
    labelsize = 18#,
)

lgd2 = Legend(
    A_legend[1,2],
    [
        MarkerElement(marker = '▉', color = ec_colors['1'], markersize = 14),
        MarkerElement(marker = '▉', color = ec_colors['2'], markersize = 14),
        MarkerElement(marker = '▉', color = ec_colors['3'], markersize = 14),
        MarkerElement(marker = '▉', color = ec_colors['4'], markersize = 14),
        MarkerElement(marker = '▉', color = ec_colors['5'], markersize = 14),
        MarkerElement(marker = '▉', color = ec_colors['6'], markersize = 14)
    ],
    [
        "EC 1. Oxidoreductases",
        "EC 2. Transferases",
        "EC 3. Hydrolases",
        "EC 4. Lyases",
        "EC 5. Isomerases",
        "EC 6. Ligases"
    ],
    "Main EC Group",
    nbanks = 2,
    orientation = :vertical,
    tellheight = true,
    tellwidth = true,
    labelsize = 18#,
    # margin = (900, -900, 300, -300)
)

# lgd3 = Legend(
#     A_legend[1,3],
#     [
#         MarkerElement(marker = '*', color = :black, points = Point2f[(0.2, 0.5), (0.5, 0.5), (0.8, 0.5)], markersize = 20),
#         MarkerElement(marker = '*', color = :black, points = Point2f[(0.33, 0.5), (0.66, 0.5)], markersize = 20),
#         MarkerElement(marker = '*', color = :black, markersize = 20),
#         MarkerElement(marker = '·', color = :black, markersize = 20),

#     ],
#     [
#         "p < 0.001, q < 0.2",
#         "p < 0.01, q < 0.2",
#         "p < 0.05, q < 0.2",
#         "p < 0.1, q < 0.2",
#     ],
#     "Significance",
#     orientation = :vertical,
#     tellheight = true,
#     tellwidth = true,
#     labelsize = 18#,
# )

## Scatterplots

scatterplot_functions = [
    "1.1.1.1: Alcohol dehydrogenase",
    "2.7.1.15: Ribokinase",
    "2.7.1.11: 6-phosphofructokinase",
    "2.7.7.85: Diadenylate cyclase"
]

for (i, f) in enumerate(scatterplot_functions)
    fnum = split(first(split(f, ':')), '.')
    srch = Regex("^" * join(fnum, raw"\.") * ":")

    df = comm2wide(filtered_functional_profiles[srch, :])
    df.sum = map(x -> sum(x)*1e6, eachrow(Matrix(df[:, 5:end])))
    # print(df.sum)

    select!(df, [:sample, :sum])
    allsamps = innerjoin(long_longitudinal_samples, df, on = :sample)

    glmform = @formula( sum ~ ageMonths + (1 | subject_id) )
    fitted_line = LinearMixedModel(glmform, allsamps)
    fit!(fitted_line)

    @show fitted_line
    @show coeftable(fitted_line)

    ax = Axis(
        BCDE_subfig[1,i],
        ylabel = rich("EC relative abundance\n( 10", superscript("-6"), " CPM)"),
        xlabel = "Age (months)",
        yticklabelsize = 18,
        xticklabelsize = 18
    )

    hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorgrid = true, minorticks = true)

    sc = scatter!(
        ax,
        allsamps.ageMonths,
        allsamps.sum,
        color = ( (sign(DataFrame(coeftable(fitted_line))[2,2]) == 1) ? (:blue) : (:red) )
    )

    abln = ablines!(
        ax,
        [ DataFrame(coeftable(fitted_line))[1,2] ],
        [ DataFrame(coeftable(fitted_line))[2,2] ],
        color = ( (sign(DataFrame(coeftable(fitted_line))[2,2]) == 1) ? (:blue) : (:red) )
    )

end


## Labeling and formatting
Label(A_subfig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, padding = (0, 15, 0, 0), halign = :right, alignmode = Inside())
Label(BCDE_subfig[1, 1, TopLeft()], "b", fontsize = 22, font = :bold, padding = (0, 15, 0, 0), halign = :right, alignmode = Inside())
Label(BCDE_subfig[1, 2, TopLeft()], "c", fontsize = 22, font = :bold, padding = (0, 15, 0, 0), halign = :right, alignmode = Inside())
Label(BCDE_subfig[1, 3, TopLeft()], "d", fontsize = 22, font = :bold, padding = (0, 15, 0, 0), halign = :right, alignmode = Inside())
Label(BCDE_subfig[1, 4, TopLeft()], "e", fontsize = 22, font = :bold, padding = (0, 15, 0, 0), halign = :right, alignmode = Inside())

colsize!(figureS5_master.layout, 1, Relative(1.0))
colgap!(figureS5_master.layout, 0.0)
rowsize!(figureS5_master.layout, 1, Relative(0.73))
rowsize!(figureS5_master.layout, 2, Relative(0.1))
rowsize!(figureS5_master.layout, 3, Relative(0.17))
# Export Figure 4

save(joinpath(outdir, "figures", "FigureS5.png"), figureS5_master)
save(joinpath(outdir, "figures", "FigureS5.eps"), figureS5_master)
save(joinpath(outdir, "figures", "FigureS5.svg"), figureS5_master)
save(joinpath(outdir, "figures", "FigureS5.pdf"), figureS5_master)
figureS5_master
