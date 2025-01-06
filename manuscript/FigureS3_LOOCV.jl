###########################################
# Figures S3 and S4 - External Validation #
###########################################

#####
# Figure S3 - LOOCV and Independent Test set Validation
#####

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
using DataToolkit
using MicrobiomeAgeModel2024

### Configurable parameters
master_colors = Dict(
    "ECHO" => "purple",
    "ECHO-RESONANCE" => "purple",
    "1kDLEAP-GERMINA" => "blue",
    "1kDLEAP-CORK" => "orange",
    "1kDLEAP-COMBINE" => "orange",
    "1kDLEAP-KHULA" => "red",
    "1kDLEAP-M4EFAD" => "darkgreen",
    "DIABIMMUNE" => "lightblue",
    "CMD" => "lightblue",
    "ENNIS" => "navyblue"
)

experiment_name = "2024AgeModelFinalSubmission"
outdir = joinpath(pwd(), "results", experiment_name)
figdir = joinpath(outdir, "figures")
deepdivemonodir, deepdivecolordir = ( joinpath(figdir, "species_monocolor_scatterplots"), joinpath(figdir, "species_colored_scatterplots") )
isdir(outdir) ? @warn("Directory $(outdir) already exists! This notebook will overwrite files already there.") : ( mkpath(outdir), mkpath(figdir), mkpath(deepdivemonodir), mkpath(deepdivecolordir) )
presence_absence = false # This argument will control whether the model will be based on abundances or binary presence/absence
## Loading data

#### UNCOMMENT ONLY ONE OF THE FOLLOWING 3 LINES TO PICK A SOURCE FOR THE ANALYSIS DATA
DataToolkit.loadcollection!("./Data_Local.toml")    ## Uncomment this line to use local files located on the "data" subfolder and the Local relative filesystem references
# DataToolkit.loadcollection!("./Data_AWS.toml")      ## Uncomment this line to use the datasets made available on the public AWS bucket
# DataToolkit.loadcollection!("./Data_Dryad.toml")    ## Uncomment this line to use the datasets published to Data Dryad (DOI: 10.5061/dryad.dbrv15f9z)

## Loading data
filtered_inputs = d"inputs_with_testdata"
bins = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

##  RF Model training

### Actual function for model training
##The following block of code will train the model on the combination of cohorts, performing crossvalidation and grid-search hyperparameter optimization. Training can take several hours if the hyperparameter grid is large. It is advised to train once and store the result on a `JLD2` object so it can be accessed with `JLD2.load`-like methods for downstream analysis and plotting. Hence, the block of code should be run only once per data update.
regression_Age_LeaveEchoOut = probe_regression_randomforest(
    "regression_Age_LeaveEchoOut",
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["ECHO-RESONANCE", "ENNIS"])),
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
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["CMD-OTHER", "CMD-DIABIMMUNE", "ENNIS"])),
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
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["1kDLEAP-KHULA", "ENNIS"])),
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
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["1kDLEAP-GERMINA", "ENNIS"])),
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
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["1kDLEAP-COMBINE", "ENNIS"])),
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
    subset(filtered_inputs, :datasource => x -> x .∉ Ref(["1kDLEAP-M4EFAD", "ENNIS"])),
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

regression_Age_LeaveEnnisOut = probe_regression_randomforest(
    "regression_Age_LeaveEnnisOut",
    subset(filtered_inputs, :datasource => x -> x .!== "ENNIS"),
    identity,
    collect(11:ncol(filtered_inputs)),
    :ageMonths;
    split_strat = "subject",
    ext_df = subset(filtered_inputs, :datasource => x -> x .== "ENNIS"),
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

@show sort(report_regression_merits(regression_Age_LeaveEnnisOut), :Val_RMSE_mean)
JLD2.@save joinpath(outdir, "regression_Age_LeaveEnnisOut.jld") regression_Age_LeaveEnnisOut

##After the code is run at least once, the results can then be loaded with:
JLD2.@load joinpath(outdir, "regression_Age_LeaveEchoOut.jld") regression_Age_LeaveEchoOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveCMDOut.jld") regression_Age_LeaveCMDOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveKhulaOut.jld") regression_Age_LeaveKhulaOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveGerminaOut.jld") regression_Age_LeaveGerminaOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveCombineOut.jld") regression_Age_LeaveCombineOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveM4EFADOut.jld") regression_Age_LeaveM4EFADOut
JLD2.@load joinpath(outdir, "regression_Age_LeaveEnnisOut.jld") regression_Age_LeaveEnnisOut

##Reporting the joint figures of merit (Supplementary Table 1):
joint_merits_df = vcat(
    report_regression_merits(regression_Age_LeaveEchoOut),
    report_regression_merits(regression_Age_LeaveGerminaOut),
    report_regression_merits(regression_Age_LeaveKhulaOut),
    report_regression_merits(regression_Age_LeaveCombineOut),
    report_regression_merits(regression_Age_LeaveM4EFADOut),
    report_regression_merits(regression_Age_LeaveCMDOut),
    report_regression_merits(regression_Age_LeaveEnnisOut)
)
println("RMSE of LOOCV: $(round(mean(joint_merits_df.Test_RMSE_mean); digits = 2)) +- $(round(conf_interval(joint_merits_df.Test_RMSE_mean); digits = 2))")
@show select(joint_merits_df, :model, :Val_RMSE_mean, :Val_RMSE_CI, :Test_RMSE_mean, :Test_RMSE_CI)

# Creating Supplementary Figure 2
supp_figure3_master = Figure(; size = (1000, 1000))

AB_subfig = GridLayout(supp_figure3_master[1,1]; alignmode = Inside())
C_subfig = GridLayout(supp_figure3_master[2,1]; alignmode = Inside())

## Panel A - LOCOCV Heatmap
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

lococv_mat = Matrix{Float64}(undef, 6, 6) ## Leacing Ennis out of the LOCOCV

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
axA = Axis(
    AB_subfig[1,1],
    xlabel = "Source Left Out",
    ylabel = "Source Metric",
    # title = "Leave-One-Source-Out",
    xticks = (1:6, ["ECHO", "Germina", "Combine", "Khula", "M4EFaD", "CMD"]),
    yticks = (1:6, ["ECHO", "Germina", "Combine", "Khula", "M4EFaD", "CMD"]),
    xticklabelrotation = pi/4,
    yticklabelrotation = pi/4,
    yticklabelsize=12,
    xticklabelsize=12,
    titlesize = 20,
    aspect = 1
)

hmA = heatmap!(axA, lococv_mat, colormap = cgrad(:plasma, rev = true), colorrange = (0.0, 4.0))

Colorbar(AB_subfig[2,1], hmA, label = "RMSE", vertical = false, width = 410)

## Scatterplot
axB = Axis(
    AB_subfig[1, 2];
    xlabel = "Sample collection age (months)",
    xticks = (2:1:18),
    ylabel = "Microbial age (months)",
    yticks = (2:1:18),
    aspect = AxisAspect(1.0),
    # title = "Predictions for independent Test set",
    yticklabelsize=16,
    xticklabelsize=16,
    titlesize = 20
)
hidedecorations!(axB, label = false, ticklabels = false, ticks = false, minorgrid = true, minorticks = true)
xlims!(axB, [1.99, 18.01]); ylims!(axB, [1.99, 18.01])

hp_idx = sort(report_regression_merits(regression_Age_LeaveEnnisOut), :Val_RMSE_mean)[1,1]

## Plotting the original data
@chain regression_Age_LeaveEnnisOut begin
    predictions_to_plot(filtered_inputs, age_bins, "val"; hp = hp_idx)
    scatter!(axB, _[:, :ageMonths],  _[:, :test_prediction], color = :gray70, marker = :circle)
end

## Plotting the Test data
ennis_trueages = subset(filtered_inputs, :datasource => x -> x .== "ENNIS").ageMonths
ennis_predictions = map(mean, eachrow(Matrix(regression_Age_LeaveEnnisOut.external_predictions[:, 2:end])))

sctest = scatter!(axB, ennis_trueages, ennis_predictions, color = :navyblue)

## ABLine
ablines!(axB, 0, 1; linestyle = :dash, linewidth=2, color = :gray )

## Annotations
@chain regression_Age_LeaveEnnisOut begin
    report_regression_merits()
    sort(:Val_RMSE_mean)
    annotations!(
        axB,
        [
            "RMSE = " * string(round(_[:, :Test_RMSE_mean][1]; digits = 2)) * " mo",
            "R = " * string(round(_[:, :Test_Cor_mean][1]; digits = 2))
        ],
        [Point(14.7, 3.0), Point(14.7, 2.3)];
        fontsize = 16,
        color = :navyblue,
        align = (:center, :bottom)
    )
end

## Legend
Legend(
    AB_subfig[2, 2],
    [
        MarkerElement(marker = :circle, color = :gray70, markersize = 16),
        MarkerElement(marker = :circle, color = :navyblue, markersize = 16),
    ],
    [
        "Original validation data",
        "Test data (Ennis2024, N=66)",
    ],
    orientation = :vertical,
    labelsize = 14,
    tellheight = false,
    tellwidth = false,
    margin=(0,0,0,0), #right, left, bottom, top
    alignmode = Inside()
)

## RMSE per bin per datasource
axC = Axis(
    C_subfig[1, :];
    xlabel = "Ground truth age bin (months)",
    xticks = (
        1:8,
        [
            "2-4",
            "4-6",
            "6-8",
            "8-10",
            "10-12",
            "12-14",
            "14-16",
            "16-18"
        ]),
    ylabel = "RMSE (months)",
    # title = "RMSE distribution along dynamic range and different cohorts",
    yticklabelsize=16,
    xticklabelsize=16,
    titlesize = 20
)

hidedecorations!(axC, label = false, ticklabels = false, ticks = false, minorgrid = true, minorticks = true)

function new_predictions_to_plot( ## changed on Leap source
    rmod::ProbeData,
    original_inputs::DataFrame,
    bins::Vector{T} where T <: Real,
    sampleset="val";
    metadata_cols = [ "study_name", "datagroup", "site", "datacolor", "datasource", "visit", "westernized_cat", "subject_id", "ageMonths", "sample"],
    hp = 1)

    combined_ploterror_df = @chain rmod begin
        get_set_predictions(sampleset; hp = hp)
        leftjoin(select(original_inputs, metadata_cols), _, on = :sample )
        dropmissing()
        insertcols!(_, 1, :l1_error => ( _.test_prediction .- _.ageMonths  ))
        insertcols!(_, 1, :absolute_error => abs.( _.test_prediction .- _.ageMonths ))
        insertcols!(_, 1, :bin_idx => map(x -> searchsortedlast(bins, x), _.ageMonths))
        insertcols!(_, 1, :model_grp => 2)
    end
    return combined_ploterror_df
end

original_predictions = new_predictions_to_plot(regression_Age_LeaveEnnisOut, regression_Age_LeaveEnnisOut.original_df, age_bins)

test_predictions = new_predictions_to_plot(regression_Age_LeaveEnnisOut, insertcols(ennis_pre_data, 1, :datagroup => "ENNIS"), age_bins, "test")
test_predictions.datacolor .= "navyblue"

boxplot_df = vcat(original_predictions, test_predictions)

dodge_map = Dict(
    "1kDLEAP-M4EFAD" => 1,
    "1kDLEAP-KHULA" => 2,
    "ECHO-RESONANCE" => 3,
    "CMD-OTHER" => 4,
    "CMD-DIABIMMUNE" => 4,
    "1kDLEAP-COMBINE" => 5,
    "1kDLEAP-GERMINA" => 6,
    "ENNIS" => 7
)

boxplot_df.dodge = [ dodge_map[el] for el in boxplot_df.datasource ]

bpC = boxplot!(
    axC,
    boxplot_df.bin_idx,
    boxplot_df.absolute_error;
    dodge = boxplot_df.dodge,
    color = boxplot_df.datacolor
    )

Legend(
    C_subfig[2, 1],
    [
        MarkerElement(marker = :circle, color = :lightblue, markersize = 14),
        MarkerElement(marker = :circle, color = :purple, markersize = 14),
        MarkerElement(marker = :circle, color = :blue, markersize = 14),
        MarkerElement(marker = :circle, color = :red, markersize = 14),
        MarkerElement(marker = :circle, color = :orange, markersize = 14),
        MarkerElement(marker = :circle, color = :darkgreen, markersize = 14),
        MarkerElement(marker = :circle, color = :navyblue, markersize = 14),
    ],
    [
        "CMD",
        "ECHO-Resonance",
        "1kDLEAP-Germina",
        "1kDLEAP-Khula",
        "1kDLEAP-Combine",
        "1kDLEAP-M4EFaD",
        "Ennis2024 (Test)",
    ],
    orientation = :vertical,
    nbanks = 4,
    labelsize = 14,
    tellheight = true,
    tellwidth = false,
    margin=(0,0,0,0), #right, left, bottom, top
    alignmode = Inside()
)
## Labeling, aligning, tweaking...
Label(AB_subfig[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, padding = (0, 10, 00, 0), halign = :right, alignmode = Inside())
Label(AB_subfig[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, padding = (0, 25, 00, 0), halign = :right, alignmode = Inside())
Label(C_subfig[1, 1, TopLeft()], "c", fontsize = 22, font = :bold, padding = (0, 10, 00, 0), halign = :right, alignmode = Inside())

rowgap!(AB_subfig, 0.0)
colgap!(AB_subfig, 0.0)
rowsize!(supp_figure3_master.layout, 1, Relative(0.6))
rowsize!(supp_figure3_master.layout, 2, Relative(0.4))

## Exporting Supplementary Figure S3
save(joinpath(outdir, "figures", "FigureS3.png"), supp_figure3_master)
save(joinpath(outdir, "figures", "FigureS3.eps"), supp_figure3_master)
save(joinpath(outdir, "figures", "FigureS3.svg"), supp_figure3_master)
save(joinpath(outdir, "figures", "FigureS3.pdf"), supp_figure3_master)
supp_figure3_master

#####
# Figure S4 - Comparing importance for different data sources
#####

collect_importances = Vector{DataFrame}(undef, length(results_vector))

for (i, res) in enumerate(results_vector)
    hp_idx = sort(report_regression_merits(regression_Age_LeaveEnnisOut), :Val_RMSE_mean)[1,1]
    imp_df = hpimportances(res, hp_idx)
    rename!(imp_df, :weightedImportance => res.name)
   collect_importances[i] = imp_df
end

concat_importances = reduce( (x,y) -> innerjoin(x,y, on = :variable), collect_importances)
concat_importances.meanAcross = map(mean, eachrow(concat_importances[:, 2:end]))
# concat_importances.meanAcrossOrder = sortperm(concat_importances.meanAcross; rev = true) ## THIS DESTROYS MICRO ORDER, DO NOT RUN
# sort!(concat_importances, :meanAcrossOrder;) ## THIS DESTROYS MICRO ORDER, DO NOT RUN
sort!(concat_importances, :meanAcross; rev=true)
concat_importances.meanAcrossOrder .= 1:nrow(concat_importances)
nfeat_toplot = 30
concat_importances = concat_importances[1:nfeat_toplot, :]

plot_impAcross_df = vcat(
    DataFrame(
        :testCohort => "ECHO-RESONANCE",
        :variable => concat_importances.variable,
        :importance => concat_importances.regression_Age_LeaveEchoOut,
        :meanAcrossOrder => concat_importances.meanAcrossOrder,
        :dodge => 1,
        :datacolor => master_colors["ECHO-RESONANCE"]
    ),
    DataFrame(
        :testCohort => "1kDLEAP-GERMINA",
        :variable => concat_importances.variable,
        :importance => concat_importances.regression_Age_LeaveGerminaOut,
        :meanAcrossOrder => concat_importances.meanAcrossOrder,
        :dodge => 2, :datacolor => master_colors["1kDLEAP-GERMINA"]
    ),
    DataFrame(
        :testCohort => "1kDLEAP-COMBINE",
        :variable => concat_importances.variable,
        :importance => concat_importances.regression_Age_LeaveCombineOut,
        :meanAcrossOrder => concat_importances.meanAcrossOrder,
        :dodge => 3,
        :datacolor => master_colors["1kDLEAP-COMBINE"]
    ),
    DataFrame(
        :testCohort => "1kDLEAP-KHULA",
        :variable => concat_importances.variable,
        :importance => concat_importances.regression_Age_LeaveKhulaOut,
        :meanAcrossOrder => concat_importances.meanAcrossOrder,
        :dodge => 4,
        :datacolor => master_colors["1kDLEAP-KHULA"]
    ),
    DataFrame(
        :testCohort => "1kDLEAP-M4EFAD",
        :variable => concat_importances.variable,
        :importance => concat_importances.regression_Age_LeaveM4EFADOut,
        :meanAcrossOrder => concat_importances.meanAcrossOrder,
        :dodge => 5,
        :datacolor => master_colors["1kDLEAP-M4EFAD"]
    ),
    DataFrame(
        :testCohort => "CMD",
        :variable => concat_importances.variable,
        :importance => concat_importances.regression_Age_LeaveCMDOut,
        :meanAcrossOrder => concat_importances.meanAcrossOrder,
        :dodge => 6,
        :datacolor => master_colors["CMD"]
    ),
)

figureS4_master = Figure(size = (1000, 1400))

ax = Axis(
    figureS4_master[1, 1:2];
    xlabel = "Proportional importance (absolute value)",
    yticks = (reverse(collect(1:nfeat_toplot)), [ replace(el, "_" => " ") for el in concat_importances.variable[1:nfeat_toplot] ]),
    # ylabel = "Predictor",
    yticklabelsize=20,
    yticklabelfont="TeX Gyre Heros Makie Italic",
    alignmode = Outside()
)

tightlimits!(ax, Top())
tightlimits!(ax, Bottom())

hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

barplot!(
    ax,
    ( 31 .+ plot_impAcross_df.meanAcrossOrder .* -1),
    plot_impAcross_df.importance,
    dodge = plot_impAcross_df.dodge,
    color = plot_impAcross_df.datacolor,
    direction=:x
)

cmd_markerplot_df = @chain plot_impAcross_df begin
    subset(:testCohort => x -> x .== "CMD")
    subset(:variable => x -> x .∈ Ref([
            "Dorea_longicatena",
            "Dorea_formicigenerans",
            "Prevotella_copri",
            "Blautia_obeum",
            "Agathobaculum_butyriciproducens"
        ])
    )
end 
khula_markerplot_df = @chain plot_impAcross_df begin
    subset(:testCohort => x -> x .== "1kDLEAP-KHULA")
    subset(:variable => x -> x .∈ Ref([
            "Eubacterium_eligens",
            "Ruminococcus_bromii"
        ])
    )
end 
markerplot_df = vcat(cmd_markerplot_df, khula_markerplot_df)

scatter!(
    ax,
    markerplot_df.importance .+ 0.005,
    ( 31 .+ markerplot_df.meanAcrossOrder .* -1),
    color = markerplot_df.datacolor,
    marker = :star5,
    markersize = 30
)

Legend(
    figureS4_master[2, 1],
    [
        MarkerElement(marker = :circle, color = :lightblue, markersize = 14),
        MarkerElement(marker = :circle, color = :purple, markersize = 14),
        MarkerElement(marker = :circle, color = :blue, markersize = 14),
        MarkerElement(marker = :circle, color = :red, markersize = 14),
        MarkerElement(marker = :circle, color = :orange, markersize = 14),
        MarkerElement(marker = :circle, color = :darkgreen, markersize = 14)
    ],
    [
        "CMD",
        "ECHO-Resonance",
        "1kDLEAP-Germina",
        "1kDLEAP-Khula",
        "1kDLEAP-Combine",
        "1kDLEAP-M4EFaD"
    ],
    rich("Data source ", rich("left out", font = :bold)),
    orientation = :vertical,
    nbanks = 3,
    labelsize = 16,
    tellheight = true,
    tellwidth = false,
    margin=(100,0,-20,20), #right, left, bottom, top
    alignmode = Inside()
)

Legend(
    figureS4_master[2, 2],
    [
        MarkerElement(marker = :star5, color = :lightblue, markersize = 14),
        MarkerElement(marker = :star5, color = :red, markersize = 14)
    ],
    [
        "Non-Western dominant",
        "Western dominant",
        # "Predominant/Exclusive of non-westernized samples",
        # "Predominant/Exclusive of westernized samples",
    ],
    "Species highlight",
    orientation = :vertical,
    nbanks = 1,
    labelsize = 16,
    tellheight = true,
    tellwidth = false,
    margin=(50,0,-20,20), #right, left, bottom, top
    alignmode = Inside()
)

rowsize!(figureS4_master.layout, 1, Relative(0.9))
rowsize!(figureS4_master.layout, 2, Relative(0.1))

figureS4_master

## Exporting Supplementary Figure S3
save(joinpath(outdir, "figures", "FigureS4.png"), figureS4_master)
save(joinpath(outdir, "figures", "FigureS4.eps"), figureS4_master)
save(joinpath(outdir, "figures", "FigureS4.svg"), figureS4_master)
save(joinpath(outdir, "figures", "FigureS4.pdf"), figureS4_master)

#####
# PPT version
#####

figureS4_PPT_master = Figure(size = (1200, 600))

ax = Axis(
    figureS4_PPT_master[1, 1];
    ylabel = "Proportional importance (absolute value)",
    xticks = (collect(1:nfeat_toplot), [ replace(el, "_" => " ") for el in concat_importances.variable[1:nfeat_toplot] ]),
    # ylabel = "Predictor",
    xticklabelsize=12,
    xticklabelrotation=-pi/4,
    xticklabelfont="TeX Gyre Heros Makie Italic"
)

tightlimits!(ax, Top())
tightlimits!(ax, Bottom())

hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorgrid = false, minorticks = false)

barplot!(
    ax,
    plot_impAcross_df.meanAcrossOrder,
    plot_impAcross_df.importance,
    dodge = plot_impAcross_df.dodge,
    color = plot_impAcross_df.datacolor
)

Legend(
    figureS4_PPT_master[2, 1],
    [
        MarkerElement(marker = :circle, color = :lightblue, markersize = 14),
        MarkerElement(marker = :circle, color = :purple, markersize = 14),
        MarkerElement(marker = :circle, color = :blue, markersize = 14),
        MarkerElement(marker = :circle, color = :red, markersize = 14),
        MarkerElement(marker = :circle, color = :orange, markersize = 14),
        MarkerElement(marker = :circle, color = :darkgreen, markersize = 14)
    ],
    [
        "CMD",
        "ECHO-Resonance",
        "1kDLEAP-Germina",
        "1kDLEAP-Khula",
        "1kDLEAP-Combine",
        "1kDLEAP-M4EFaD"
    ],
    rich("Data source ", rich("left out", font = :bold)),
    orientation = :vertical,
    nbanks = 3,
    labelsize = 16,
    tellheight = true,
    tellwidth = false,
    # margin=(0,0,-20,20), #right, left, bottom, top
    alignmode = Inside()
)

colsize!(figureS4_PPT_master.layout, 1, Relative(0.85))
rowsize!(figureS4_PPT_master.layout, 1, Relative(0.75))
rowsize!(figureS4_PPT_master.layout, 2, Relative(0.25))

figureS4_PPT_master

## Exporting Supplementary Figure S3
save(joinpath(outdir, "figures", "FigureS4_PPT.png"), figureS4_PPT_master)
save(joinpath(outdir, "figures", "FigureS4_PPT.eps"), figureS4_PPT_master)
save(joinpath(outdir, "figures", "FigureS4_PPT.svg"), figureS4_PPT_master)
save(joinpath(outdir, "figures", "FigureS4_PPT.pdf"), figureS4_PPT_master)
