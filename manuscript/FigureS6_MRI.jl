###########################################
# Figures S6 and S7 - Outcomes #
###########################################

#####
# Figure S6 - MRI and age model
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
using MLJDecisionTreeInterface
using DataToolkit

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

experiment_name = "2024AgeModelRevisions"
outdir = joinpath(pwd(), "results", experiment_name)
figdir = joinpath(outdir, "figures")
deepdivemonodir, deepdivecolordir = ( joinpath(figdir, "species_monocolor_scatterplots"), joinpath(figdir, "species_colored_scatterplots") )
isdir(outdir) ? @warn("Directory $(outdir) already exists! This notebook will overwrite files already there.") : ( mkpath(outdir), mkpath(figdir), mkpath(deepdivemonodir), mkpath(deepdivecolordir) )
presence_absence = false # This argument will control whether the model will be based on abundances or binary presence/absence
## Loading data

DataToolkit.loadcollection!("./Data_Local.toml")    ## Uncomment this line to use local files located on the "data" subfolder and the Local relative filesystem references
filtered_inputs = d"filtered_taxonomic_inputs"

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

# Getting the ECHO predictions

hp_idx = sort(report_regression_merits(regression_Age_LeaveEnnisOut), :Val_RMSE_mean)[1,1]
age_bins = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]

echo_predictions =  @chain regression_Age_LeaveEnnisOut begin
    predictions_to_plot(filtered_inputs, age_bins, "val"; hp = hp_idx)
    subset!(:datasource => x -> x .== "ECHO-RESONANCE")
end

# ## Loading brain data //Version 0.1 - old brain data

# # brain_df = CSV.read(inputfiles("brain_normalized.csv"), DataFrame)
# fact = brain_df."White-matter" .+ brain_df."Gray-matter"
# for f in names(brain_df, Not(["subject", "timepoint"]))
#     # brain_df[!, f] ./= fact

#     brain_df[!, f] .= log10.((brain_df[!, f] ./ fact) .+ 1e-10)
# end

# resmdata = @chain CSV.read("/home/guilherme/Repos/Resonance/resonance_mdata.csv", DataFrame) begin
#     select!([ :subject, :timepoint, :sample])
#     dropmissing()
# end

# brain_df = innerjoin(resmdata, brain_df, on = [ :subject => :subject, :timepoint => :timepoint ])

# echo_predictions.sample = map(x -> split(x, '_')[1], echo_predictions.sample)

# joined_brain_df = innerjoin(echo_predictions, brain_df, on = :sample)



brain_df = @chain CSV.read("data/2024-12_EchoNewBrainVolumes.csv", DataFrame) begin
    select!(Not(:subject))
    rename!( :StudyID => :subject, :Timepoint => :timepoint)
end



resmdata = @chain CSV.read("/home/guilherme/Repos/Resonance/resonance_mdata.csv", DataFrame) begin
    select!([ :subject, :timepoint, :sample])
    dropmissing()
end

brain_df = innerjoin(resmdata, brain_df, on = [ :subject => :subject, :timepoint => :timepoint ])

echo_predictions.sample = map(x -> split(x, '_')[1], echo_predictions.sample)
joined_brain_df = innerjoin(echo_predictions, brain_df, on = :sample)

samps_to_remove = [ "SEQ02403" ]
subset!(joined_brain_df, :sample => x -> x .∉ Ref(samps_to_remove))

CSV.write("new_brain_inputs.csv", joined_brain_df)

#####
# Building the stats
#####

# function myrf(data, target)
#     # Split data into features (X) and target (y)
#     y, X = unpack(data, ==(target))
#     train, test = partition(collect(1:length(y)), 0.75, shuffle=true)

#     model = MLJDecisionTreeInterface.RandomForestRegressor()
#     mach = machine(model, X, y)
#     fit!(mach, rows=train)

#     # Predict on both train and test sets
#     y_train_pred = MLJ.predict(mach, rows=train)
#     y_test_pred = MLJ.predict(mach, rows=test)

#     # Compute correlations
#     train_corr = cor(y[train], y_train_pred)
#     test_corr = cor(y[test], y_test_pred)

#     # Evaluate RMS on test set
#     rms_test = rms(y_test_pred, y[test])

#     return (
#         model=mach,
#         train_predictions=y_train_pred,
#         test_predictions=y_test_pred,
#         train_correlation=train_corr,
#         test_correlation=test_corr,
#         test_rms=rms_test
#     )
# end

# brain_segments = names(select(brain_df, Not([:subject, :timepoint, :sample])))
brain_segments = [
    "left_cerebral_white_matter",
    "left_cerebral_cortex",
    "left_lateral_ventricle",
    "left_inferior_lateral_ventricle",
    "left_cerebellum_white_matter",
    "left_cerebellum_cortex",
    "left_thalamus",
    "left_caudate",
    "left_putamen",
    "left_pallidum",
    "3rd_ventricle",
    "4th_ventricle",
    "brain-stem",
    "left_hippocampus",
    "left_amygdala",
    "csf",
    "left_accumbens_area",
    "left_ventral_DC",
    "right_cerebral_white_matter",
    "right_cerebral_cortex",
    "right_lateral_ventricle",
    "right_inferior_lateral_ventricle",
    "right_cerebellum_white_matter",
    "right_cerebellum_cortex",
    "right_thalamus",
    "right_caudate",
    "right_putamen",
    "right_pallidum",
    "right_hippocampus",
    "right_amygdala",
    "right_accumbens_area",
    "right_ventral_DC",
    "ctx-lh-bankssts",
    "ctx-lh-caudalanteriorcingulate",
    "ctx-lh-caudalmiddlefrontal",
    "ctx-lh-cuneus",
    "ctx-lh-entorhinal",
    "ctx-lh-fusiform",
    "ctx-lh-inferiorparietal",
    "ctx-lh-inferiortemporal",
    "ctx-lh-isthmuscingulate",
    "ctx-lh-lateraloccipital",
    "ctx-lh-lateralorbitofrontal",
    "ctx-lh-lingual",
    "ctx-lh-medialorbitofrontal",
    "ctx-lh-middletemporal",
    "ctx-lh-parahippocampal",
    "ctx-lh-paracentral",
    "ctx-lh-parsopercularis",
    "ctx-lh-parsorbitalis",
    "ctx-lh-parstriangularis",
    "ctx-lh-pericalcarine",
    "ctx-lh-postcentral",
    "ctx-lh-posteriorcingulate",
    "ctx-lh-precentral",
    "ctx-lh-precuneus",
    "ctx-lh-rostralanteriorcingulate",
    "ctx-lh-rostralmiddlefrontal",
    "ctx-lh-superiorfrontal",
    "ctx-lh-superiorparietal",
    "ctx-lh-superiortemporal",
    "ctx-lh-supramarginal",
    "ctx-lh-frontalpole",
    "ctx-lh-temporalpole",
    "ctx-lh-transversetemporal",
    "ctx-lh-insula",
    "ctx-rh-bankssts",
    "ctx-rh-caudalanteriorcingulate",
    "ctx-rh-caudalmiddlefrontal",
    "ctx-rh-cuneus",
    "ctx-rh-entorhinal",
    "ctx-rh-fusiform",
    "ctx-rh-inferiorparietal",
    "ctx-rh-inferiortemporal",
    "ctx-rh-isthmuscingulate",
    "ctx-rh-lateraloccipital",
    "ctx-rh-lateralorbitofrontal",
    "ctx-rh-lingual",
    "ctx-rh-medialorbitofrontal",
    "ctx-rh-middletemporal",
    "ctx-rh-parahippocampal",
    "ctx-rh-paracentral",
    "ctx-rh-parsopercularis",
    "ctx-rh-parsorbitalis",
    "ctx-rh-parstriangularis",
    "ctx-rh-pericalcarine",
    "ctx-rh-postcentral",
    "ctx-rh-posteriorcingulate",
    "ctx-rh-precentral",
    "ctx-rh-precuneus",
    "ctx-rh-rostralanteriorcingulate",
    "ctx-rh-rostralmiddlefrontal",
    "ctx-rh-superiorfrontal",
    "ctx-rh-superiorparietal",
    "ctx-rh-superiortemporal",
    "ctx-rh-supramarginal",
    "ctx-rh-frontalpole",
    "ctx-rh-temporalpole",
    "ctx-rh-transversetemporal",
    "ctx-rh-insula"
]

variables = String[]
age_cors = Float64[]
age_pval = Float64[]
pred_cors = Float64[]
pred_pval = Float64[]
error_cors = Float64[]
error_pval = Float64[]
jointmodel_cors = Float64[]
rfmodel_traincors = Float64[]
rfmodel_testcors = Float64[]

for this_segment in brain_segments

    push!(variables, this_segment)

    joined_brain_df[:, this_segment] = joined_brain_df[:, this_segment] ./ joined_brain_df[:,"total_intracranial"]
    joined_brain_df[:, this_segment] = -log10.(joined_brain_df[:, this_segment] .+ 1e-10)

    @show sort(select(joined_brain_df, ["sample", this_segment]), Symbol(this_segment); rev=true)

    push!(age_cors, cor(joined_brain_df.ageMonths, joined_brain_df[:, this_segment]))
    agelmod = lm(@formula(a ~ b), DataFrame(:a => joined_brain_df[:, this_segment], :b => joined_brain_df.ageMonths))
    push!(age_pval, coeftable(agelmod).cols[4][2])

    push!(pred_cors, cor(joined_brain_df.test_prediction, joined_brain_df[:, this_segment]))
    predlmod = lm(@formula(a ~ b), DataFrame(:a => joined_brain_df[:, this_segment], :b => joined_brain_df.test_prediction))
    push!(pred_pval, coeftable(predlmod).cols[4][2])

    push!(error_cors, cor(joined_brain_df.l1_error, joined_brain_df[:, this_segment]))
    errorlmod = lm(@formula(a ~ b), DataFrame(:a => joined_brain_df[:, this_segment], :b => joined_brain_df.l1_error))
    push!(error_pval, coeftable(errorlmod).cols[4][2])

    lmod_df = DataFrame(:a => joined_brain_df[:, this_segment], :b => joined_brain_df.ageMonths , :c => joined_brain_df.l1_error)
    linpredmod = lm(@formula(a ~ b + c), lmod_df)
    push!(jointmodel_cors, cor(GLM.predict(linpredmod, lmod_df), joined_brain_df[:, this_segment]))

    rf_df = DataFrame(:smp => joined_brain_df.sample, :a => joined_brain_df[:, this_segment], :b => joined_brain_df.ageMonths , :c => joined_brain_df.l1_error)

    try
        regrf = probe_regression_randomforest(
            this_segment,
            rf_df,
            identity,
            [3, 4],
            :a;
            split_strat = nothing,
            unique_col = :smp,
            n_folds = 4,
            n_replicas = 10,
            n_rngs = 4,
            tuning_space = (; #PRODUCTION
                maxnodes_range = [ 3, 5 ],
                nodesize_range = [ 7, 9 ],
                min_samples_split = [ 2, 3 ],
                sampsize_range = [ 0.8 ],
                mtry_range = [ -1, 0 ],
                ntrees_range = [ 128 ]
            ),
            verbose = false
        )   
        # @show rfres = myrf(lmod_df, :a)
        @show report_regression_merits(regrf)
        push!(rfmodel_traincors, report_regression_merits(regrf).Train_Cor_mean[1])
        push!(rfmodel_testcors, report_regression_merits(regrf).Val_Cor_mean[1])
    catch
        push!(rfmodel_traincors, -1.0)
        push!(rfmodel_testcors, -1.0)
    end

end

correlations_df = DataFrame(
    "variable" => variables,
    "age_cor" => age_cors,
    "age_cor_sq" => age_cors .^2,
    "age_pval" => age_pval,
    "pred_cor" => pred_cors,
    "pred_pval" => pred_pval,
    "l1_cor" => error_cors,
    "l1_pval" => error_pval,
    "jointmodel_cor" => jointmodel_cors,
    "rf_traincor" => rfmodel_traincors,
    "rf_testcor" => rfmodel_testcors
)

correlations_df.lindiff = abs.(correlations_df.jointmodel_cor) .- abs.(correlations_df.age_cor)
correlations_df.rfdiff = abs.(correlations_df.rf_testcor) .- abs.(correlations_df.age_cor)
correlations_df.modiff = abs.(correlations_df.rf_testcor) .- abs.(correlations_df.jointmodel_cor)
## Build the predictive models

sort(correlations_df, :modiff)

## Group 1 - Very good by itself, not much added by mbiome age
# ctx-lh-fusiform
# ctx-rh-entorhinal
# ctx-lh-transversetemporal
# ctx-lh-frontalpole
# left_thalamus
# right_thalamus

## Group 2 - Not really good by itself, mbiome age perhaps improves?
# ctx-lh-parsorbitalis

## Group 4 - Behaves better with Age GAP
# ctx-lh-bankssts
# ctx-lh-precuneus
# ctx-rh-paracentral

# Group 5 - really shines on RF
# ctx-lh-parsorbitalis
# ctx-rh-superiorfrontal
# right_putamen
# ctx-rh-caudalanteriorcingulate
# ctx-rh-precentral
# ctx-rh-superiorparietal

# regions_toplot = [
#     "left-cuneus",
#     "Brain-stem",
#     "left-caudal-middle-frontal",
#     "left-pericalcarine",
#     "right-cuneus",
#     "left-precentral",
#     "right-precuneus",
#     "right-accumbens-area"
# ]

regions_toplot = [
    "ctx-lh-fusiform",
    "left_thalamus",
    "right_cerebellum_cortex",
    "ctx-rh-precentral",
    "ctx-rh-paracentral",
    "ctx-lh-parsorbitalis",
    "right_putamen",
]

idxes_toplot = [ 1,2,3,4,5,6
    # (1,1),
    # (1,2),
    # (2,1),
    # (2,2),
    # (3,1),
    # (3,2),
    # (4,1),
    # (4,2)
]

# Creating Supplementary Figure 6
supp_figure6_master = Figure(; size = (1800, 1000))

for (i, j) in enumerate(regions_toplot)
    ax1 = Axis(supp_figure6_master[1,i], title = j, xlabel = "Age in Months", ylabel = "log(rel.vol.)")
    scatter!(ax1, joined_brain_df.ageMonths, joined_brain_df[:, j], color = "red" )

    ax2 = Axis(supp_figure6_master[2,i], xlabel = "Predicted Age", ylabel = "log(rel.vol.)")
    scatter!(ax2, joined_brain_df.test_prediction, joined_brain_df[:, j], color = "orange" )

    ax3 = Axis(supp_figure6_master[3,i], xlabel = "Age gap (residual)", ylabel = "log(rel.vol.)")
    scatter!(ax3, joined_brain_df.l1_error, joined_brain_df[:, j], color = "purple" )

    ax4 = Axis(supp_figure6_master[4,i], xlabel = "age + gap (linear)", ylabel = "log(rel.vol.)")
    lmod_df = DataFrame(:a => joined_brain_df[:, j], :b => joined_brain_df.ageMonths , :c => joined_brain_df.l1_error)
    linpredmod = lm(@formula(a ~ b + c), lmod_df)
    scatter!(ax4, GLM.predict(linpredmod, lmod_df), joined_brain_df[:, j], color = "blue" )

    ax5 = Axis(supp_figure6_master[5,i], xlabel = "age + gap (nonlinear)", ylabel = "log(rel.vol.)")
    rf_df = DataFrame(:smp => joined_brain_df.sample, :a => joined_brain_df[:, j], :b => joined_brain_df.ageMonths , :c => joined_brain_df.l1_error)
    show(rf_df, allrows=true)
    regrf = probe_regression_randomforest(
            j,
            rf_df,
            identity,
            [3, 4],
            :a;
            split_strat = nothing,
            unique_col = :smp,
            n_folds = 4,
            n_replicas = 50,
            n_rngs = 4,
            tuning_space = (; #PRODUCTION
                maxnodes_range = [ 5 ],
                nodesize_range = [ 7 ],
                min_samples_split = [3 ],
                sampsize_range = [ 0.8 ],
                mtry_range = [ -1 ],
                ntrees_range = [ 200 ]
            ), verbose=false
        )   
        preds = get_set_predictions(regrf, "val"; unique_col = :smp)
        plot_Df = innerjoin(preds, rf_df, on = :sample =>:smp)
    scatter!(ax5, plot_Df.test_prediction, plot_Df.a, color = "green" )

end

supp_figure6_master
save("brain_agemodel_newdata.png", supp_figure6_master)


#####
# Eczema/ECHO
#####

echo_eczema_mdata = @chain CSV.read("/home/guilherme/Repos/HMOEczema2023/data/00-demographics_RESONANCE.csv", DataFrame) begin
    
end