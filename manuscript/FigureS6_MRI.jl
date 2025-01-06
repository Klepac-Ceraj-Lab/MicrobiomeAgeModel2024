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

### Loading taxonomic profiles from all the cohorts
## This line will evoke the auxiliary notebook that contains the code to load data from all cohorts
include("/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/notebooks/allcohorts_data_loading_nofeed.jl")
combined_inputs.richness = map(x -> sum(x .> 0.0), eachrow(Matrix(combined_inputs[:, 11:ncol(combined_inputs)-1])))

### Filtering and preparing the data
##Computing samples and features before and after filters:

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

## Loading brain data

brain_df = CSV.read(inputfiles("brain_normalized.csv"), DataFrame)
fact = brain_df."White-matter" .+ brain_df."Gray-matter"
for f in names(brain_df, Not(["subject", "timepoint"]))
    # brain_df[!, f] ./= fact

    brain_df[!, f] .= log10.((brain_df[!, f] ./ fact) .+ 1e-10)
end

resmdata = @chain CSV.read("/home/guilherme/Repos/Resonance/resonance_mdata.csv", DataFrame) begin
    select!([ :subject, :timepoint, :sample])
    dropmissing()
end

brain_df = innerjoin(resmdata, brain_df, on = [ :subject => :subject, :timepoint => :timepoint ])

echo_predictions.sample = map(x -> split(x, '_')[1], echo_predictions.sample)

joined_brain_df = innerjoin(echo_predictions, brain_df, on = :sample)

#####
# Building the stats
#####

function myrf(data, target)
    # Split data into features (X) and target (y)
    y, X = unpack(data, ==(target))
    train, test = partition(collect(1:length(y)), 0.75, shuffle=true)

    model = MLJDecisionTreeInterface.RandomForestRegressor()
    mach = machine(model, X, y)
    fit!(mach, rows=train)

    # Predict on both train and test sets
    y_train_pred = MLJ.predict(mach, rows=train)
    y_test_pred = MLJ.predict(mach, rows=test)

    # Compute correlations
    train_corr = cor(y[train], y_train_pred)
    test_corr = cor(y[test], y_test_pred)

    # Evaluate RMS on test set
    rms_test = rms(y_test_pred, y[test])

    return (
        model=mach,
        train_predictions=y_train_pred,
        test_predictions=y_test_pred,
        train_correlation=train_corr,
        test_correlation=test_corr,
        test_rms=rms_test
    )
end

brain_segments = names(select(brain_df, Not([:subject, :timepoint, :sample])))

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
            n_replicas = 50,
            n_rngs = 4,
            tuning_space = (; #PRODUCTION
                maxnodes_range = [ 5 ],
                nodesize_range = [ 7 ],
                min_samples_split = [3 ],
                sampsize_range = [ 0.8 ],
                mtry_range = [ -1 ],
                ntrees_range = [ 200 ]
            )
        )   
        # @show rfres = myrf(lmod_df, :a)
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

regions_toplot = [
    "left-cuneus",
    "Brain-stem",
    "left-caudal-middle-frontal",
    "left-pericalcarine",
    "right-cuneus",
    "left-precentral",
    "right-precuneus",
    "right-accumbens-area"
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
            )
        )   
        preds = get_set_predictions(regrf, "val"; unique_col = :smp)
        plot_Df = innerjoin(preds, rf_df, on = :sample =>:smp)
    scatter!(ax5, plot_Df.test_prediction, plot_Df.a, color = "green" )

end

supp_figure6_master
save("brain_agemodel.png", supp_figure6_master)


#####
# Eczema/ECHO
#####

echo_eczema_mdata = @chain CSV.read("/home/guilherme/Repos/HMOEczema2023/data/00-demographics_RESONANCE.csv", DataFrame) begin
    
end