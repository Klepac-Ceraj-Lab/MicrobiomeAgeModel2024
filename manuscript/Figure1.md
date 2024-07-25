# Figure 1 - General Statistics, Community analysis

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
using Statistics
using CairoMakie
using Leap
using CategoricalArrays
using GLM
using StatsBase
using StableRNGs
```

### Configurable parameters
```julia
master_colors = Dict(
    "ECHO-RESONANCE" => "purple",
    "1kDLEAP-GERMINA" => "blue",
    "1kDLEAP-COMBINE" => "orange",
    "1kDLEAP-KHULA" => "red",
    "1kDLEAP-M4EFAD" => "darkgreen",
    "CMD-OTHER" => "lightblue",
    "CMD-DIABIMMUNE" => "lightblue",
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

##  Summary Tables

### Number of unique samples and subjects before prevalence filtering
```julia
replace!(combined_inputs.datasource, "CMD-DIABIMMUNE" => "CMD")
replace!(combined_inputs.datasource, "CMD-OTHER" => "CMD")

println("Number of unique stool samples: $(length(unique(combined_inputs.sample)))")
println("Number of unique subjects: $(length(unique(combined_inputs.subject_id)))")

println("Mean host age at sample collection: $(round(mean(combined_inputs.ageMonths); digits = 2))")
println("SD of host age at sample collection: $(round(Statistics.std(combined_inputs.ageMonths); digits = 2))")

countries_present = unique(combined_inputs.site)
println("$(length(countries_present)) countries represented in our dataset: $(countries_present)")
lmic_present = [ "BRA", "ZAF", "BGD", "SLV" ]
lmic_samples = subset(combined_inputs, :site => x -> x .∈ Ref(lmic_present))
println("$(nrow(lmic_samples)) samples are from LMICs, or $(nrow(lmic_samples)*100/nrow(combined_inputs)) %")
LEAP_samples = subset(combined_inputs, :datasource => x -> x .∈ Ref(["1kDLEAP-GERMINA", "1kDLEAP-COMBINE", "1kDLEAP-KHULA", "1kDLEAP-M4EFAD"]))
println("$(nrow(LEAP_samples)) samples are from 1kD-LEAP, or $(nrow(LEAP_samples)*100/nrow(combined_inputs)) %")

println("Mean host age at sample collection for LEAP samples: $(round(mean(LEAP_samples.ageMonths); digits = 2))")
println("SD of host age at sample collection for LEAP samples: $(round(Statistics.std(LEAP_samples.ageMonths); digits = 2))")

println("Proportion of LEAP samples from LMICs: $(((nrow(LEAP_samples) - sum(LEAP_samples.datasource .== "1kDLEAP-COMBINE"))/nrow(LEAP_samples)))")
```

```julia
datasource_summary_table = combine(
    groupby(combined_inputs, :datasource),
    :subject_id => (x -> length(unique(x))) => :Unique_subjects,
    :sample => (x -> length(unique(x))) => :Unique_samples,
    :ageMonths => mean => :Mean,
    :ageMonths => std => :Std,
    :ageMonths => ( x -> "$(round(mean(x); digits = 2)) ($(round(Statistics.std(x); digits = 2)))" ) => :Formatted_Mean_Std,
    )
datasource_summary_table.color = [ master_colors[el] for el in datasource_summary_table.datasource ]
@show sort!(datasource_summary_table, :Unique_subjects)
```

## Exporting Table 1
```julia
table1 = combine(
    groupby(combined_inputs, :study_name),
    :subject_id => (x -> length(unique(x))) => :Unique_subjects,
    :sample => (x -> length(unique(x))) => :Unique_samples,
    :ageMonths => mean => :Mean,
    :ageMonths => std => :Std,
    :ageMonths => ( x -> "$(round(mean(x); digits = 2)) ($(round(Statistics.std(x); digits = 2)))" ) => :Formatted_Mean_Std,
    )
@show sort!(table1, :study_name)
CSV.write("manuscript/Table1.csv", table1)
```

# Creating Master Figure 1
```julia
figure1_master = Figure(; size = (1000, 1000))

AB_Subfig = GridLayout(figure1_master[1,1], alignmode=Outside()) 
C_Subfig = GridLayout(figure1_master[1,2], alignmode=Outside())
DE_Subfig  = GridLayout(figure1_master[2,1:2], alignmode=Outside())
```

## World map of data sources
```julia
fig1_panelA = rotr90(load("manuscript/assets/Figure1_PanelA_WorldMap.PNG"))
axA = Axis(AB_Subfig[1,1], alignmode = Inside())
hidedecorations!(axA); hidespines!(axA)
image!(axA, fig1_panelA)
```

## Figure 1, Legend between A and B
```julia
Legend(
    AB_Subfig[2, 1],
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
        "1kDLEAP-Germina",
        # "1kDLEAP-KhulaMW",
        "1kDLEAP-Khula",
        "1kDLEAP-Combine",
        "1kDLEAP-M4EFaD"
    ],
    orientation = :vertical,
    nbanks = 3,
    labelsize = 14,
    tellheight = true,
    tellwidth = false
)
```

## Figure 1, Panel B - Stacked histogram of ages
```julia
# Define the number of bins and bin edges
nbins = 16
bin_edges = range(2, stop = 18, length = nbins + 1)

# Create a new column in the DataFrame to hold the bin information
combined_inputs.bin = cut(combined_inputs.ageMonths, bin_edges; labels = 1:nbins)

# Initialize a dictionary to store the binned counts
binned_counts = Dict((category, bin) => 0 for category in unique(combined_inputs.datasource), bin in 1:nbins)

# Bin the data and count occurrences
for row in eachrow(combined_inputs)
    binned_counts[(row.datasource, row.bin)] += 1
end

# Prepare the data for plotting
categories = [ "CMD", "ECHO-RESONANCE", "1kDLEAP-KHULA", "1kDLEAP-COMBINE", "1kDLEAP-GERMINA", "1kDLEAP-M4EFAD" ] #equivalent to `unique(combined_inputs.datasource)`, but hardcoded in this order for aesthetic purposes
bins = 1:nbins
bar_heights = zeros(length(categories), nbins)

for (i, category) in enumerate(categories)
    for bin in bins
        bar_heights[i, bin] = binned_counts[(category, bin)]
    end
end

# Cumulative sums for stacking
cumulative_sums = cumsum(bar_heights, dims = 1)

# Set up the figure and axis
axB = Axis(
    AB_Subfig[3, 1],
    ylabel = "Number of samples",
    xlabel = "Age in months",
    xticks = (collect(1:nbins+1) .- 0.5, string.(floor.(Int64, bin_edges))),
    alignmode = Outside()
)
axB.rightspinevisible = false
axB.topspinevisible = false
hidedecorations!(axB, label = false, ticklabels = false, ticks = false, grid = true, minorgrid = true, minorticks = true)
ylims!(axB, (0, 800))
xlims!(axB, (0.49, 16.51))
# Plot each category as a stacked bar
for bin in bins
    barplot!(axB,
        [bin for _ in eachrow(bar_heights)],
        bar_heights[:, bin],
        stack = [bin for _ in eachrow(bar_heights)],
        color = [ master_colors[el] for el in categories ],
        strokewidth = 0.0,
        width = 1.2
    )
end
```

## Figure 1, Panel B - Pie chart
```julia
fig = Figure(; size = (800, 500))
ax1 = Axis(fig[1,1], title = "Participants", autolimitaspect = 1, titlesize = 30)
ax2 = Axis(fig[1,2], title = "Samples", autolimitaspect = 1, titlesize = 30)
xlims!(ax1, (-5.5, +5.5)) 
ylims!(ax1, (-5.5, +5.5)) 
xlims!(ax2, (-5.5, +5.5)) 
ylims!(ax2, (-5.5, +5.5)) 
hidedecorations!(ax1); hidespines!(ax1)
hidedecorations!(ax2); hidespines!(ax2)

pie1 = pie!(ax1,
   datasource_summary_table.Unique_subjects, color = datasource_summary_table.color,
    radius = 4, inner_radius = 1.5, strokecolor = :white, strokewidth = 5
    )

# Calculate the angles for each segment
subject_angles = cumsum([0; datasource_summary_table.Unique_subjects]) / sum(datasource_summary_table.Unique_subjects) * 2π
# Calculate the positions for the annotations
subject_positions = (subject_angles[1:end-1] .+ subject_angles[2:end]) / 2
# Add annotations to the pie chart
for (i, pos) in enumerate(subject_positions)
    x = 4.8 * cos(pos)
    y = 4.8 * sin(pos)
    text!(ax1, x, y, text = string(datasource_summary_table.Unique_subjects[i]), align = (:center, :center), fontsize=32, color=master_colors[datasource_summary_table.datasource[i]])
end

pie2 = pie!(ax2,
   datasource_summary_table.Unique_samples, color =datasource_summary_table.color,
    radius = 4, inner_radius = 1.5, strokecolor = :white, strokewidth = 5
    )
# Calculate the angles for each segment
subject_angles = cumsum([0; datasource_summary_table.Unique_samples]) / sum(datasource_summary_table.Unique_samples) * 2π
# Calculate the positions for the annotations
subject_positions = (subject_angles[1:end-1] .+ subject_angles[2:end]) / 2
# Add annotations to the pie chart
for (i, pos) in enumerate(subject_positions)
    x = 4.8 * cos(pos)
    y = 4.8 * sin(pos)
    text!(ax2, x, y, text = string(datasource_summary_table.Unique_samples[i]), align = (:center, :center), fontsize=32, color=master_colors[datasource_summary_table.datasource[i]])
end

Legend(
    fig[2, :],
    [ MarkerElement(marker = :circle, color = master_colors[el], markersize = 14) for el in keys(master_colors) ],
    collect(keys(master_colors)),
    orientation = :vertical,
    nbanks = 4,
    tellheight = true,
    tellwidth = true
)
save(joinpath(outdir, "figures", "Figure1_PanelB_ElementII.png"), fig)
save(joinpath(outdir, "figures", "Figure1_PanelB_ElementII.eps"), fig)
save(joinpath(outdir, "figures", "Figure1_PanelB_ElementII.svg"), fig)
```
![Samples Pie Chart](../results/2024AgeModelManuscript/figures/Figure1_PanelB_ElementII.png)

## Figure 1, Panel C - Methodology Workflow
```julia
fig1_panelC = rotr90(load("manuscript/assets/Figure1_PanelC_Workflow.PNG"))
axC = Axis(C_Subfig[1,1])
hidedecorations!(axC); hidespines!(axC)
image!(axC, fig1_panelC)
```

## PERMANOVAS
```julia
spedm = Distances.pairwise(BrayCurtis(), Matrix(combined_inputs[:, 11:end-2]), dims=1)

lt4idx = combined_inputs.ageMonths .< 4.0
lt8idx = combined_inputs.ageMonths .< 8.0
commlabels = ["Taxa"]
mdlabels = [ "Age", "Country", "Data\nSource"]

pmn_all = permanovas(
    [ spedm ], 
    [
        combined_inputs.ageMonths,
        combined_inputs.site,
        combined_inputs.datasource,
    ]; commlabels, mdlabels
)

pmn_lt4  = permanovas(
    [ spedm[lt4idx, lt4idx] ],
    [
        combined_inputs.ageMonths[lt4idx],
        combined_inputs.site[lt4idx],
        combined_inputs.datasource[lt4idx],
    ]; commlabels, mdlabels
)

pmn_lt8 = permanovas(
    [ spedm[lt8idx, lt8idx] ],
    [
        combined_inputs.ageMonths[lt8idx],
        combined_inputs.site[lt8idx],
        combined_inputs.datasource[lt8idx],
    ]; commlabels, mdlabels
)

pmn_all.label .= "all samples"
pmn_lt4.label .= "< 4mo"
pmn_lt8.label .= "< 8mo"

## 2.2 Combined
fig = Figure(;size = (400,300))
ax = Axis(
    fig[1, 1];
    xticklabelsize = 16,
    yticklabelsize = 16,
    title = "PERMANOVAs",
)

plot_permanovas!(ax, vcat(pmn_all[[3,2,1],:], pmn_lt4[[3,2,1],:], pmn_lt8[[3,2,1],:]))

save(joinpath(outdir, "figures", "FigureSX_PERMANOVAs.png"), fig)
```

## NMDS (Principal Coordinate Analysis)

### Block of code to perform the Ordination
```julia
using MultivariateStats

if presence_absence
    MDS_results = StatsBase.fit(PCA, spedm; maxoutdim = 20)
    MDS_columns = DataFrame(:MDS1 => MDS_results.proj[:,1], :MDS2 => MDS_results.proj[:,2], :MDS3 => MDS_results.proj[:,3], :MDS4 => MDS_results.proj[:,4], :MDS5 => MDS_results.proj[:,5])
    MDS_variances = MDS_results.prinvars ./ sum(MDS_results.prinvars)
else
    MDS_results = StatsBase.fit(MDS, spedm; maxoutdim = 20, distances=true)
    MDS_columns = DataFrame(:MDS1 => MDS_results.U[:,1], :MDS2 => MDS_results.U[:,2], :MDS3 => MDS_results.U[:,3], :MDS4 => MDS_results.U[:,4], :MDS5 => MDS_results.U[:,5])
    MDS_variances = MDS_results.λ ./ sum(MDS_results.λ)
end

stress(MDS_results)
```

### Add NMDS to Panel D of Figure 1, Color by Site
```julia
axD = Axis(
    DE_Subfig[1,1],
    xlabel = "MDS1 ("*string(round(100*MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*MDS_variances[2]; digits = 2))*"%)",
    # title = "By data source",
    aspect = AxisAspect(1.2),
    alignmode = Outside()
)
hidedecorations!(axD, label = false, ticklabels = false, ticks = false, minorgrid = true, minorticks = true)
scD = scatter!(
    axD,
    MDS_columns[:,1],
    MDS_columns[:,2];
    color = [ (el, 0.6) for el in combined_inputs.datacolor ]
)

# save(joinpath(outdir, "figures", "Figure1_PanelDE_PCoA.png"), DE_Subfig[1,1])
# save(joinpath(outdir, "figures", "Figure1_PanelDE_PCoA.svg"), DE_Subfig[1,1])
# save(joinpath(outdir, "figures", "Figure1_PanelDE_PCoA.eps"), DE_Subfig[1,1])

```

### Add NMDS to Panel D of Figure 1, Color by Age
```julia
axE = Axis(
    DE_Subfig[1,2],
    xlabel = "MDS1 ("*string(round(100*MDS_variances[1]; digits = 2))*"%)",
    ylabel = "MDS2 ("*string(round(100*MDS_variances[2]; digits = 2))*"%)",
    # title = "By age in months",
    aspect = AxisAspect(1.2),
    alignmode = Outside()
)
hidedecorations!(axE, label = false, ticklabels = false, ticks = false, minorgrid = true, minorticks = true)
scE = scatter!(
    axE,
    MDS_columns[:,1],
    MDS_columns[:,2];
    color = combined_inputs.ageMonths,
    colormap = :viridis
)

Colorbar(DE_Subfig[1, 3], scE, tellheight = false, alignmode = Inside())

# save(joinpath(outdir, "figures", "Figure1_PanelDE_PCoA.png"), fig)
# save(joinpath(outdir, "figures", "Figure1_PanelDE_PCoA.svg"), fig)
# save(joinpath(outdir, "figures", "Figure1_PanelDE_PCoA.eps"), fig)
```

## Add labels
```julia
Label(AB_Subfig[1, 1, TopLeft()], "A", fontsize = 22,font = :bold, padding = (0, 5, 0, 0), halign = :right, alignmode = Inside())
Label(AB_Subfig[3, 1, TopLeft()], "B", fontsize = 22,font = :bold, padding = (0, 5, 0, 0), halign = :right, alignmode = Inside())
Label(C_Subfig[1, 1, TopLeft()], "C", fontsize = 22,font = :bold, padding = (0, 5, 0, 0), halign = :right, alignmode = Inside())
Label(DE_Subfig[1, 1, TopLeft()], "D", fontsize = 22, font = :bold, padding = (0, 5, 0, 0), halign = :right, alignmode = Inside())
Label(DE_Subfig[1, 2, TopLeft()], "E", fontsize = 22,font = :bold, padding = (0, 5, 0, 0), halign = :right, alignmode = Inside())
```

## Fix layout
```julia

colgap!(figure1_master.layout, 0)
rowgap!(figure1_master.layout, 0)
colgap!(AB_Subfig, 0)
rowgap!(AB_Subfig, 0)
colgap!(C_Subfig, 0)
rowgap!(C_Subfig, 0)
colgap!(DE_Subfig, 10)
rowgap!(DE_Subfig, 0)

colsize!(figure1_master.layout, 2, Relative(0.38))
rowsize!(figure1_master.layout, 2, Relative(0.4))
rowsize!(AB_Subfig, 1, Relative(0.50))
rowsize!(AB_Subfig, 2, Relative(0.15))
rowsize!(AB_Subfig, 3, Relative(0.35))
```

# Export Figure 1

```julia
save(joinpath(outdir, "figures", "Figure1.png"), figure1_master)
figure1_master
```