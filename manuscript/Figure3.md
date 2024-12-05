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
using CSV
using DataToolkit
using MicrobiomeAgeModel2024
```

### Configurable parameters and notebook set-up
```julia
outdir, figdir, deepdivemonodir, deepdivecolordir = setup_outdir(; experiment_name = "2024AgeModelRevisions")
presence_absence = false # This argument controls whether the analysis will be based on continous relative abundances or binary presence/absence of species.
```
#### UNCOMMENT ONLY ONE OF THE FOLLOWING 3 LINES TO PICK A SOURCE FOR THE ANALYSIS DATA
```julia
# DataToolkit.loadcollection!("./Data_Local.toml")    ## Uncomment this line to use local files located on the "data" subfolder and the Local relative filesystem references
DataToolkit.loadcollection!("./Data_AWS.toml")      ## Uncomment this line to use the datasets made available on the public AWS bucket
# DataToolkit.loadcollection!("./Data_Dryad.toml")    ## Uncomment this line to use the datasets published to Data Dryad (DOI: 10.5061/dryad.dbrv15f9z)
```

## Loading data
```julia
regression_Age_FullCV = d"cv_results"["regression_Age_FullCV"]
taxonomic_profiles = regression_Age_FullCV.original_df

bins = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
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

## To plot the dendrogram
# save(joinpath(outdir, "figures", "prevalences_dendrogram.png"), StatsPlots.plot(hcl_taxa_prevalences))

correlations_df = [
    (; 
        species = important_bugs[i],
        cor1 = cor(cmd_prevalence_matrix[i, :], echo_prevalence_matrix[i, :] ),
        cor2 = cor(cmd_prevalence_matrix[i, :], khula_prevalence_matrix[i, :] ),
        cor3 = cor(khula_prevalence_matrix[i, :], echo_prevalence_matrix[i, :] ),
    ) for i in eachindex(important_bugs)
]
correlations_df = DataFrame(correlations_df)
correlations_df.geomeans = ([ maximum([0.0, el]) for el in correlations_df.cor1 .* correlations_df.cor2 .* correlations_df.cor3 ]) .^ 0.33
# sort!(correlations_df, :geomeans)
```

# Creating Master Figure 3
```julia
figure3_master = Figure(; size = (1600, 1000))
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
    title = "Baltic samples",
    xticks = (eachindex(bins), ["2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14"]),
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
    title = "North American samples",
    xticks = (eachindex(bins), ["2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14"]),
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
    title = "South African samples",
    xticks = (eachindex(bins), ["2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14"]),
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

Colorbar(figure3_master[1,5], hmA, label = "Prevalence",  ticks = 0.0:0.2:1.0, ticklabelsize = 24, labelsize = 24)

## Adding the HCA pdendrogram1

# Node struct
struct DNode{N}
    idx::Int
    position::Point{N, Float32}
    children::Union{Tuple{Int,Int}, Nothing}
end

function DNode(idx::Int, point::Point{N}, children::Union{Tuple{Int,Int}, Nothing}) where N
    return DNode{N}(idx, point, children)
end

"""
    dendrogram(x, y; kwargs...)

Draw a [dendrogram](https://en.wikipedia.org/wiki/Dendrogram),
with leaf nodes specified by `x` and `y` coordinates,
and parent nodes identified by `merges`.
# Arguments
- `x`: x positions of leaf nodes
- `y`: y positions of leaf nodes (default = 0)
# Keywords
- `merges`: specifies connections between nodes (see below)
- `treestyle`: one of `:tree`, `:box`.  Overload `dendrogram_connectors(::Val{:mystyle}, parent, child1, child2)` to define a new style.
"""
@recipe(Dendrogram, nodes) do scene
    Theme(
        weights = Makie.automatic,
        branch_shape = :box,
        linewidth = Makie.inherit(scene, :linewidth, 1.0),
        color = Makie.inherit(scene, :color, :black),
        colormap = Makie.inherit(scene, :colormap, :tableau_10),
        colorrange = Makie.automatic,
        orientation = :vertical,
        groups = nothing,
        cycle = [:color => :patchcolor],
        inspectable = Makie.inherit(scene, :inspectable, false),
        xautolimits = Makie.inherit(scene, :xautolimits, true),
        yautolimits = Makie.inherit(scene, :yautolimits, true),
    )
end

function recursive_dendrogram_points(node, nodes, ret_points = Point2f[], ret_colors = []; color=:black, branch_shape=:tree, groups=nothing)
    isnothing(node.children) && return nothing
    child1 = nodes[node.children[1]]
    child2 = nodes[node.children[2]]
   
    l = dendrogram_connectors(Val(branch_shape), node, child1, child2)
    
    # even if the inputs are 2d, the outputs should be 3d - this is what `to_ndim` does.
    append!(ret_points, Makie.to_ndim.(Point3f, l, 0))
    push!(ret_points, Point3f(NaN)) # separate segments

    if isnothing(groups)
        cgroup = 0
    else
        gs = recursive_leaf_groups(node, nodes, groups)
        @debug gs
        cgroup = length(unique(gs)) == 1 ? first(gs) : 0
    end

    @debug cgroup maxlog=2

    append!(ret_colors, [cgroup for _ in 1:length(l)])
    push!(ret_colors, NaN32) # separate segments

    recursive_dendrogram_points(child1, nodes, ret_points, ret_colors; branch_shape, groups)
    recursive_dendrogram_points(child2, nodes, ret_points, ret_colors; branch_shape, groups)
    return ret_points, ret_colors
end


function Makie.plot!(plot::Dendrogram{<: Tuple{<: Dict{<: Integer, <: Union{DNode{2}, DNode{3}}}}})
    args = @extract plot (color, groups)

    points_vec = Observable{Vector{GeometryBasics.Point{2, Float32}}}([[0,0]])
    colors_vec = Observable{Any}([NaN])

    length(plot[1][])>1 && lift(plot[1], plot.branch_shape, plot[:color]) do nodes, branch_shape, color
        # this pattern is basically first updating the values of the observables,
        points_vec.val, colors_vec.val = recursive_dendrogram_points(nodes[maximum(keys(nodes))], nodes; color, branch_shape, groups=groups.val)
        # then propagating the signal, so that there is no error with differing lengths.
        notify(points_vec); notify(colors_vec)
    end

   
    lines!(plot, points_vec; color = colors_vec, colormap = plot.colormap, colorrange = plot.colorrange, linewidth = plot.linewidth, inspectable = plot.inspectable, xautolimits = plot.xautolimits, yautolimits = plot.yautolimits) 
end


# branching styles

function dendrogram_connectors(::Val{:tree}, parent, child1, child2)
    return [child1.position, parent.position, child2.position]
end

function dendrogram_connectors(::Val{:box}, parent::DNode{2}, child1::DNode{2}, child2::DNode{2})
    yp = parent.position[2]
    x1 = child1.position[1]
    x2 = child2.position[1]

    return Point2f[(x1, child1.position[2]), (x1, yp), (x2, yp), (x2, child2.position[2])]
end

function dendrogram_connectors(::Val{:box}, parent::DNode{3}, child1::DNode{3}, child2::DNode{3})
    yp = parent.position[2]
    x1 = child1.position[1]
    x2 = child2.position[1]

    return Point3f[
        (x1, child1.position[2], child1.position[3]), 
        (x1, yp, (parent.position[3] + child1.position[3])./2), 
        (x2, yp, (parent.position[3] + child2.position[3])./2), 
        (x2, child2.position[2], child2.position[3])
    ]
end


# convert utils

function find_merge(n1, n2; height=1)
    newx = min(n1[1], n2[1]) + abs((n1[1] - n2[1])) / 2
    newy = max(n1[2], n2[2]) + height
    return Point2f(newx, newy)
end

function find_merge(n1::DNode{2}, n2::DNode{2}; height=1, index=max(n1.idx, n2.idx)+1)
    newx = min(n1.position[1], n2.position[1]) + abs((n1.position[1] - n2.position[1])) / 2
    newy = max(n1.position[2], n2.position[2]) + height

    return DNode{2}(index, Point2f(newx, newy), (n1.idx, n2.idx))
end

function find_merge(n1::DNode{3}, n2::DNode{3}; height=1, index=max(n1.idx, n2.idx)+1)
    newx = min(n1.position[1], n2.position[1]) + abs((n1.position[1] - n2.position[1])) / 2
    newy = max(n1.position[2], n2.position[2]) + height
    newz = min(n1.position[3], n2.position[3]) + abs((n1.position[3] - n2.position[3])) / 2

    return DNode{3}(index, Point3f(newx, newy, newz), (n1.idx, n2.idx))
end

function Makie.convert_arguments(::Type{<: Dendrogram}, leaves::Vector{<: Point}, merges::Vector{<: Tuple{<: Integer, <: Integer}})
    nodes = Dict(i => DNode(i, n, nothing) for (i,n) in enumerate(leaves))
    nm = maximum(keys(nodes))

    for m in merges
        nm += 1
        nodes[nm] = find_merge(nodes[m[1]], nodes[m[2]]; index = nm)
    end
    return (nodes,)
end


function hcl_nodes(hcl; useheight=false)
    nleaves = length(hcl.order)
    nodes = Dict(i => DNode(i, Point2f(x, 0), nothing) for (i,x) in enumerate(invperm(hcl.order)))
    nm = maximum(keys(nodes))

    for (m1, m2) in eachrow(hcl.merges)
        nm += 1
        
        m1 = m1 < 0 ? -m1 : m1 + nleaves
        m2 = m2 < 0 ? -m2 : m2 + nleaves
        nodes[nm] = find_merge(nodes[m1], nodes[m2]; index=nm)            
    end

    return nodes
end

function recursive_leaf_groups(node, nodes, groups)
    if isnothing(node.children)
        return [groups[node.idx]]
    else
        return vcat(
            recursive_leaf_groups(nodes[node.children[1]], nodes, groups),
            recursive_leaf_groups(nodes[node.children[2]], nodes, groups)
            )
        end
end

nodes = hcl_nodes(hcl_taxa_prevalences)
(points, colors) = recursive_dendrogram_points(nodes[maximum(keys(nodes))], nodes; color=:black, branch_shape=:box)

axDend = Axis(
    figure3_master[1,4],
    xticks = (eachindex(bins), ["2", "", "4", "", "6", "", "8", "", "10", "", "12", "", "14"]),
    yticks = (eachindex(important_bugs), replace.(important_bugs[prevalence_order], "_" => " ")),
    yticklabelfont="TeX Gyre Heros Makie Italic",
    yticklabelsize=24,
    xticklabelsize=24,
    ylabelsize=24,
    xlabelsize=24,
    titlesize = 24
)
hidedecorations!(axDend)
hidespines!(axDend)
tightlimits!(axDend)
# ylims!(axDend, (0.5, 34.5))
ylims!(axDend, (34.5, 0.5))

lines!(axDend, Point2f.([[x[2],x[1]] for x in points]); color=:black)
# lines!(axDend, Point2f.([[x[2],-x[1]] for x in points]); color=:black)
```

## Add labels
```julia
Label(figure3_master[1, 1, TopLeft()], "a", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(figure3_master[1, 2, TopLeft()], "b", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
Label(figure3_master[1, 3, TopLeft()], "c", fontsize = 22, font = :bold, padding = (0, 5, 5, 0), halign = :right, alignmode = Inside())
```

## Fix layout
```julia
colsize!(figure3_master.layout, 1, Relative(0.29))
colsize!(figure3_master.layout, 2, Relative(0.29))
colsize!(figure3_master.layout, 3, Relative(0.29))
colsize!(figure3_master.layout, 4, Relative(0.13))
```

# Export Figure 2

```julia
save(joinpath(outdir, "figures", "Figure3.png"), figure3_master)
save(joinpath(outdir, "figures", "Figure3.eps"), figure3_master)
save(joinpath(outdir, "figures", "Figure3.svg"), figure3_master)
save(joinpath(outdir, "figures", "Figure3.pdf"), figure3_master)
figure3_master
```