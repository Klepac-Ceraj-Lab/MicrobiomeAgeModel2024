# Auxiliary scripts to generate the `ecs.arrow` file for reproducibility purposes

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
using MicrobiomeAgeModel2024

### Configurable parameters and notebook set-up
DataToolkit.loadcollection!("./Data_Local.toml")    ## Uncomment this line to use local files located on the "data" subfolder and the Local relative filesystem references

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

function my_write_gfs_arrow(gfs_rootdir, arrow_outdir; kind="genefamilies", names=false, stratified=false, sample_filter=nothing)
    stripper = "_$kind" * (names ? "_rename.tsv" : ".tsv")
    filt = Regex(string(raw"SEQ\d+_S\d+_", kind))
    df = DataFrame(file = filter(f-> contains(f, filt), readdir(gfs_rootdir, join=true)))
    df.sample = map(s-> replace(s, stripper => ""), basename.(df.file))
    df.seqid = map(s-> replace(s, r"_S\d+$"=> ""), df.sample)
    isnothing(sample_filter) || subset!(df, AsTable(["sample", "seqid"]) => ByRow(s-> s.sample ∈ sample_filter || s.seqid ∈ sample_filter))
    @debug "Found $(nrow(df)) files"

    knead = Leap.load(Leap.ReadCounts())
    leftjoin!(df, select(knead, "sample_uid"=>"sample", 
                                AsTable(["final pair1", "final pair2"])=> ByRow(row-> row[1]+row[2]) =>"read_depth");
                            on="sample"
    )
    
    @info "getting features"
    feats = mapreduce(union, eachrow(df)) do row
        fs = CSV.read(row.file, DataFrame; header=["feature", "value"], skipto=2, select=[1])[!,1]
        stratified || filter!(f->!contains(f, '|'), fs) # skip stratified feats
        Set(fs)
    end
    featuremap = Dict(f=> i for (i,f) in enumerate(feats))
    
    isdir(arrow_outdir) || mkpath(arrow_outdir)

    @info "writing arrow file"
    open(joinpath(arrow_outdir, "$kind.arrow"), "w") do io
        tbls = Tables.partitioner(eachrow(df)) do row
            @debug "writing $(row.sample)"

            sdf = CSV.read(row.file, DataFrame; header=["feature", "value"], skipto=2)
            stratified || subset!(sdf, "feature"=> ByRow(f-> !contains(f, '|'))) # skip stratified feats
            sdf.sample .= row.sample
            sdf.sidx .= rownumber(row)
            sdf.fidx = ThreadsX.map(f-> featuremap[f], sdf.feature)

            sdf
        end

        Arrow.write(io, tbls; metadata=("features" => join(feats, '\n'), 
                                        "samples"  => join(df.sample, '\n'),
                                        "files"    => join(df.file, '\n'),
                                        "reads"    => join(df.read_depth, '\n')
                                        )
        )                                
    end

    return nothing
end


## Loading data
regression_Age_FullCV = d"cv_results"["regression_Age_FullCV"]
taxonomic_profiles = regression_Age_FullCV.original_df

bins = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

## Loading data
functional_mdata = @chain taxonomic_profiles begin
    khula_bugvec = subset(:datasource => x -> x .== "1kDLEAP-KHULA")
    select([:subject_id, :sample, :ageMonths, :visit])
end 

# functional_profiles = myload(ECProfiles(), open(dataset("functional_profiles"), DataToolkit.FilePath).path; timepoint_metadata = functional_mdata) # this can take a bit
filtered_functional_profiles = relativeabundance(filter(f -> name(f) == "UNMAPPED" || hastaxon(f), functional_profiles))
ecs_mdata = DataFrame(get(filtered_functional_profiles))
CSV.write("ecs_metadata.csv", ecs_mdata)

## Isolating the ECS on a specific folder with symlinks
origin_path = "/grace/sequencing/processed/mgx/humann/rename"
link_folder_path = "/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/notebooks/isolate_ecs"

# Ensure the folder for symbolic links exists
if !isdir(link_folder_path)
    mkdir(link_folder_path)
end

# Iterate over each prefix
for prefix in ecs_mdata.sample
    # Get all files in the folder
    files = readdir(origin_path)

    # Filter files that match the current prefix
    matching_files = filter(f -> startswith(f, prefix), files)

    # Create symbolic links for each matching file
    for file in matching_files
        original_file = joinpath(origin_path, file)
        symlink_path = joinpath(link_folder_path, file)

        try
            # Create symbolic link
            symlink(original_file, symlink_path)
            println("Symbolic link created for: $original_file -> $symlink_path")
        catch e
            println("Failed to create symbolic link for: $original_file. Error: $e")
        end
    end
end

## Writing arrow file
my_write_gfs_arrow(
    "/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/notebooks/isolate_ecs",
    "/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/notebooks";
    kind="ecs",
    names=true,
    stratified=true,
    sample_filter=nothing
)

## Testing the new file
functional_profiles = myload(ECProfiles(), "/home/guilherme/.julia/dev/MicrobiomeAgeModel2024/notebooks/ecs.arrow"; timepoint_metadata = functional_mdata) # this can take a bit

t1_samples = subset(ecs_mdata, :visit => x -> x .== "3mo")
println("For t1 samples, M = $(round(mean(t1_samples.ageMonths); digits = 2)), SD = $(round(std(t1_samples.ageMonths); digits = 2)), range = $(extrema(t1_samples.ageMonths))")
t2_samples = subset(ecs_mdata, :visit => x -> x .== "6mo")
println("For t2 samples, M = $(round(mean(t2_samples.ageMonths); digits = 2)), SD = $(round(std(t2_samples.ageMonths); digits = 2)), range = $(extrema(t2_samples.ageMonths))")
t3_samples = subset(ecs_mdata, :visit => x -> x .== "12mo")
println("For t3 samples, M = $(round(mean(t3_samples.ageMonths); digits = 2)), SD = $(round(std(t3_samples.ageMonths); digits = 2)), range = $(extrema(t3_samples.ageMonths))")

longitudinal_samples = innerjoin(t1_samples, t3_samples, on = :subject_id, makeunique=true)