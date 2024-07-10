#####
# This script is a centralized data loading file that avoids mismatches netween data loading blocks between different notebooks and reduces the number of lines
# Author: BOTTINO, G. F. 
#####
using CSV
using DataFrames
using BiobakeryUtils
using Chain

#####
# 1. CuratedMetagenomicData
#####

cmd_mdata = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/CMD/2024-06-09-cmd_samplemeta.csv", DataFrame; stringtype = String) begin
# cmd_mdata = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/cmd_samplemeta_old.csv", DataFrame; stringtype = String) begin
    select!(["study_name", "subject_id", "NCBI_accession", "infant_age", "country"])
    transform!(:infant_age => (x -> x ./ 30.5) => :infant_age)
    rename!(:infant_age => :ageMonths, :country => :site, :NCBI_accession => :sample)
    subset(:ageMonths => (x -> 2.0 .< x .< 18.0))
    transform!(:subject_id => (x -> "cmd-" .* string.(x)) => :subject_id)
end
cmd_profiles = @chain Leap.load_raw_metaphlan("/vassar/guilherme/cmd_sequences/pipeline_outputs/processed/metaphlan";replace_pattern = r"\w+_profile") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
end
cmd_pre_data = @chain innerjoin(cmd_mdata, cmd_profiles, on = :sample => :sample_base ) begin
   # TODO: remember outerjoin and 5 mismatches all just renaming samples 
    dropmissing()
    insertcols!(2, :westernized_cat => 2)    
    insertcols!(2, :datasource => "CMD")
    insertcols!(2, :datacolor => "lightblue")
    insertcols!(2, :visit => "NA")
end

print(setdiff(cmd_profiles.sample_base, cmd_mdata.sample))
print(setdiff(cmd_mdata.sample, cmd_profiles.sample_base))

#####
# 2. DIABIMMUNE
#####

diabimmune_mdata = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/DIABIMMUNE/diabimmune_combined_metadata.csv", DataFrame; stringtype = String) begin
    select!(["study_name", "subject_id", "sample", "age_at_collection", "country"])
    rename!(:age_at_collection => :ageMonths, :country => :site)
    transform!(:ageMonths => (x -> x ./ 30.5) => :ageMonths)
    subset(:ageMonths => x -> 2.0 .< x .< 18.0)
    transform!(:subject_id => (x -> "diabimmune-" .* string.(x)) => :subject_id)
end
diabimmune_profiles = @chain Leap.load_raw_metaphlan("/vassar/guilherme/cmd_sequences/pipeline_outputs/processed/metaphlan";replace_pattern = r"G\w+_profile") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
end
diabimmune_pre_data = @chain innerjoin(diabimmune_mdata, diabimmune_profiles, on = :sample => :sample_base ) begin
   # TODO: remember outerjoin and 5 mismatches all just renaming samples 
    dropmissing()
    insertcols!(2, :westernized_cat => 2)    
    insertcols!(2, :datasource => "DIABIMMUNE")
    insertcols!(2, :datacolor => "lightblue")
    insertcols!(2, :visit => "NA")
end

#####
# 3. ECHO RESONANCE
#####
echo_pre_data = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/ECHO/ECHO_clock_inputs.csv", DataFrame) begin
    select(Not(:new_seqid))
    rename!(:subject => :subject_id)
    insertcols!(1, :study_name => "BonhamK_2023")
    insertcols!(2, :westernized_cat => 1)
    subset(:ageMonths => x -> 2.0 .< x .< 18.0)
    insertcols!(2, :visit => "NA")
    insertcols!(2, :datasource => "ECHO")
    insertcols!(2, :datacolor => "purple")
    insertcols!(2, :site => "USA")
    transform!(:subject_id => (x -> "resonance-" .* string.(x) ) => :subject_id; renamecols = false)
end

#####
# 4. 1kD LEAP BRAINRISE
#####
brainrise_ages = CSV.read("/home/guilherme/Repos/Leap/ext_data/BRAINRISE/BRAINRISE_age_timepoints_months.csv", DataFrame)
brainrise_pre_data = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/BRAINRISE/BRAINRISE_combined_taxonomic_profiles.csv", DataFrame) begin
    leftjoin(brainrise_ages; on = [:subject_id => :id_estudo, :visit => :timepoint ])
    # transform!(:sample .=> (x -> string.(x)) .=> [:sample, :subject_id]; renamecols = false )
    dropmissing()
    insertcols!(1, :study_name => "NasponiliN_2024")
    transform!(:subject_id => (x -> "brainrise-" .* string.(x) ) => :subject_id; renamecols = false)
    sort!(:ageMonths)
    insertcols!(2, :datasource => "1kDLEAP-BRAINRISE")
    insertcols!(2, :datacolor => "blue")
    insertcols!(2, :site => "BRA")
    insertcols!(2, :westernized_cat => 2)
    subset(:ageMonths => x -> 2.0 .< x .< 18.0)
end

#####
# 5. 1kD LEAP CORK
#####
combine_mdata = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/CORK/manifest.csv", DataFrame; stringtype = String) begin
    select!(["individualID", "specimenID", "samplingAge", "timepoint_in_days"])
    rename!(:specimenID => :sample, :samplingAge => :ageMonths, :individualID => :subject_id, :timepoint_in_days => :visit)
    transform!(:subject_id => (x -> "cork-" .* string.(x)) => :subject_id)
end
combine_profiles = @chain Leap.load_raw_metaphlan("/home/guilherme/Repos/Leap/ext_data/CORK/metaphlan_profiles";replace_pattern = r"_profile") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
    transform!(:sample_base => (col -> map(el -> match(r"(\D+\d+)_(\D+_)?[a-zA-Z0-9.-]+", el).match, col)) => :sample_base )
end
combine_pre_data = @chain innerjoin(combine_mdata, combine_profiles, on = :sample => :sample_base ) begin
    #TODO: remember outerjoin and 5 mismatches all just renaming samples 
    dropmissing()
    insertcols!(1, :study_name => "BrittonR_2024")
    insertcols!(2, :westernized_cat => 2)
    transform!( :ageMonths => (col -> col ./ 30.5) => :ageMonths )
    subset(:ageMonths => x -> 2.0 .< x .< 18.0)
    insertcols!(2, :datasource => "1kDLEAP-COMBINE")
    insertcols!(2, :datacolor => "orange")
    insertcols!(2, :site => "IRL")
end

## 1.5. 1kD LEAP KHULA
khula_ages = CSV.read("/home/guilherme/Repos/Leap/processed_data/khula_ages_visit_reltable.csv", DataFrame)
khula_pre_data = @chain Leap.load_raw_metaphlan(;replace_pattern = r"SEQ0\d+_S\d+_profile") begin
    filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
    filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
    comm2wide()
    select(Not([:sample, :file]))
    leftjoin(khula_ages, _, on = :specimenID => :sample_base )
    #transform!(:samplingAge => (x -> x ./ 30.41) => :samplingAge; renamecols=false) #not necessary anymore, ages already in Months from Airtable.
    # select(Not([:individualID]))
    rename!(:specimenID => :sample, :samplingAge => :ageMonths, :individualID => :subject_id)
    dropmissing()
    insertcols!(2, :westernized_cat => 2)
    subset(:ageMonths => x -> 2.0 .< x .< 18.0)
    insertcols!(1, :study_name => "DonaldK_2024")
    insertcols!(2, :datasource => "1kDLEAP-KHULA")
end

## 1.6. 1kD LEAP M4EFaD
m4efad_mdata = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/M4EFaD/M4EFaD_Biospecimen_Metadata_04062024.csv", DataFrame; stringtype = String) begin
    select!(["individualID", "specimenID", "samplingAge"])
    rename!(:samplingAge => :ageMonths, :individualID => :subject_id)
    transform!(:ageMonths => (x -> x ./ 30.5) => :ageMonths)
    subset(:ageMonths => x -> 2.0 .< x .< 18.0)
    transform!(:subject_id => (x -> "m4efad-" .* string.(x)) => :subject_id)
    innerjoin(CSV.read("/home/guilherme/Repos/Leap/ext_data/M4EFaD/Biospecimen_sample_mapping.csv", DataFrame), on = :specimenID => :Sample_Name)
    select!(Not(:specimenID))
    rename!(:Seq_ID => :sample)
end

m4efad_pre_data = @chain CSV.read("/home/guilherme/Repos/Leap/ext_data/M4EFaD/attic/M4EFaD_all_taxonomic_profiles.csv", DataFrame) begin
    rename!(:ageMonths => :visit)
    innerjoin(m4efad_mdata, _ ; on = :sample => :sample)
    insertcols!(1, :study_name => "OsullivanJ_2024")
    insertcols!(2, :westernized_cat => 2)
    insertcols!(2, :datasource => "1kDLEAP-M4EFAD")
    insertcols!(2, :datacolor => "darkgreen")
    insertcols!(2, :site => "BGD")
    dropmissing()
end

### Removing MW and MAM samples
subset!(khula_pre_data, :site => x -> x .!= "Malawi")
malnourished_samples = "m4efad-" .* subset(CSV.read("/home/guilherme/Repos/Leap/ext_data/M4EFaD/attic/m4efad_ages.csv", DataFrame; stringtype = String), :Condition => x -> x .== "MAM").sample_id
subset!(m4efad_pre_data, :subject_id => x -> x .∉ Ref(malnourished_samples))

### Building the pooled dataset
combined_inputs = @chain vcat(echo_pre_data, brainrise_pre_data, combine_pre_data, khula_pre_data, m4efad_pre_data, diabimmune_pre_data, cmd_pre_data; cols=:union) begin
    transform!(_, names(_) .=> (x -> replace(x, missing => 0.0)) .=> names(_); renamecols=false)
end

combined_inputs.subject_id = String.(string.(combined_inputs.subject_id))
combined_inputs.sample = String.(combined_inputs.sample)
combined_inputs.site = String.(combined_inputs.site)
combined_inputs.visit = String.(string.(combined_inputs.visit))
combined_inputs.datasource = String.(string.(combined_inputs.datasource))
combined_inputs.datacolor = String.(string.(combined_inputs.datacolor))
combined_inputs.westernized_cat = Int64.(combined_inputs.westernized_cat)
insertcols!(combined_inputs, 2, :datagroup => map( x -> begin
        if ( x ∈ ["CMD", "DIABIMMUNE"] )
            return "CMD"
        elseif ( x == "ECHO" )
            return "ECHO"
        else
            return "LEAP"
        end
    end, combined_inputs.datasource)
)

combined_inputs = unique(combined_inputs, :sample)

for i in 1:nrow(combined_inputs)
    for j in 11:ncol(combined_inputs)
        if combined_inputs[i,j] < 0.0 
            combined_inputs[i,j] = 0.0
        elseif (presence_absence & (combined_inputs[i,j] > 0.0))
            combined_inputs[i,j] = 1.0
        end
    end
end

### Calculating Shannon Diversity and Richness
combined_inputs.richness = map(x -> sum(x .> 0.0), eachrow(Matrix(combined_inputs[:, 11:ncol(combined_inputs)])))
subset!(combined_inputs, :richness => x -> x .>= 5)
combined_inputs.shannon_index = map(x -> Microbiome.shannon(collect(x)), eachrow(combined_inputs[:, 11:ncol(combined_inputs)-1]))
select!(combined_inputs, Not(:richness))

### Exporting pooled dataset
CSV.write(joinpath(outdir, "combined_inputs.csv"), combined_inputs)