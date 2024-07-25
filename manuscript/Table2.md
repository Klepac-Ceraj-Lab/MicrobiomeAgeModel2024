## Import the necessary packages
```julia
using XLSX
using CSV
using DataFrames
using Chain
using CategoricalArrays
using Statistics
```

## Loading the Sample metadata and Demographics
```julia
demo_xlsx = XLSX.readxlsx("khula_demographics.xlsx")
demo_df = DataFrame(Tables.columntable(XLSX.gettable(demo_xlsx.workbook.sheets[1]; infer_eltypes=true)))
demo_df.subject_id = "khula-" .* demo_df.subject_id

khula_samples = @chain CSV.read("results/2024AgeModelManuscript/combined_inputs.csv", DataFrame) begin
    khula_bugvec = subset(:datasource => x -> x .== "1kDLEAP-KHULA")
    select([:subject_id, :sample, :ageMonths, :visit])
end 
```

### Keep only demographics for the participants with metagenoms in this study
```julia
initial_rows = nrow(demo_df)
println("Initial number of rows in Demographics table: $(initial_rows)")
subset!(demo_df, :subject_id => x -> x .∈ Ref(khula_samples.subject_id))
total_rows = nrow(demo_df)
println("Final number of rows in Demographics table: $(total_rows)")
```

## Renaming some columns to make our like easier
_(It ended up being only one...)_
```julia
@chain demo_df begin
    rename!("child_sex_0female_1male_UPDATEDfeb2024" => "child_sex")
end
```
Remember that 0 is female, 1 is male!

## Ensuring order for the Income
```julia
demo_df.householdIncomebyCurrency = string.(demo_df.householdIncomebyCurrency)
income_levels = [ "Unknown", "Less than R1000 per month", "R1000-R5000 per month", "R5000-R10 000 per month", "More than R10 000 per month" ]
demo_df.householdIncomebyCurrency = CategoricalArray(demo_df.householdIncomebyCurrency, ordered=true, levels=income_levels)

## Selecting the columns of interest
```julia
columns_of_interest = [
    "country",
    "child_sex",
    "ethnicity",
    "spokenLanguage",
    "householdIncomebyCurrency",
    "prepostnatal",
    "mom_edu_en"
]
```

# Initialize an empty DataFrame to store the combined summary data
```julia
combined_summary = DataFrame()
```

# Summarize each column and calculate percentages
```julia
for col in columns_of_interest
    summary = combine(groupby(demo_df, col), nrow => :count)
    summary.percentage = 100 * summary[:, :count] ./ total_rows
    summary.column .= string(col)
    rename!(summary, col => :value)
    (col == "householdIncomebyCurrency") && (summary.value = CategoricalArray(summary.value, ordered=true, levels=income_levels))
    sort!(summary, :value)
    
    # Append the summary to the combined DataFrame
    combined_summary = vcat(combined_summary, summary)
end

# Rearrange columns for better readability
combined_summary = select(combined_summary, :column, :value, :count, :percentage)
```

## Dealing with EPDS
```julia
epds_data = demo_df[:, "updated_epds_score_en"]
mean_sd = DataFrame(
    column = "EPDS",
    value = "Mean (SD)",
    count = mean(epds_data),
    percentage = std(epds_data)
)
    
median_min_max = DataFrame(
    column = "EPDS",
    value = "Median [Min, Max]",
    count = median(epds_data),
    percentage = "[$(minimum(epds_data)), $(maximum(epds_data))]"
)
    
combined_summary = vcat(combined_summary, mean_sd, median_min_max)
```

# Display the combined summary table
```julia
using PrettyTables
pretty_table(combined_summary)
CSV.write("manuscript/Table2.csv", combined_summary)
```