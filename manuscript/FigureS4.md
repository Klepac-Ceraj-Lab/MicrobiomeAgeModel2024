# This is the leftover code from previously-planned Figure S4 that did not make it into the final submission.

## Data for Supplementary Figure 4 - Not included in final manuscript!
```julia
function myload(::UnirefProfiles; timepoint_metadata = load(Metadata()))
    comm = Leap.read_arrow(inputfiles("genefamilies.arrow"); featurefunc = genefunction)
    insert!(comm, timepoint_metadata; namecol=:sample)
    return comm[:, timepoint_metadata.sample]
end

unirefs = myload(UnirefProfiles(); timepoint_metadata = functional_mdata) # this can take a bit

previndex_ecs = vec(prevalence(filtered_functional_profiles[:, 1:length(samples(filtered_functional_profiles))]) .> 0.05)
prevfiltered_ecs = filtered_functional_profiles[previndex_ecs, 1:length(samples(filtered_functional_profiles))]
previndex_unirefs = vec(prevalence(unirefs[:, 1:length(samples(unirefs))]) .> 0.05)
prevfiltered_unirefs = unirefs[previndex_unirefs, 1:length(samples(unirefs))]

unirefs_n_features = unirefs.abundances.m
unirefs_n_samples = unirefs.abundances.n
unirefs_n_cells = unirefs_n_features * unirefs_n_samples
unirefs_n_nzval = length(unirefs.abundances.nzval)
unirefs_nzval_proportion = unirefs_n_nzval/unirefs_n_cells
unirefs_featbysamp = unirefs_n_features / unirefs_n_samples
unirefs_logfeatbysamp = log10(unirefs_featbysamp)

prevfiltered_unirefs_n_features = prevfiltered_unirefs.abundances.m
prevfiltered_unirefs_n_samples = prevfiltered_unirefs.abundances.n
prevfiltered_unirefs_n_cells = prevfiltered_unirefs_n_features * prevfiltered_unirefs_n_samples
prevfiltered_unirefs_n_nzval = length(prevfiltered_unirefs.abundances.nzval)
prevfiltered_unirefs_nzval_proportion = prevfiltered_unirefs_n_nzval/prevfiltered_unirefs_n_cells
prevfiltered_unirefs_featbysamp = prevfiltered_unirefs_n_features / prevfiltered_unirefs_n_samples
prevfiltered_unirefs_logfeatbysamp = log10(prevfiltered_unirefs_featbysamp)

ecs_n_features = filtered_functional_profiles.abundances.m
ecs_n_samples = filtered_functional_profiles.abundances.n
ecs_n_cells = ecs_n_features * ecs_n_samples
ecs_n_nzval = length(filtered_functional_profiles.abundances.nzval)
ecs_nzval_proportion = ecs_n_nzval/ecs_n_cells
ecs_featbysamp = ecs_n_features / ecs_n_samples
ecs_logfeatbysamp = log10(ecs_featbysamp)

prevfiltered_ecs_n_features = prevfiltered_ecs.abundances.m
prevfiltered_ecs_n_samples = prevfiltered_ecs.abundances.n
prevfiltered_ecs_n_cells = prevfiltered_ecs_n_features * prevfiltered_ecs_n_samples
prevfiltered_ecs_n_nzval = length(prevfiltered_ecs.abundances.nzval)
prevfiltered_ecs_nzval_proportion = prevfiltered_ecs_n_nzval/prevfiltered_ecs_n_cells
prevfiltered_ecs_featbysamp = prevfiltered_ecs_n_features / prevfiltered_ecs_n_samples
prevfiltered_ecs_logfeatbysamp = log10(prevfiltered_ecs_featbysamp)

taxa_n_features = ncol(taxonomic_profiles[:, 11:end-1])
taxa_n_samples = nrow(taxonomic_profiles[:, 11:end-1])
taxa_n_cells = taxa_n_features * taxa_n_samples
taxa_n_nzval = sum(vec(Matrix(taxonomic_profiles[:, 11:end-1])) .> 0.0)
taxa_nzval_proportion = taxa_n_nzval/taxa_n_cells
taxa_featbysamp = taxa_n_features / taxa_n_samples
taxa_logfeatbysamp = log10(taxa_featbysamp)
```