## Master colors utilized throughout the manuscript
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

## Function to setup the output directory under the "results" folder
function setup_outdir(; experiment_name = "MicrobiomeAge2024_Reproduction")
    outdir = joinpath(pwd(), "results", experiment_name)
    figdir = joinpath(outdir, "figures")
    deepdivemonodir, deepdivecolordir = ( joinpath(figdir, "species_monocolor_scatterplots"), joinpath(figdir, "species_colored_scatterplots") )
    isdir(outdir) ? @warn("Directory $(outdir) already exists! This notebook will overwrite files already there.") : ( mkpath(outdir), mkpath(figdir), mkpath(deepdivemonodir), mkpath(deepdivecolordir) )
    return(outdir, figdir, deepdivemonodir, deepdivecolordir)
end

