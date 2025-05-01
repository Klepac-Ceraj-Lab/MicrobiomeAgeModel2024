#####
# Microbiome Age Model Server logic
# Authors: Guilherme Fahur Bottino, Kevin S. Bonham, Vanja Klepac-Ceraj
#####

using HTTP, JSON3, JLD2, Dates

using Chain, MultivariateStats, Distances, UUIDs
using Diversity, Random, KernelDensity, Statistics
using CairoMakie, CategoricalArrays, GLM, StatsBase
using StableRNGs, Polynomials, CSV, DataToolkit

using BiobakeryUtils, Microbiome
using Leap
using MicrobiomeAgeModel2024

const REQUESTS_ROOT = "runs"  # ← change to wherever you want

# 1) load your trained RegressionProbeData
JLD2.@load "results/2025MaaSDev/AgeModel_FullCV_Results.jld"
# JLD2.@load "AgeModel_FullCV_Results.jld"

function _dispatch_tsv(path::String)
    # find the real header row and split it
    header_line = ""
    open(path) do io
        for raw in eachline(io)
            s = strip(raw)
            # the actual column‐names line can be commented or not
            if startswith(s, "#clade_name") || startswith(s, "clade_name")
                header_line = startswith(s, "#") ? s[2:end] : s
                break
            end
        end
    end

    header_cols = split(header_line, '\t')

    if "additional_species" in header_cols ## Meaning, a single-sample
        sample_name = split(readlines(path)[3])[end] ## Assuming standard MP3 output
        mp3_profiles = @chain BiobakeryUtils.metaphlan_profile(path; sample = sample_name) begin
            filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
            filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
            comm2wide()
            # select(Not([:sample_base, :file]))    
        end
        return mp3_profiles
    else ## Meaning, a multi-sample merged file
        mp3_profiles = @chain BiobakeryUtils.metaphlan_profiles(path) begin
            filter(t-> !ismissing(taxrank(t)), _[:, samplenames(_)])
            filter(t-> taxrank(t) == :species, _[:, samplenames(_)])
            comm2wide()
            # select(Not([:file]))
        end
        return mp3_profiles
    end

end

function prediction_handler(req::HTTP.Request)

    this_uuid = string(UUIDs.uuid4())

    try

        # 0. sanity check on REQ method
        if req.method != "POST"
            return HTTP.Response(
                405,
                ["Content-Type" => "application/json"];
                body = JSON3.write(Dict(
                    :message => "Only POST requisitions allowed",
                    :uuid => this_uuid
                ))
            )
        end
        
        # 1. parse all the parts
        parts = try

            pts = HTTP.parse_multipart_form(req)
            println(pts) # for debugging

            pts

        catch e

            @error "Parsing failed for multipart request" exception = e uuid = this_uuid
            return HTTP.Response(
                500,
                ["Content-Type" => "application/json"];
                body = JSON3.write(Dict(
                    :message => "Parsing failed for multipart request on job $(this_uuid) at $(Dates.now()).",
                    :uuid => this_uuid
                ))
            )
            
        end

        # 2. separate text fields vs. file uploads
        params, uploads = try

            p = Dict{String,String}()
            u = Dict{String,Tuple{String,Vector{UInt8}}}()

            for part in parts
                raw = read(part.data)                # Vector{UInt8}
                if isnothing(part.filename)
                    # text field
                    p[part.name] = String(raw)
                else
                    # file field
                    u[part.name] = (part.filename, raw)
                end
            end
            
            # println(params) # for debugging
            # println(keys(uploads)) # for debugging
            (p,u)

        catch e
            @error "Argument handling failed for multipart request" exception = e uuid = this_uuid
            return HTTP.Response(
                500,
                ["Content-Type" => "application/json"];
                body = JSON3.write(Dict(
                    :message => "Argument handling failed for multipart request on job $(this_uuid) at $(Dates.now()).",
                    :uuid => this_uuid
                ))
            )
        end

        # 3. generating the UUID for the folder and creating the folder
        # try
            this_dir = joinpath(REQUESTS_ROOT, this_uuid)
            mkpath(this_dir)

            # 4. save params.json
            open(joinpath(this_dir, "params.json"), "w") do io
                JSON3.write(io, params; indent=2)
            end

            # 5. save the uploaded files
            for (_field, (fname, bytes)) in uploads
                open(joinpath(this_dir, fname), "w") do io
                    write(io, bytes)
                end
            end
        # catch
        #     return HTTP.Response(
        #         500,
        #         ["Content-Type" => "application/json"];
        #         body = JSON3.write(Dict(
        #             :message => "File I/O failed for write operations on job $(this_uuid) at $(Dates.now()).",
        #             :uuid => this_uuid
        #         ))
        #     )
        # end

        # 6. Parseing the Metaphlan profiles
        # try 
            mp3_profiles = _dispatch_tsv(joinpath(this_dir, "mp3_profiles.tsv"))
            mp3_profiles.Shannon_index = map(x -> Microbiome.shannon(collect(x)), eachrow(mp3_profiles[:, 2:ncol(mp3_profiles)]))
            print(mp3_profiles) # for debugging
        # catch
        #     return HTTP.Response(
        #         500,
        #         ["Content-Type" => "application/json"];
        #         body = JSON3.write(Dict(
        #             :message => "Failed to read MetaPhlAn outputs for job $(this_uuid) at $(Dates.now()). Check internal file format.",
        #             :uuid => this_uuid
        #         ))
        #     )
        # end

        # 7. Prediction time!!!
        # try
            preds = predict_regression_runtime(regression_Age_FullCV, mp3_profiles)
            print(preds) # for debugging
            CSV.write(joinpath(this_dir, "predictions.csv"), preds)
        # catch
        #     return HTTP.Response(
        #         500,
        #         ["Content-Type" => "application/json"];
        #         body = JSON3.write(Dict(
        #             :message => "Failed to run prediction for $(nrow(mp3_profiles)) samples from job $(this_uuid) at $(Dates.now()).",
        #             :uuid => this_uuid
        #         ))
        #     )
        # end

        # 8. If all goes well, respond with 200 and success message.
        return HTTP.Response(
            200,
            ["Content-Type" => "application/json"];
            body = JSON3.write(Dict(
                :message => "Successful run for $(nrow(mp3_profiles)) samples from job $(this_uuid) at $(Dates.now()).",
                :uuid => this_uuid
            ))
        )

    catch e
        @error "Argument handling failed for multipart request" exception = e uuid = this_uuid
        return HTTP.Response(
            500,
            ["Content-Type" => "application/json"];
            body = JSON3.write(Dict(
                :message => "Error on job $(this_uuid) at $(Dates.now()).",
                :uuid => this_uuid
            ))
        )
    end
end

### download endpoint
function download_handler(req::HTTP.Request)
    # only GET allowed
    if req.method != "GET"
        return HTTP.Response(405; body="Only GET allowed")
    end

    target = String(req.target)
    print(target)
    parts = split(target, '/')
    print(parts)
    if length(parts) != 4 || parts[3] != "download"
        return HTTP.Response(400; body="Bad target path")
    end
    uuid = parts[end]

    # locate the CSV
    csv_path = joinpath(REQUESTS_ROOT, uuid, "predictions.csv")
    if !isfile(csv_path)
        return HTTP.Response(404,
            ["Content-Type" => "application/json"];
            body    = JSON3.write(Dict(
              "message" => "No predictions file found for UUID $uuid",
              "uuid"    => uuid
            ))
        )
    end

    # stream it back with download headers
    csv_text = read(csv_path, String)
    return HTTP.Response(
        200,
        [
            "Content-Type"        => "text/csv",
            "Content-Disposition" => "attachment; filename=\"predictions_$uuid.csv\""
        ];
        body = csv_text
    )
end

## Define REST endpoints
application_router = HTTP.Router()

HTTP.@register(application_router, "POST", "/age_model_v1/prediction", prediction_handler)
HTTP.@register(application_router, "GET", "/age_model_v1/download", download_handler)


HTTP.serve(application_router, "0.0.0.0", 1025)