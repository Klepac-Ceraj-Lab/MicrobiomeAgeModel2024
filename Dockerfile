# start from the official Julia image
FROM julia:1.10

# set working directory
WORKDIR /app

# 1) copy project files
COPY Project.toml Manifest.toml ./

# 2) copy your local dev packages into the image
#    (make sure to run docker build from a context that includes .julia/dev)
COPY .julia/dev/Leap              /app/dev/Leap
COPY .julia/dev/BiobakeryUtils    /app/dev/BiobakeryUtils

# 3) develop and install dependencies
RUN julia --project=. -e 'using Pkg; \
      Pkg.develop(path="dev/Leap"); \
      Pkg.develop(path="dev/BiobakeryUtils"); \
      Pkg.instantiate()'

# 4) copy your server code
COPY start_server.jl ./

# 5) copy the saved model file into the container
#    replace the source path with wherever your .jld2 actually lives
COPY .julia/dev/MicrobiomeAgeModel2024/results/2025MaaSDev/AgeModel_FullCV_Results.jld ./

# 6) expose the port your Julia server listens on
EXPOSE 1025

# 7) default command: launch your HTTP server
CMD ["julia", "--project=.", "start_server.jl"]