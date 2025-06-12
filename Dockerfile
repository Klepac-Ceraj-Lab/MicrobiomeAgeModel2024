# Start from the official Julia image
FROM julia:1.10

# Add the Github Token sourced from dotenv as a build-arg 
ARG GITHUB_TOKEN
ENV GITHUB_TOKEN=${GITHUB_TOKEN}

RUN apt-get update && apt-get -y install git

RUN git config --global \
  url."https://${GITHUB_TOKEN}@github.com/".insteadOf "https://github.com/"

# Set working directory
WORKDIR /api

# Copy project files
COPY . /api

# Develop and install dependencies
RUN julia --project=. -e 'import Pkg; \
  Pkg.develop(Pkg.PackageSpec( \
    url="https://github.com/Klepac-Ceraj-Lab/Leap.jl.git" \
  )); \
  Pkg.develop(Pkg.PackageSpec( \
    url="https://github.com/Klepac-Ceraj-Lab/BiobakeryUtils.jl.git", \
  )); \
  Pkg.precompile();'

# Expose the port Julia server listens on
EXPOSE 1025

# Default command: launch HTTP server
CMD ["julia", "--project=.", "start_server.jl"]