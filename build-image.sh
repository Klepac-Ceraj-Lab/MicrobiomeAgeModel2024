#!/bin/bash

source .env

docker build \
  --build-arg GITHUB_TOKEN=$GITHUB_TOKEN \
  . -t microbiomeagemodelv1:latest

docker run -d -p 1025:1025 microbiomeagemodelv1:latest