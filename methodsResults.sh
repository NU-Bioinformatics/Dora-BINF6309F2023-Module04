#!/usr/bin/env bash
# methodsResults.sh

module load singularity
singularity run -B "/scratch:/scratch,/courses:/courses"

# Render RMarkdown to HTML
/shared/container_repository/rstudio/rocker-geospatial-4.2.1.sif Rscript -e "rmarkdown::render('methodsResults.Rmd')"
