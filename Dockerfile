FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('shiny', 'shinythemes', 'DT', 'ggplot2', 'openxlsx', 'tibble', 'e1071', 'lmomco', 'MASS', 'survival', 'fitdistrplus', 'goftest', 'scales'), repos='https://cran.rstudio.com/')"

# Copy app files
COPY . /srv/shiny-server/
WORKDIR /srv/shiny-server/

# Expose Hugging Face Spaces port
EXPOSE 7860

# Run app on port 7860
CMD ["R", "-e", "shiny::runApp(host='0.0.0.0', port=7860)"]
