# Base image
FROM rocker/shiny:4.4.1

# General updates
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev libssl-dev libharfbuzz-dev libfribidi-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install required R packages manually
RUN Rscript -e 'install.packages(c("shiny","circlize","RColorBrewer","dplyr","colorspace","BiocManager"), dependencies = TRUE)'
RUN Rscript -e 'BiocManager::install(c("ComplexHeatmap","InteractiveComplexHeatmap"), ask = FALSE, update = FALSE)'

# Copy the app files (app.R + data)
RUN rm -rf /srv/shiny-server/*
COPY /app/ /srv/shiny-server/

# Ensure that the expected user is present in the container
RUN if id shiny &>/dev/null && [ "$(id -u shiny)" -ne 999 ]; then \
        userdel -r shiny; \
        id -u 999 &>/dev/null && userdel -r "$(id -un 999)"; \
    fi; \
    useradd -u 999 -m -s /bin/bash shiny; \
    chown -R shiny:shiny /srv/shiny-server/ /var/lib/shiny-server/ /var/log/shiny-server/

USER shiny
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
