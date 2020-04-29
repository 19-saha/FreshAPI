FROM rocker/r-base

RUN apt-get update -qq && apt-get install -y \
  bash \
  curl \
  pandoc \
  pandoc-citeproc \
  git-core \
  libssl-dev \
  libcurl4-gnutls-dev


RUN R -e "install.packages(c('plumber', 'caret', 'jsonlite', 'stats', 'kernlab', 'rmarkdown', 'tinytex', 'openxlsx', 'prospectr', 'readxl'))"

RUN R -e "tinytex::install_tinytex()"

COPY plumber-api.R /plumber-api.R
COPY plumber-api_router.R /plumber-api_router.R
COPY models.xlsx /models.xlsx
COPY functions.R /functions.R
COPY models /models
COPY Test Data /Test Data
COPY MSI_Report.Rmd /MSI_Report.Rmd
COPY FTIR_CB_Report.Rmd /FTIR_CB_Report.Rmd
COPY FTIR_CTF_Report.Rmd /FTIR_CTF_Report.Rmd

VOLUME /Fresh-API

ENTRYPOINT ["R", "-e", "source('plumber-api_router.R')"]




