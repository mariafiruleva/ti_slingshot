FROM dynverse/dynwrapr:v0.1.0

ARG GITHUB_PAT

RUN apt-get update && apt-get install -y libcgal-dev libglu1-mesa-dev libgsl-dev

RUN R -e 'install.packages("irlba")' && \
	R -e 'devtools::install_cran(c("Rcpp", "RcppEigen", "RSpectra"), verbose = TRUE, quiet = TRUE)' && \
	R -e 'devtools::install_github("kstreet13/slingshot")' && \
	rm -rf /tmp/* # clean up temp files

COPY definition.yml run.R example.sh /code/

ENTRYPOINT ["/code/run.R"]
