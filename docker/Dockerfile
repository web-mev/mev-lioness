FROM mambaorg/micromamba

# for building the conda env
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# Build the environment based on the yml file
RUN micromamba install -y -n base -f /tmp/environment.yml \
  && micromamba clean --all --yes

# Note that for the WDL run, these files will be addressed directly.
# However, failure of this step in the Docker build (e.g. as part of CI tools)
# will indicate that the files are missing and something is wrong.
# So this serves as a double-check.
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN curl --silent --output /dev/null https://webmev-public.s3.us-east-2.amazonaws.com/tissues_motif.ensg.tsv
RUN curl --silent --output /dev/null https://webmev-public.s3.us-east-2.amazonaws.com/tissues_motif.symbol.tsv
RUN curl --silent --output /dev/null https://webmev-public.s3.us-east-2.amazonaws.com/tissues_ppi.tsv

ADD panda.py /usr/local/bin/
ADD lioness.py /usr/local/bin/
ADD merge_lioness.py /usr/local/bin/
ADD determine_scatter.py /usr/local/bin/