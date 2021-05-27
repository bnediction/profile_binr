FROM colomoto/colomoto-docker:2021-04-01

USER root

RUN conda install -y colomoto::profile_binr==0.1.2

RUN rm -rf /notebook/*
COPY --chown=user:user *.ipynb /notebook/

USER user
