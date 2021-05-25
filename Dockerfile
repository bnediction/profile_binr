FROM colomoto/colomoto-docker:2021-04-01

USER root

RUN conda install -y -c conda-forge \
    'r-base>=4.0'\
    r-bigmemory\
    r-diptest\
    r-doSNOW\
    r-dplyr\
    r-foreach\
    r-glue\
    r-magrittr\
    r-mclust\
    r-moments\
    r-tibble\
    r-tidyr

    #r-docstring\

COPY . /src
RUN pip install --no-deps /src

USER user
