FROM ubuntu:24.10

# Install system dependencies
RUN apt-get update && apt-get -y --no-install-recommends install build-essential \
    gfortran \
    ca-certificates \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy sievgene pack
ARG YYMMDD="240118"
RUN mkdir -p /opt/mypresto
COPY sievgene_pack_${YYMMDD}.tar.gz /opt/mypresto/sievgene_pack_${YYMMDD}.tar.gz
WORKDIR /opt/mypresto
RUN tar xzvf sievgene_pack_${YYMMDD}.tar.gz

# Compile sievgene
WORKDIR /opt/mypresto/sievgene_pack_${YYMMDD}
RUN bin/install.sh
RUN bash bin/test_all.sh
RUN bash bin/test_all.sh screening

# Copy the launching script run_dok.sh
COPY run_dok.sh ./run_dok.sh 

ENTRYPOINT ["/bin/bash"]



