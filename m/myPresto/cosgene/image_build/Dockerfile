FROM ubuntu:24.10

# Install system dependencies
RUN apt-get update && apt-get -y --no-install-recommends install build-essential \
    gfortran \
    ca-certificates \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy cosgene pack
ARG YYMMDD="240724"
RUN mkdir -p /opt/mypresto
COPY cosgene_pack_${YYMMDD}.tar.gz /opt/mypresto/cosgene_pack_${YYMMDD}.tar.gz
WORKDIR /opt/mypresto
RUN tar xzvf cosgene_pack_${YYMMDD}.tar.gz

# Compile cosgene
WORKDIR /opt/mypresto/cosgene_pack_${YYMMDD}
RUN bin/install.sh

# Copy the launching script run_dok.sh
COPY run_dok.sh ./run_dok.sh 

ENTRYPOINT ["/bin/bash"]



