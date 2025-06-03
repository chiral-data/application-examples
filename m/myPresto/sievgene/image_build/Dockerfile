FROM ubuntu:24.10

RUN apt-get update && apt-get -y install build-essential gfortran

ARG YYMMDD="240118"

RUN mkdir -p /opt/mypresto
COPY sievgene_pack_${YYMMDD}.tar.gz /opt/mypresto/sievgene_pack_${YYMMDD}.tar.gz
WORKDIR /opt/mypresto
RUN tar xzvf sievgene_pack_${YYMMDD}.tar.gz

WORKDIR /opt/mypresto/sievgene_pack_${YYMMDD}
RUN bin/install.sh
RUN bash bin/test_all.sh
RUN bash bin/test_all.sh screening

ENTRYPOINT ["/bin/bash"]

