FROM chiral.sakuracr.jp/mypresto_sievgene:latest

# Install system dependencies and "wget"
RUN apt-get update && apt-get install -y build-essential
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# enter the working directory
WORKDIR /opt/mypresto/sievgene_pack_FreeBSD_240118

# Copy the launching script run.sh from the project folder @common
COPY run.sh ./run.sh 

ENTRYPOINT ["/bin/sh"]
