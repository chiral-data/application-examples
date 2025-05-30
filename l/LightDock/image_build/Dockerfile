FROM python 

RUN apt-get update && \
    apt-get install -y python3-numpy

# Install system dependencies
# RUN apt-get update && apt-get install -y build-essential
# RUN apt-get update && apt-get install -y --no-install-recommends \
#     wget \
#     && apt-get clean \
#     && rm -rf /var/lib/apt/lists/*

# # Create a working directory
# WORKDIR /workspace

# # Copy the launching script run.sh from the project folder @common
# COPY run_dok.sh /workspace/run_dok.sh

# RUN pip install lightdock

# Install by compiling
RUN git clone https://github.com/lightdock/lightdock.git
WORKDIR /lightdock
RUN ./setup.sh

ENTRYPOINT ["/bin/sh"]
