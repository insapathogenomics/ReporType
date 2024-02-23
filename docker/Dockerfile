# Use the official Ubuntu base image
FROM ubuntu:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary packages
RUN apt-get update && \
    apt-get install -y \
    git \
    build-essential \
    python3 \
    python3-pip \
	wget \
    bzip2

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Add Conda to the system PATH
ENV PATH="/opt/conda/bin:/app/ReporType:${PATH}"

# Clone the ReporType repository
RUN git clone https://github.com/insapathogenomics/ReporType /app/ReporType

# Set working directory
WORKDIR /app/ReporType
RUN mkdir -p data_shared

COPY show_databases.sh  .

## install ReportType
RUN chmod +x *.sh; ./install.sh

### Activate ReportType (need to downgrade pulp because of Conda)
RUN conda init bash
RUN eval "$(/opt/conda/bin/conda shell.bash hook)" && conda activate ReporType && pip install --force-reinstall -v "pulp==2.7" && chmod +x databases_install.sh && ./databases_install.sh
 
COPY ReporType /usr/bin
RUN chmod u+x /usr/bin/ReporType

## keep the container running
#ENTRYPOINT ["bash"]
#CMD ["-f","/dev/null"]
#CMD ["bash"]
