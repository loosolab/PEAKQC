FROM condaforge/mambaforge

LABEL maintainer="Jan Detleffsen <jan.detleffsen@mpi-bn.mpg.de>"

COPY . /tmp/

# Set the time zone (before installing any packages)
RUN echo 'Europe/Berlin' > apt-get install -y tzdata

# update container
RUN apt-get update --assume-yes

# install git to check for file changes
RUN apt-get install -y git

# install other dependencies
RUN apt-get install -y make

# update mamba
RUN mamba update -n base mamba && \
    mamba --version

# install enviroment
RUN mamba env update -n base -f /tmp/peakqc_env.yml

# install sctoolbox
RUN pip install "/tmp/" --group "/tmp/pyproject.toml:test"

# clear tmp
RUN rm -r /tmp/*

# Generate an ssh key
RUN apt-get install -y openssh-client && \
    mkdir .ssh && \
    ssh-keygen -t ed25519 -N "" -f .ssh/id_ed25519
