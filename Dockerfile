FROM condaforge/mambaforge

LABEL maintainer="Jan Detleffsen <jan.detleffsen@mpi-bn.mpg.de>"

COPY . /tmp/

# Set the time zone (before installing any packages)
RUN echo 'Europe/Berlin' > apt-get install -y tzdata

# update container
RUN apt-get update --assume-yes

# install git to check for file changes
RUN apt-get install -y git

# create a non-root user
RUN useradd --no-log-init -r -g users user

# setup home directory
WORKDIR /home/user

# change owner and permissions to all users
RUN chown -R :users /opt && \
    chown -R :users /tmp && \
    chown -R :users /home/user && \
    chmod -R 775 /opt && \
    chmod -R 775 /tmp && \
    chmod -R 775 /home/user

# switch to non-root user for safer operations
USER user

# update mamba
RUN mamba update -n base mamba && \
    mamba --version

# install enviroment
RUN mamba env update -n base -f /tmp/peakqc_env.yml

# install sctoolbox
RUN pip install "/tmp/" && \
    pip install pytest && \
    pip install pytest-cov && \
    pip install pytest-html && \
    pip install pytest-mock

# change user to cleanup and install ssh
USER root 

# clear tmp
RUN rm -r /tmp/*

# Generate an ssh key
RUN apt-get install -y openssh-client && \
    mkdir .ssh && \
    ssh-keygen -t ed25519 -N "" -f .ssh/id_ed25519

# Switch to non root default user 
USER user
