# Modified (along with supplemental files and scripts) from the Dockerimage used by FEniCS:
# https://bitbucket.org/fenics-project/docker
#
#
# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.
#
#
# Builds a base Docker image that all other TerraFERMA images are
# built from. We create the user tfuser, set up the base environment
# and default entry point and cmd.
#
# The built image can be found at:
#
#   https://github.com/orgs/TerraFERMA/packages/container/package/base
#

# build off Ubuntu 20.04LTS (Focal) with a few fixes for docker
FROM phusion/baseimage:focal-1.2.0
MAINTAINER terraferma <terraferma@lists.columbia.edu>

# Get Ubuntu updates
USER root
RUN apt-get update && \ 
    apt-get upgrade -y -o Dpkg::Options::="--force-confold" -o Dpkg::Options::="--force-confdef" && \
    apt-get -y install locales sudo && \
    apt-get clean && \
    echo "C.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Set locale environment
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV LANGUAGE=C.UTF-8

COPY docker/focal/base/set-home-permissions.sh /etc/my_init.d/set-home-permissions.sh

# Set up user so that we do not run as root
# See https://github.com/phusion/baseimage-docker/issues/186
# Disable forward logging
# Add script to set up permissions of home directory on myinit
# Run ldconfig so that /usr/local/lib is in the default search
# path for the dynamic linker.
RUN useradd -m -s /bin/bash -G sudo,docker_env,adm tfuser && \
    echo "tfuser:docker" | chpasswd && \
    echo "tfuser ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers && \
    echo "PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\W\[\033[00m\]\$ '" >> /home/tfuser/.bashrc && \
    echo "cat /home/tfuser/WELCOME" >> /home/tfuser/.bashrc && \
    rm /etc/my_init.d/10_syslog-ng.init && \
    chmod +x /etc/my_init.d/set-home-permissions.sh && \
    ldconfig

ENV INSTALL_PREFIX  /usr/local/terraferma
ENV INCLUDE         $INSTALL_PREFIX/include:$INCLUDE
ENV LD_LIBRARY_PATH $INSTALL_PREFIX/lib:$LD_LIBRARY_PATH
ENV LIBRARY_PATH    $INSTALL_PREFIX/lib:$LIBRARY_PATH
ENV PATH            $INSTALL_PREFIX/bin:$PATH
ENV PKG_CONFIG_PATH $INSTALL_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH
ENV PYTHONPATH      $INSTALL_PREFIX/lib/python3.8/site-packages:$PYTHONPATH
ENV MANPATH         $INSTALL_PREFIX/share/man:$MANPATH
 
RUN mkdir $INSTALL_PREFIX && \
    chgrp adm -R $INSTALL_PREFIX && \
    chmod g+w -R $INSTALL_PREFIX && \
    chmod a+w /tmp

USER tfuser
ENV TF_HOME /home/tfuser
RUN touch $TF_HOME/.sudo_as_admin_successful && \
    mkdir $TF_HOME/shared && \
    mkdir -p $INSTALL_PREFIX/lib/python3.8/site-packages
VOLUME /home/tfuser/shared

# Print a default welcome message on startup
COPY docker/focal/base/WELCOME $TF_HOME/WELCOME

WORKDIR $TF_HOME
USER root
ENTRYPOINT ["/sbin/my_init","--quiet","--","/sbin/setuser","tfuser","/bin/bash","-l","-c"]
CMD ["/bin/bash","-i"]

