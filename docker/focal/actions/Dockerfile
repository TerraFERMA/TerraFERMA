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
# Builds a Docker image with a full installation of TerraFERMA.  
# The built image can be found at:
#
#   https://hub.docker.com/r/terraferma/dev
#

# build using the TerraFERMA development environment
FROM terraferma/dev-env:fenics-2019.1.0-focal
MAINTAINER terraferma <terraferma@lists.columbia.edu>

USER root

ENV TF_PREFIX /usr/local/terraferma 
RUN mkdir $TF_PREFIX && \
    chgrp adm -R $TF_PREFIX && \
    chmod g+w -R $TF_PREFIX

COPY . $TF_HOME
RUN chown -R tfuser $TF_HOME

USER tfuser
WORKDIR $TF_HOME

RUN export PYTHONPATH=$FENICS_PREFIX/lib/python3.8/site-packages && \
    export PETSC_DIR=$PETSC_PREFIX && \
    export PATH=$FENICS_PREFIX/bin:$PATH && \
    export LD_LIBRARY_PATH=$FENICS_PREFIX/lib:$LD_LIBRARY_PATH && \
    export PKG_CONFIG_PATH=$FENICS_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH && \
    export SPUD_DIR=$SPUD_PREFIX && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=${TF_PREFIX} \
          -DCMAKE_BUILD_TYPE=Release  \
          ../. && \
    make && \
    make install

# install an environment module for terraferma
COPY docker/focal/actions/terraferma.configmodule /usr/local/share/environment-modules/terraferma

COPY docker/focal/actions/bash.load-modules $TF_HOME/.bash.load-modules
RUN cat $TF_HOME/.bash.load-modules >> $TF_HOME/.profile
COPY docker/focal/actions/WELCOME $TF_HOME/WELCOME

USER root

