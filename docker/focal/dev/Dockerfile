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
#   https://github.com/orgs/TerraFERMA/packages/container/package/dev
#

# build using the TerraFERMA development environment
ARG TAG=fenics-2019.1.0-focal
FROM ghcr.io/terraferma/dev-env:${TAG}
MAINTAINER terraferma <terraferma@lists.columbia.edu>

ARG BRANCH=master-2019.1.0
ARG REPO=https://github.com/TerraFERMA/TerraFERMA.git
ARG DEBUG=0

USER tfuser
WORKDIR /tmp

# install terraferma
RUN git clone -b $BRANCH -q $REPO && \
    cd TerraFERMA && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
          -DCMAKE_BUILD_TYPE=`[[ $DEBUG -eq 1 ]] && echo "Debug" || echo "Release"` \
          ../. && \
    make && \
    make install && \
    rm -rf /tmp/*

ENV TF_CMAKE_PATH       $INSTALL_PREFIX/share/terraferma/cpp
ENV CMAKE_MODULE_PATH   $INSTALL_PREFIX/share/buckettools/cmake/modules:$CMAKE_MODULE_PATH
ENV DIAMOND_CONFIG_PATH $INSTALL_PREFIX/share/buckettools/diamond:$DIAMOND_CONFIG_PATH

# link the installed tests and tutorials directories to the home directory
RUN ln -s $INSTALL_PREFIX/share/terraferma/tests $TF_HOME/tests && \
    ln -s $INSTALL_PREFIX/share/terraferma/tutorials $TF_HOME/tutorials
COPY docker/focal/dev/WELCOME $TF_HOME/WELCOME

WORKDIR $TF_HOME
USER root

