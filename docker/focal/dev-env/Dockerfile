# Modified from the Dockerimage used by FEniCS:
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
# Builds a Docker image with the correct environment for compiling 
# TerraFERMA.  The built image can be found at:
#
#   https://hub.docker.com/r/terraferma/dev-env
#

# build from the base image
FROM terraferma/base:fenics-2019.1.0-focal
MAINTAINER terraferma <terraferma@lists.columbia.edu>

USER root
WORKDIR /tmp

# install non-python dependencies - splitting these up reduces the size of individual image layers
RUN apt-get update && \
    apt-get -y --with-new-pkgs \
        -o Dpkg::Options::="--force-confold" -o Dpkg::Options::="--force-confdef" upgrade && \
    apt-get -y install \
                       git \
                       gcc \
                       g++ \
                       gfortran \
                       cmake \
                       flex \
                       libopenblas-dev \
                       libeigen3-dev \
                       libboost-dev \
                       libboost-filesystem-dev \
                       libboost-program-options-dev \
                       libboost-iostreams-dev \
                       libboost-timer-dev \
                       libhdf5-dev \
                       libhdf5-openmpi-dev \
                       gmsh \
                       wget \
                       gir1.2-gtksource-3.0 \
                       gir1.2-gtk-3.0 \
                       adwaita-icon-theme-full \
                       2to3 \
                       tcl \
                       ack \
                       vim \
                       vtk7 \
                       evince \
                       python3 \
                       python3-ipdb \
                       python3-numpy \
                       python3-scipy \
                       python3-sympy \
                       python3-matplotlib \
                       python3-ply \
                       python3-pytest \
                       python3-setuptools \
                       python3-h5py \
                       python3-mpi4py \
                       python3-dev \
                       python3-future \
                       python3-gi \
                       python3-lxml \
                       cython3 \
                       python3-cairocffi \
                       python3-cairo \
                       python3-gi-cairo \
                       python3-pip \
                       python3-meshio \
                       python3-pybind11 \
                       python3-pkgconfig \
                       python3-vtk7 && \
    apt-get clean && \
    pip3 install ttictoc && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# set environment variables for OPENBLAS to make single threaded
ENV OPENBLAS_NUM_THREADS=1
ENV OPENBLAS_VERBOSE=0

USER tfuser

# download and install petsc from source
RUN git clone -b v3.14.6 -q https://gitlab.com/petsc/petsc.git && \
    cd petsc && \
    python3 ./configure --COPTFLAGS="-O3" \
                        --CXXOPTFLAGS="-O3" \
                        --FOPTFLAGS="-O3" \
                        --with-shared-libraries=1 \
                        --with-openblas=1 \
                        --with-debugging=0 \
                        --download-suitesparse=1 \
                        --download-hypre=1 \
                        --download-ml=1 \
                        --download-parmetis=1 \
                        --download-metis=1 \
                        --download-blacs=1 \
                        --download-scalapack=1 \
                        --download-mumps=1 \
                        --prefix=$INSTALL_PREFIX && \
    make all && \
    make install && \
    rm -rf /tmp/*
ENV PETSC_DIR $INSTALL_PREFIX

# download and install all the fenics components from source
# ufl
RUN git clone -b tferma-master-2019.1.0 -q https://bitbucket.org/tferma/ufl.git && \
    cd ufl && \
    python3 setup.py install --prefix=$INSTALL_PREFIX && \
    rm -rf /tmp/*
# fiat
RUN git clone -b tferma-master-2019.1.0 -q https://bitbucket.org/tferma/fiat.git && \
    cd fiat && \
    python3 setup.py install --prefix=$INSTALL_PREFIX && \
    rm -rf /tmp/*
# dijitso
RUN git clone -b tferma-master-2019.1.0 -q https://bitbucket.org/tferma/dijitso.git && \
    cd dijitso && \
    python3 setup.py install --prefix=$INSTALL_PREFIX && \
    rm -rf /tmp/*
# ffc
RUN git clone -b tferma-master-2019.1.0 -q https://bitbucket.org/tferma/ffc.git && \
    cd ffc && \
    python3 setup.py install --prefix=$INSTALL_PREFIX && \
    rm -rf /tmp/*
# dolfin
RUN git clone -b tferma-master-2019.1.0 -q https://bitbucket.org/tferma/dolfin.git && \
    cd dolfin && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DDOLFIN_ENABLE_CHOLMOD:BOOL=OFF \
          -DDOLFIN_ENABLE_UMFPACK:BOOL=OFF \
          -DDOLFIN_ENABLE_TRILINOS:BOOL=OFF \
          -DDOLFIN_ENABLE_SLEPC:BOOL=OFF  \
          -DMETIS_LIBRARY=${PETSC_DIR}/lib/libmetis.so \
          -DPARMETIS_INCLUDE_DIRS=${PETSC_DIR}/include \
          -DPARMETIS_LIBRARY=${PETSC_DIR}/lib/libparmetis.so \
          -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
          ../. && \
    make && \
    make install && \
    cd ../python && \
    python3 setup.py install --prefix=$INSTALL_PREFIX && \
    rm -rf /tmp/*

# download and install spud from source
RUN git clone -q https://github.com/FluidityProject/spud.git && \
    cd spud && \
    ./configure --prefix=$INSTALL_PREFIX && \
    make && \
    make install && \
    rm -rf /tmp/*
ENV SPUD_DIR $INSTALL_PREFIX

# supress some spurious gtk errors
ENV NO_AT_BRIDGE 1

# welcome at startup
COPY WELCOME $TF_HOME/WELCOME

WORKDIR $TF_HOME
USER root

