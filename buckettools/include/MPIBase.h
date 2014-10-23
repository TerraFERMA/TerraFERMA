// Copyright (C) 2013 Columbia University in the City of New York and others.
//
// Please see the AUTHORS file in the main source directory for a full list
// of contributors.
//
// This file is part of TerraFERMA.
//
// TerraFERMA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TerraFERMA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.


#ifndef __MPI_BASE_H
#define __MPI_BASE_H

#include <dolfin.h>
#ifdef HAS_MPI
#include <mpi.h>

namespace buckettools
{

  void mpi_error(const int &mpierror,                                // translate mpi error codes into TF errors
                 const std::string &filename, 
                 const int &line,
                 const int &accept=MPI_SUCCESS);

  #define mpi_err(mpierr) do {mpi_error(mpierr, __FILE__, __LINE__);} while(0)

  #define mpi_err_accept(mpierr, accept) do {mpi_error(mpierr, __FILE__, __LINE__, accept);} while(0)

}
#endif

#endif

