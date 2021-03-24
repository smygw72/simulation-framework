// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "transpose_blocks.h"

#include <cassert>

template <typename T>
IGL_INLINE void igl::transpose_blocks(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & A,
  const size_t k,
  const size_t dim,
  Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & B)
{
  // Eigen matrices must be 2d so dim must be only 1 or 2
  assert(dim == 1 || dim == 2);
  // Output is not allowed to be input
  assert(&A != &B);


  // block height, width, and number of blocks
  int m,n;
  if(dim == 1)
  {
    m = A.rows()/k;
    n = A.cols();
  }else// dim == 2
  {
    m = A.rows();
    n = A.cols()/k;
  }

  // resize output
  if(dim == 1)
  {
    B.resize(n*k,m);
  }else//dim ==2
  {
    B.resize(n,m*k);
  }

  // loop over blocks
  for(int b = 0;b<(int)k;b++)
  {
    if(dim == 1)
    {
      B.block(b*n,0,n,m) = A.block(b*m,0,m,n).transpose();
    }else//dim ==2
    {
      B.block(0,b*m,n,m) = A.block(0,b*n,m,n).transpose();
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template void igl::transpose_blocks<double>(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, size_t, size_t, Eigen::Matrix<double, -1, -1, 0, -1, -1>&);
#endif