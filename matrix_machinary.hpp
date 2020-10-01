#pragma once

#include <dlib/matrix.h> 
#include <dlib/optimization.h>


namespace mm {

using namespace dlib; 


using matrix  = dlib::matrix<double, 0, 0>;
using cvector = dlib::matrix<double, 0, 1>; 
using rvector = dlib::matrix<double, 1, 0>;

matrix join_column_units(const matrix& x) {

  matrix t(x.nr(), x.nc() + 1);

  t = subm(x, range(0,x.nr()-1), range(0,x.nc()));
  set_colm(t, x.nc()) = 1.0;

  return t;
}


} // mm

