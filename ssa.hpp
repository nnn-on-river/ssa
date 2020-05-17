#pragma once

#include "dlib\matrix\matrix_la.h"

#include <cassert>

namespace nnn {

template <class Type>
using matrix_type = dlib::matrix<Type>;

using matrix = matrix_type<double>;


class ssa_model { 

  public: 

    ssa_model(const matrix& data, size_t n, size_t l) : l_(l), k_(n - l + 1) {

      assert( l > 1 && n > l); 

      serie_.set_size(1, n); 

      for(size_t i = 0; i < n; i++)
        serie_(i) = data(i);

      build_model();
    } 

    const matrix& sequence()  const { return serie_;  }
    const matrix& traectory() const { return traectory_; }

    const matrix& u_basis()  const { return u_base_; }
    const matrix& v_basis()  const { return v_base_; }
    const matrix& w_matrix() const { return w_; }

    const matrix& lambdas() const { return lambdas_; }

    matrix composition(size_t a_mount) const { // the sum of the first components

      assert( lambdas_.nr() >= a_mount); 
  
      return composition(0, a_mount);
    }

    matrix composition(size_t index, size_t a_mount) const {

      assert( lambdas_.nr() >= index + a_mount); 

      matrix rr;      

      for(size_t ix=index; ix<(index + a_mount); ix++) {

        matrix ui = colm(u_base_, ix);
        matrix vi = colm(v_base_, ix);

        rr +=  ui * trans(vi);
      }

      rr *= traectory_; 
      return do_restore(rr);

    }

    matrix component(size_t ix) const { 

      assert(ix < lambdas_.nr());

      matrix ui = colm(u_base_, ix);
      matrix vi = colm(v_base_, ix);

      matrix temp =  ui * trans(vi);
      temp *= traectory_; 

      return do_restore(temp);
    }


  private: 

    size_t l_;
    size_t k_;

    matrix serie_;
    matrix traectory_;

    matrix u_base_;
    matrix v_base_;
    matrix w_;

    matrix lambdas_;


  private:
    void build_model();

    void do_traectory();
    size_t do_ranking();
    matrix do_restore( matrix& mm) const;
};

void ssa_model::build_model() {

  do_traectory();

  matrix s_matrix = traectory_ * trans( traectory_);
  svd( s_matrix, u_base_, w_, v_base_);

  lambdas_ = diag(w_);

  do_ranking();
}  


matrix ssa_model::do_restore( matrix& mm) const {  

  size_t LL = l_ < k_ ? l_ : k_;
  size_t KK = k_ > l_ ? k_ : l_;
  size_t NN = l_ + k_ - 1;

  matrix xx;
  xx.set_size(1, NN);

  if(l_ > k_)
    mm = trans(mm);

  double x; size_t j;
  for( size_t i = 0; i < LL-1 ; i++)  {
    x = 0.0;
    for( j = 0; j < i+1; j++)
      x += mm(i-j, j);

    xx(i) = x / (j);
  }

  for( size_t i = LL-1; i < KK ; i++)  {
    x = 0.0;
    for( j = 0; j < LL ; j++)  
      x += mm(j, i-j);
    xx(i) = x / (j);
  }

  for( size_t i = KK; i < NN ; i++) {
    x = 0.0;
    for( j = i-KK+1; j < NN-KK+1 ; j++)
      x += mm(j, i-j);
    xx(i) = x / (NN - i);
  }

  return xx;
}

void ssa_model::do_traectory() {    
  
  traectory_.set_size(l_, k_);
  for( size_t i = 0; i < l_; i++)
    for(size_t j = 0; j < k_; j++) {
      double vv = serie_(i+j);
      traectory_(i, j) = vv;
    }
}

size_t ssa_model::do_ranking() { 

  matrix x, y;
  size_t positive = 0;

  if( lambdas_(0) > 0.0)
    positive = 1; 

  for(size_t i = 1; i < lambdas_.size(); i++) {

    double vv = lambdas_(i);
    if(vv > 0.0 )
      positive++;

    for(size_t j = 0; j < i; j++)  {

      if(vv > lambdas_(j)) { 

        double t = lambdas_(j);

        lambdas_(j) = lambdas_(i);
        lambdas_(i) = t;

        x = colm(u_base_, i);
        y = colm(u_base_, j);

        set_colm(u_base_, i) = y;
        set_colm(u_base_, j) = x;

        x = colm(v_base_, i);
        y = colm(v_base_, j);

        set_colm(v_base_, i) = y;
        set_colm(v_base_, j) = x;

      }
    }
  }

  //std::cout << "positive: " << positive << "\n";
  return positive;
}

}
