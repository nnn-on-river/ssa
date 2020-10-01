#pragma once

#include <cassert>

#include "..\matrix_machinary.hpp"

// for forecasting,  q  signal noise ratio
// q -> 0.0  L -> N/2;
// q -> ***  L = o(N)  L -> sqrt(N);

namespace nnn {

class ssa { 

public: 

  ssa(const mm::cvector& data, size_t l) { 

    size_t n = data.size(); 
    l_ = l; k_ = n-l+1;
 
    assert( l > 1 && n > l); 

    build_traectory(data);
    build_model();
  }

  template<typename Serie>
  ssa(const Serie& data, size_t l) { 

    size_t n = data.size(); 
    l_ = l; k_ = n-l+1;
 
    assert( l > 1 && n > l); 

    build_traectory(data);
    build_model();
  } 
 

  size_t l_parameter() const { return l_; }
  size_t k_parameter() const { return k_; }
  size_t n_parameter() const { return l_+k_-1; }

  const mm::matrix& traectory() const { return traectory_; }

  const mm::matrix&  u_basis() const { return u_base_; }
  const mm::matrix&  v_basis() const { return v_base_; }
  const mm::cvector& lambdas() const { return lambdas_; }
  
  const mm::cvector& rebuild(size_t l);  

  mm::cvector lr_formula(size_t c_count = 0) const { 
    assert(c_count < u_base_.nc());

    size_t count = c_count == 0 ? u_base_.nc()-1 : c_count;
    return build_lrf(u_base_, count);
  }

  mm::cvector compose(size_t total) const { // the sum of the first components
    assert( lambdas_.nr() >= total); 
    return compose(0, total);
  }

  mm::cvector compose(size_t index, size_t amount) const {

    assert( lambdas_.nr() >= index + amount); 

    mm::matrix rr;      

    for(size_t i = index; i < (index + amount); i++) {

      mm::cvector ui = mm::colm(u_base_, i);
      mm::cvector vi = mm::colm(v_base_, i);

      rr +=  ui * mm::trans(vi);
    }

    rr *= traectory_; 
    return do_restore(rr);

  }

  mm::cvector component(size_t i) const { 
    assert(i < lambdas_.size());

    mm::cvector ui = mm::colm(u_base_, i);
    mm::cvector vi = mm::colm(v_base_, i);

    mm::matrix t =  ui * mm::trans(vi);

    t *= traectory_; 

    return do_restore(t);
  }


  mm::cvector correlation_weights_vector() const {
    return build_weights_vector();
  }

  double w_correlation(const mm::cvector& x, const mm::cvector& z, const mm::cvector& w) const;

private: 

  size_t l_;
  size_t k_;

  mm::matrix  traectory_;

  mm::matrix  u_base_;
  mm::matrix  v_base_;
  mm::cvector lambdas_;
//  mm::matrix  w_;

private:
  void build_model();
  mm::cvector build_lrf(const mm::matrix& u, size_t cc) const;

  template<typename Serie>
  void build_traectory(const Serie& d);

  void do_ranking();

  mm::cvector do_restore( mm::matrix& mx) const;

  mm::cvector build_weights_vector() const; 

};

void ssa::build_model() {

  mm::matrix w_;
  mm::matrix s_matrix = traectory_ * mm::trans(traectory_);
  mm::svd(s_matrix, u_base_, w_, v_base_);

  lambdas_ = mm::diag(w_);

  do_ranking();
//std::cout << "w:\n" << w_;

}  

mm::cvector ssa::do_restore(mm::matrix& mx) const {  

  size_t LL = l_ < k_ ? l_ : k_;
  size_t KK = k_ > l_ ? k_ : l_;
  size_t NN = l_ + k_ - 1;

  mm::cvector xx(NN);

  if(l_ > k_)
    mx = mm::trans(mx);

  double x; size_t j;
  for( size_t i = 0; i < LL-1 ; i++)  {
    x = 0.0;
    for( j = 0; j < i+1; j++)
      x += mx(i-j, j);

    xx(i) = x / (j);
  }

  for( size_t i = LL-1; i < KK ; i++)  {
    x = 0.0;
    for( j = 0; j < LL ; j++)  
      x += mx(j, i-j);
    xx(i) = x / (j);
  }

  for( size_t i = KK; i < NN ; i++) {
    x = 0.0;
    for( j = i-KK+1; j < NN-KK+1 ; j++)
      x += mx(j, i-j);
    xx(i) = x / (NN - i);
  }

  return xx;
}

template<typename Serie>
void ssa::build_traectory(const Serie& serie) {    

  auto bb = serie.begin();
  traectory_.set_size(l_, k_);

  for( size_t i = 0; i < l_; i++) {

    for(size_t j = 0; j < k_; j++) 
      traectory_(i, j) = *bb++;
    
    bb = serie.begin(); 
    bb += (i+1);
  }  
}

void ssa::do_ranking() { 

  for(size_t i = 1; i < lambdas_.size(); i++) {
    for(size_t j = 0; j < i; j++)  {

      if(lambdas_(i) > lambdas_(j)) { 

        double t = lambdas_(j);
        lambdas_(j) = lambdas_(i);
        lambdas_(i) = t;

        mm::cvector x = mm::colm(u_base_, i);
        mm::set_colm(u_base_, i) = mm::colm(u_base_, j);
        mm::set_colm(u_base_, j) = x;

        x = mm::colm(v_base_, i);
        mm::set_colm(v_base_, i) = mm::colm(v_base_, j);
        mm::set_colm(v_base_, j) = x;

      }
    }
  }
} 

const mm::cvector& ssa::rebuild(size_t l ) { 

  size_t n = l_ + k_ - 1;
  assert( l > 1 && n > l);

  mm::cvector data(n);

  long k; 
  for(k=0; k<k_; k++) 
    data(k) = traectory_(0, k);

  for(long i=1; i<l_; i++) 
    data(k++) = traectory_(i, k_-1); 
  
  l_ = l; 
  k_ = n-l+1; 

  build_traectory(data);
  build_model();

  return lambdas_;
}  

mm::cvector ssa::build_lrf(const mm::matrix& u, size_t c_count) const { 

  double p; 
  double vp = 0.0;

  mm::cvector lrf(u.nr()-1, 1);
  lrf = 0.0;

  for(long i = 0; i < c_count; i++) { 

    mm::cvector x = mm::colm(u, i);
    mm::cvector z = mm::subm(x, dlib::rectangle(0, 0, 0, x.size()-2));

    p = x(x.size() - 1);
    vp += p*p; 

    z *= p;
    lrf += z;
  }

  double f = 1.0 / (1.0 - vp);
  lrf *= f;

  return lrf;
} 

mm::cvector ssa::build_weights_vector() const { 

  size_t n = l_ + k_ - 1;
  mm::cvector w(n);

  for(long i=0; i < n; i++) {
  
    if( i >= 0 && i < k_)
      w(i) = i+1;
    else if( i > k_ && i <= l_)
      w(i) = k_;
    else if( i >= l_ && i < l_+k_)
      w(i) = n - i;
    else
      w(i) = 0.0;  

  }
  
  return w;
}
 
double ssa::w_correlation(const mm::cvector& x, const mm::cvector& z, const mm::cvector& w) const { 

  double no = 0.0;
  double d0 = 0.0;
  double d1 = 0.0;

  for(long i=0; i < x.size(); i++) { 

    double u = w(i);

    no += u*x(i)*z(i);
    d0 += u*x(i)*x(i);
    d1 += u*z(i)*z(i);
  }       

  double r = no / std::sqrt(d0*d1);
  return r;
}

} // vi
