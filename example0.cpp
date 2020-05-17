#include <iostream>

// #include "algo/generators/mackey_glass_generator.hpp"

#include "ssa.hpp"

using namespace nnn;

void print(const matrix& mat) {

  for( auto vv : mat)
    std::cout << vv << " ";
  std::cout << "\n";               
}

int main() {

  std::cout << std::fixed << std::setprecision(2);

  matrix sq(1, 15); 

  for(int i=0; i<15; i++)
    sq(i) = 1*std::sin(i*0.1) * std::exp(i*0.01);
                          // 8
  nnn::ssa_model ssa(sq, 15, 5); 

  std::cout << "l: " << ssa.traectory().nr() << " k: "  << ssa.traectory().nc() << "\n";

  const matrix& rr = ssa.sequence();
  print(rr);

//  matrix sm(ssa.traectory());
//  sm *= trans(sm);
//  std::cout << sm << "\n"; 

//  std::cout << ssa.traectory() << "\n"; 
//  std::cout << ssa.u_basis() << "\n";
//  std::cout << ssa.v_basis() << "\n";
//  std::cout << ssa.w_matrix() << "\n";

  const matrix& ld = ssa.lambdas();
  std::cout << std::fixed << std::setprecision(6); 
  std::cout << ld;

  std::cout << std::fixed << std::setprecision(2); 

  for(size_t ix=0; ix<ld.nr(); ix++) {
  
    matrix mt = ssa.component(ix);
    print(mt);
  }

  matrix comp = ssa.composition(2);  // first 2, 
  print(comp);

  comp = ssa.composition(2, ld.nr()-2);
  print(comp);

        
  return 0;
}

