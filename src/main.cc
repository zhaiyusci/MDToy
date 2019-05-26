#include"molecule.h"
#include"atom.h"
#include<iostream>
int main(){
  MDToy::Molecule ch4("ch4.xyz");
  std::cout << ch4.repr()<<std::endl;
  std::cout << ch4.centerOfMass() <<std::endl;

  Eigen::Vector3d r;
  r << 1.0, 2.0, 3.0;
  ch4.translate(r);
  std::cout << ch4.repr()<<std::endl;
  std::cout << ch4.centerOfMass() <<std::endl;


  MDToy::Atom h("H", 0.0, 0.0, 0.0);
  std::cout << h.repr()<<std::endl;
  h.translate(r);
  std::cout << h.repr()<<std::endl;

  
  return 0;
}

