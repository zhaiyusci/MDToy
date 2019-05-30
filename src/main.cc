#include"molecule.h"
#include"atom.h"
#include<iostream>
#include<fstream>
int main(){

  std::ifstream file;
  file.open("ch4.xyz");
  MDToy::Molecule ch4(file);
  std::cout << ch4.repr()<<std::endl;
  std::cout << ch4.centerOfMass() <<std::endl;
  file.close();

  Eigen::Vector3d r;
  r << 1.0, 2.0, 3.0;
  ch4.translate(r);
  std::cout << ch4.repr()<<std::endl;
  std::cout << ch4.centerOfMass() <<std::endl;

  // file.exceptions(std::istream::failbit|std::istream::badbit|std::istream::eofbit);
  file.open("ch4.md.xyz");
  for(int i=0; i!=50; ++i){
    ch4.updateCoordinates(file);
    std::cout << ch4.repr()<<std::endl;
    std::cout << ch4.centerOfMass() <<std::endl;
  }
  file.close();


  MDToy::Atom h("H", 0.0, 0.0, 0.0);
  std::cout << h.repr()<<std::endl;
  h.translate(r);
  std::cout << h.repr()<<std::endl;


  return 0;
}

