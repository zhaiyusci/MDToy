#include"molecule.h"
#include"atom.h"
#include"kabsch.h"
#include<iostream>
#include<fstream>
#include<Eigen/Geometry>
int main(){

  std::ifstream file;
  file.open("ch4.xyz");
  MDToy::Molecule ch4(file);
  std::cout << "origin"<<std::endl;
  std::cout << ch4.repr()<<std::endl;
  std::cout << ch4.center() <<std::endl;
  file.close();

  // return 0;

  ch4.setCenter0();
  std::cout << "reset the C.O.M"<<std::endl;
  std::cout << ch4.repr()<<std::endl;
  std::cout << ch4.center() <<std::endl;
  std::ofstream output;

  MDToy::Molecule ch4ref=ch4;
  output.open("ref.xyz");
  output<< ch4ref.repr()<<std::endl;
  output.close();

  // file.exceptions(std::istream::failbit|std::istream::badbit|std::istream::eofbit);
  file.open("ch4.md.xyz");
  output.open("out.xyz");

  std::ofstream output2;
  output2.open("out2.xyz");
  for(int i=0; i!=22; ++i){
    ch4.updateCoordinates(file);
    ch4.setCenter0();
    Eigen::Matrix3d m;
    m = Eigen::AngleAxisd(1.0, Eigen::Vector3d::UnitZ())
      *Eigen::AngleAxisd(1.0, Eigen::Vector3d::UnitY())
      *Eigen::AngleAxisd(1.0, Eigen::Vector3d::UnitZ());
    ch4.transform(m);
    output2<< ch4.repr() << std::endl;

    Eigen::Matrix3d U= MDToy::kabsch(ch4, ch4ref);
    std::cout << U << std::endl;
    output<< ch4.transform(U).repr() << std::endl;
  }
  file.close();
  output.close();
  output2.close();



  return 0;
}

