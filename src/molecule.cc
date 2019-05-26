#include"molecule.h"
#include<fstream>
#include<iostream>
#include<string>
#include<sstream>

namespace MDToy{
  Molecule::Molecule(const char *filename){
    std::ifstream file(filename);
    std::string line;
    getline(file, line);
    std::stringstream ss(line);
    int n;
    double x,y,z;
    ss>> n;
    getline(file, line);
    for(int i =0; i!=n; ++i){
      getline(file,line);
      std::string elem;
      std::stringstream(line) >> elem >> x >> y >> z;
      list_.push_back(Atom(elem.c_str(), x, y, z));
    }
    file.close();
    m_=0.0;
    for(auto && i: list_){
      m_+=i.mass();
    }
    updateCOM();
    return;
  }

  void Molecule::updateCOM(){
    // set the mass and c o m
    cm_=Eigen::Vector3d::Zero();
    for(auto && i: list_){
      cm_+=i.mass()*i.xyz();
    }
    cm_/=m_;
  }

  std::string Molecule::repr(){
    std::stringstream ss;
    ss<< list_.size() <<std::endl <<std::endl;
    for(auto && i : list_){
      ss<< i.repr() <<std::endl;
    }
    return ss.str();
  }

  Eigen::Vector3d Molecule::centerOfMass() {
    return cm_;
  }

  Molecule & Molecule::translate (const Eigen::Vector3d & r){
    for( auto && i : list_){
      i.translate(r);
    }
    cm_+=r;
    return *this;
  }

  Molecule & Molecule::transform(const Eigen::Matrix3d & r){
    for(auto && i:list_){
      i.transform(r);
    }
    auto&& q=r*cm_;
    cm_=q;
    return *this;
  }

  Molecule & Molecule::setCenterOfMass0(){
    translate(cm_);
    return *this;
  }

}
