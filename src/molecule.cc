#include"molecule.h"
#include<fstream>
#include<iostream>
#include<string>
#include<sstream>

namespace MDToy{
  Molecule::Molecule(std::istream& ins){
    std::string line;
    std::getline(ins, line);
    std::stringstream ss(line);
    int n;
    double x,y,z;
    ss>> n;
    std::getline(ins, line);
    for(int i =0; i!=n; ++i){
      std::getline(ins,line);
      std::string elem;
      std::stringstream(line) >> elem >> x >> y >> z;
      list_.push_back(Atom(elem, x, y, z));
    }
    m_=0.0;
    for(auto && i: list_){
      m_+=i.mass();
    }
    updateCOM();
    return;
  }

  void Molecule::updateCoordinates(std::istream& ins){
    try{
      if(!ins) throw std::runtime_error("End of file.");
      std::string line;
      std::getline(ins, line);
      if(!ins) throw std::runtime_error("End of file.");
      std::stringstream ss(line);
      int n;
      double x,y,z;
      ss>> n;
      std::getline(ins, line);
      if(n==list_.size()){
        for(int i =0; i!=n; ++i){
          std::getline(ins,line);
          std::string elem;
          std::stringstream(line) >> elem >> x >> y >> z;
          if(elem==list_[i].element()){
            list_[i].xyz(x,y,z);
          }
        }
      }
    }
    catch(std::runtime_error err){
      std::cerr <<"Line " << __LINE__ << ": " << err.what() << std::endl;
      // std::abort();
      // std::terminate();
      throw;
      // return;
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

    auto r=ss.str();
    r.erase(r.end()-1);
    return r;
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
