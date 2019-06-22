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
      weight_.push_back(1.0);
    }

    phi_=&massweighted;

    updateCenter();

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
      throw;
    }

    updateCenter();

    return;
  }

  void Molecule::updateCenter(const WeightOfAtom &f, std::vector<double> weight){
    phi_=&f;
    weight_=weight;
    double w;
    center_=Eigen::Vector3d::Zero();
    int ii=0;
    for(auto && i : list_){
      w+=f(i)*weight[ii];
      center_+=i.xyz()*f(i)*weight[ii];
    }
    center_/=w;
  }

  void Molecule::updateCenter(){
    updateCenter(*phi_, weight_ );
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

  Molecule & Molecule::translate (const Eigen::Vector3d & r){
    for( auto && i : list_){
      i.translate(r);
    }
    center_+=r;
    return *this;
  }

  Molecule & Molecule::transform(const Eigen::Matrix3d & r){
    for(auto && i:list_){
      i.transform(r);
    }
    auto&& q=r*center_;
    center_=q;
    return *this;
  }

  Molecule & Molecule::setCenter0(){
    translate(-center_);
    // center_=Eigen::Vector3d::Zero();
    return *this;
  }

  Eigen::MatrixXd Molecule::xyz()const{
    Eigen::MatrixXd a(3,list_.size());

    int c=0; // column number
    for(auto && i: list_){
      a.col(c)=i.xyz();
      c++;
    }
    return a;
  }

  const std::vector< Atom >& Molecule::atoms() const{
    return list_;
  }

  const Atom & Molecule::atoms(int i) const{
    return list_[i];
  }

  const std::vector<double>& Molecule::weight() const{
    return weight_;
  }

  double Molecule::weight(int i) const{
    return weight_[i];
  }

  const WeightOfAtom &Molecule::phi() const{
    return *phi_;
  }

  const Eigen::Vector3d & Molecule::center(){
    return center_;
  }

}
