#include"atom.h"
#include<fmt/format.h>
#include<boost/regex.hpp>
#include<iostream>
#include<sstream>
#include<fstream>
#include<algorithm>
#include<chrono>

namespace MDToy{
  Atom::Atom(const std::string& e, double x, double y, double z){
    elem_=e;
    xyz_ << x, y, z;
    setMass();
  }

  void Atom::setMass(){
    auto [atomicNumber_, massNumber_, mass_]=NISTAtom(elem_.c_str());
    return;
  }

  std::vector< std::tuple< int, double, double> > Atom::NISTIsotopes(const std::string& e){
    int atomicn;
    int massn;
    double mass;
    static std::string nistdata;
    // std::cout << nistdata.size() << std::endl;
    std::string re=e;
    re="_[^_]*? "+re+" [^_]*?_";

    // read in NIST data
    if (nistdata.size() == 0){
      std::cout << "Loading data... "<< std::endl;
      std::cout << "Please cite:" <<std::endl;
      std::cout << "    Wang , Audi et al., Chin. Phys. C 36, 1603 (2012) doi:10.1088/1674-1137/36/12/003"  <<std::endl;
      std::ifstream in("../data/atomicweight_nist2015.shortest.txt");
      nistdata=static_cast<std::stringstream const &>(std::stringstream() << in.rdbuf()).str();
    }

    std::vector< std::tuple<int, double, double> > isotopes;
    boost::regex element(re.c_str(), boost::regex_constants::ECMAScript);
    boost::smatch ss;
    if(boost::regex_search(nistdata, ss, element) ){
      // std::cout<< "Found: " << ss.size() << " " << ss.str(0) << std::endl;
      std::stringstream elemtext(ss.str(0));
      std::string line;
      while(getline(elemtext, line)){
        // std::cout << line << std::endl;
          boost::smatch iso;
        if(boost::regex_match(line, iso, 
              boost::regex(".{8}([[:digit:]][[:digit:][:space:]]{2})  ([[:digit:].]+)[^[:space:]]* {1,10}([[:digit:].]+)[^[:space:]]*.*", 
                boost::regex_constants::ECMAScript))){
          // std::cout << iso[0] << "   " << iso[1] << "   " << iso[2] << "   " << iso[3] << std::endl;
          int mn;
          std::stringstream(iso[1])>> mn;
          double m;
          std::stringstream(iso[2])>> m;
          double comp;
          std::stringstream(iso[3])>> comp;
          isotopes.push_back(std::make_tuple(mn, m, comp));
        }else if(boost::regex_match(line, iso,
              boost::regex(".{8}([[:digit:]][[:digit:][:space:]]{2})  ([[:digit:].]+)[^[:space:]]*.*", 
                boost::regex_constants::ECMAScript))){
          // std::cout << iso[0] << "   " << iso[1] << "   " << iso[2]  << std::endl;
          int mn;
          std::stringstream(iso[1])>> mn;
          double m;
          std::stringstream(iso[2])>> m;
          double comp=0.0;
          isotopes.push_back(std::make_tuple(mn, m, comp));
        }
      }
    }else{
      std::cerr<< "Error: " << " Element " << e << " not found!" << std::endl;
      isotopes.push_back( std::make_tuple(0,0.0,0.0));
    }
    return isotopes;
  }

  std::tuple<int, double, double> Atom::NISTAtom(const std::string& e){

      std::vector< std::tuple<int, double, double> > isotopes(NISTIsotopes(e));
      std::sort(isotopes.begin(), isotopes.end(), 
          [](std::tuple<int, double, double> a , std::tuple<int, double, double> b) -> bool{ return std::get<2>(a) > std::get<2>(b);});
        // std::cout << std::get<0>(isotopes[0]) << "   " << std::get<1>(isotopes[0]) << "   " <<std::get<2>(isotopes[0]) << std::endl;
      return isotopes[0];
  }
  std::tuple<int, double, double> Atom::NISTAtom(const std::string& e, int mn){

      std::vector< std::tuple<int, double, double> > isotopes(NISTIsotopes(e));
      for(auto i : isotopes){
        if(std::get<0>(i) == mn) return i;
      }
      std::cerr<< "Error: " << " Isotope " << e << "-" << mn << " not found!" << std::endl;
      return isotopes[0];
  }

  double Atom::mass() {
    return mass_;
  }

  int Atom::massNumber(){
    return massNumber_;
  }

  int Atom::atomicNumber(){
    return atomicNumber_;
  }

  Eigen::Vector3d & Atom::xyz(){
    return xyz_;
  }

  Eigen::Vector3d & Atom::xyz(double x, double y, double z){
    xyz_<< x, y, z;
    return xyz_;
  }

  std::string Atom::element(){
    return elem_;
  }

  std::string Atom::repr(){
    return fmt::format("{:3s} {:18.8f} {:18.8f} {:18.8f}", elem_, xyz_[0], xyz_[1], xyz_[2]);
  }

  Atom & Atom::translate(const Eigen::Vector3d & r){
    xyz_+=r;
    return *this;
  }

  Atom & Atom::transform(const Eigen::Matrix3d & r){
    auto&& q=r*xyz_;
    xyz_=q;
    return *this;
  }

}
