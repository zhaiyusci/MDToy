#ifndef _MDTOY_ATOM_H
#define _MDTOY_ATOM_H
#include<string>
#include<Eigen/Eigen>
#include<tuple>

namespace MDToy{
  class Atom{
    private:
      std::string elem_; // type of element of atom, case sensitive
      double mass_; // atomic mass, in Da
      int atomicNumber_; // atomic number, for D it is 1
      int massNumber_; // mass number, for D it is 2
      Eigen::Vector3d xyz_; // position in Cartesian coordinates, in Angstrom
      static std::vector< std::tuple< int, double, double> > NISTIsotopes(const std::string&); // return the atomic number, mass number and mass list from NIST database
      void setMass();
    public:
      Atom(const std::string&, double, double, double); // initialization with element symbol and Cartesian coodinates
      Atom & translate(const Eigen::Vector3d &); // translation
      Atom & transform(const Eigen::Matrix3d &); // maybe rotation presented by an rotation matrix, however we do not check it 
      double mass() const; // get its mass
      int atomicNumber() const; // get its atomic number
      int massNumber() const; // get its mass number
      std::string element() const; // get what kind of element it is
      const Eigen::Vector3d & xyz() const; // get the Cartesian coordinates
      double xyz(int i) const; // get the Cartesian coordinates
      const Eigen::Vector3d & xyz(double, double, double); // set the Cartesian coordinates
      static std::tuple< int, double, double>  NISTAtom(const std::string&); // return the atomic number, mass number and mass from NIST database with most composition
      static std::tuple< int, double, double>  NISTAtom(const std::string&, int); // return the atomic number, mass number and mass from NIST database 
      std::string repr();

  };

  class WeightOfAtom{
    public:
      virtual double operator()(const Atom&)const=0;
      virtual std::string what() const=0;
  };
}
#endif // _MDTOY_ATOM_H
