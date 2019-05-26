#ifndef _ATOM_H
#define _ATOM_H
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
      static std::vector< std::tuple< int, double, double> > NISTIsotopes(const char*); // return the atomic number, mass number and mass list from NIST database
      void setMass();
    public:
      Atom(const char *, double, double, double); // initialization with element symbol and Cartesian coodinates
      Atom & translate(const Eigen::Vector3d &); // translation
      Atom & transform(const Eigen::Matrix3d &); // maybe rotation presented by an rotation matrix, however we do not check it 
      double mass(); // get its mass
      int atomicNumber(); // get its atomic number
      int massNumber(); // get its mass number
      std::string element(); // get what kind of element it is
      Eigen::Vector3d & xyz(); // get the Cartesian coordinates
      Eigen::Vector3d & xyz(double, double, double); // set the Cartesian coordinates
      static std::tuple< int, double, double>  NISTAtom(const char*); // return the atomic number, mass number and mass from NIST database with most composition
      static std::tuple< int, double, double>  NISTAtom(const char*, int); // return the atomic number, mass number and mass from NIST database 
      std::string repr();

  };
}
#endif // _ATOM_H
