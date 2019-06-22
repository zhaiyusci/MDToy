#ifndef _MDTOY_MOLECULE_H
#define _MDTOY_MOLECULE_H
#include"atom.h"
#include<Eigen/Eigen>

namespace MDToy{
  class Molecule{
    private:
      std::vector< Atom > list_; // atom list
      Eigen::MatrixXd mode_; // normal or local mode
      Eigen::Vector3d center_; // center by some weight ...
      // double m_; // total mass
      const WeightOfAtom (*phi_) ;
      std::vector<double> weight_;
    public:
      const Eigen::Vector3d & center() ; // get center of mass
      Molecule(std::istream &); // initialize the class from xyz formate file
      void updateCoordinates(std::istream &); // initialize the class from xyz formate file
      Molecule & setCenter0(); 
      Molecule & updateCenter(const WeightOfAtom &, std::vector<double>); // calc the center using the above weights
      Molecule & updateCenter(); // calc using the old weight
      Molecule & transform(const Eigen::Matrix3d &);
      Molecule & translate(const Eigen::Vector3d &);
      void vibrateAlongMode(double );
      std::string repr();
      Eigen::MatrixXd xyz()const; // return a 3*N-sized matrix for some molecule-based calculation
      const std::vector< Atom >& atoms() const;
      const Atom & atoms(int i) const;
      const std::vector<double> &weight() const ;
      double weight(int)const;
      const WeightOfAtom& phi() const;

  };

  class MassWeighted: public WeightOfAtom{
    public:
      double operator() (const Atom& a) const{
        return a.mass();
      }
      std::string what() const{
        return "MassWeighted";
      }
  };
  const MassWeighted massweighted;
}
#endif // _MDTOY_MOLECULE_H
