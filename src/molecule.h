#ifndef _MOLECULE_H
#define _MOLECULE_H
#include"atom.h"
#include<Eigen/Eigen>

namespace MDToy{
  class Molecule{
    private:
      std::vector< Atom > list_; // atom list
      Eigen::MatrixXd mode_; // normal or local mode
      Eigen::Vector3d cm_; // center of mass
      double m_; // total mass
      void updateCOM(); // update center of mass
    public:
      Eigen::Vector3d centerOfMass() ; // get center of mass
      Molecule(std::istream &); // initialize the class from xyz formate file
      void updateCoordinates(std::istream &); // initialize the class from xyz formate file
      Molecule & setCenterOfMass0(); // get the center of mass to origin
      Molecule & transform(const Eigen::Matrix3d &);
      Molecule & translate(const Eigen::Vector3d &);
      void vibrateAlongMode(double );
      std::string repr();

  };
}
#endif // _MOLECULE_H
