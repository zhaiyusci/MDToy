#ifndef _MDTOY_KABSCH_H
#define _MDTOY_KABSCH_H
#include"molecule.h"
namespace MDToy{
  extern Eigen::Matrix3d kabsch(Molecule, Molecule);
}
#endif // _MDTOY_KABSCH_H

