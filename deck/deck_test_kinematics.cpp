#include <iostream>
#include <fstream>
#include <cmath>
#include "deck_kinematics.h"
#include "deck_pi_n.h"

namespace dk = deck_kinematics;
namespace dp = deck_pi_n;
namespace pn = pion_nucleon_scattering;

int main()
{
  pn::SAID_pion_nucleon_pw ( );
  double s  =  340.0;
  double t  = -0.01;
  double Wi =  1.0;
  double Wf =  5.0;
  int Npts  =  100;
  //  dk::plot_Eabd_Pabd ( s, t, Wi, Wf, Npts );
  double s1 = 0.77;

  double z = -1.0;
  double phi = 2.0 * 3.1415926535;

  dp::plot_M_pp_p_mm (s, t, s1, z, phi, Wi, Wf, Npts);



  //  dk::plot_psi_tR ( s, t, s1, Wf, Wi, Npts );
  return 0;
}
