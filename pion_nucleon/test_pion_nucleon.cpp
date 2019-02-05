#include "pion_nucleon_scattering.h"

namespace pna = pion_nucleon_scattering;

int main()
{
  double xi = -5;
  double xf = 5;
  int nsteps = 100000;
  int n = 2;
  //  pna::plot_legendre_derivative (xi,xf,nsteps,n);

  pna::SAID_pion_nucleon_pw ( );


  double Pi = 0.1;
  double Pf = 1000.0;
  
  double si = 1.25;
  double sf = 6.0;
  double t = 0.0;



  pna::plot_pi_N_ABC ( si, sf, t, nsteps);

  //  pna::plot_pi_N_cross ( Pi, Pf, nsteps);

  //  pna::plot_gamma ( xi, xf, nsteps );


  return 0;
}
