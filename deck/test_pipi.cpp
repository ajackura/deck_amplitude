#include "pion_pion_scattering.h"

namespace pp = pion_pion_scattering;

int main ()
{
  double mi = 2.0 * pp::m_pi;
  double mf = 2.5;
  int Npts  = 1000;
  int L = 0;
  int I = 0;
  pp::plot_pipi_amp ( L, I, mi, mf, Npts );
  return 0;
}
