#include "deck_amplitudes.h"
#include "utilities.h"
#include "deck_observables.h"

namespace d = deck_amplitudes;
namespace pn = pion_nucleon_scattering;
namespace u = utilities;
namespace dko = deck_observables;

int main()
{

  pn::SAID_pion_nucleon_pw ( );

  const double E_beam   = 10.0; //190                                        
  const double Plab     = 190.0;//16.0;
  const double s       = d::sq_m_pi + d::sq_m_p + 2.0 * d::m_p * std::sqrt ( d::sq_m_pi + Plab * Plab );//E_beam;

  const double t        = -0.01;//-0.05;//-0.1; //-0.1


  double W_i = 3.0 * d::m_pi + 0.0001;//d::m_pi + d::m_rho + 0.0001;
  double W_f = 2.5;
  int Npts = 60;

  double cos_th_i = -1.0;
  double cos_th_f =  1.0;
  int nstep_y = 200;

  int nstep = 500;
  

  int J;
  int M;
  int L;
  int S;

  //  J = 0;
  //  M = 0;
  //  L = 1;
  //  S = 1;

  std::cin >> J;
  std::cin >> M;
  std::cin >> L;
  std::cin >> S;

  /////////////////  d::out_deck ( );

  //  double s = 380.0;
  //  double t = -0.01;
  //double s1 = d::m_rho * d::m_rho;
  //  double W = 1.5;
  int mu = 1;
  int mu_p = 1;

  //  double W = 2.30916;
  double W = 1.1;
  double s1 = 0.820 * 0.820;
  //  std::cout << d::Famp ( W, s1, s, t, J, M, L, S, mu, mu_p );
  //d::plot_1d_projection ( s, t, s1, mu, mu_p, J, M, L, S, W_i, W_f, Npts );

  dko::plot_1d_intensity ( s, t, mu, mu_p, J, M, L, S, W_i, W_f, Npts );

  double mag;
  double phi;
  //  cd amp = d::Famp ( W, s1, s, t, J, M, L, S, mu, mu_p );
  //  u::convert_complex ( amp, mag, phi );
  //  std::cerr << J << " " << M << " " << L << " " << S
  //            << " " <<mag << " " << phi << std::endl;
  //  std::cout << J << " " << M << " " << L << " " << S 
  //	    << " " <<mag << " " << phi << std::endl;

  

  return 0;
}
