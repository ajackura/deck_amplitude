#include "deck.h"
#include "pion_nucleon_scattering.h"
#include "deck_pi_n.h"

namespace d = deck;
namespace dp = deck_pi_n;
namespace pn = pion_nucleon_scattering;
int main()
{

  pn::SAID_pion_nucleon_pw ( );

  double M3pi_i = d::m_pi + d::m_rho + 0.0001;
  double M3pi_f = 7.0;
  int nstep_x = 200;

  double cos_th_i = -1.0;
  double cos_th_f =  1.0;
  int nstep_y = 200;

  int nstep = 500;
  
  d::out_deck ( );

  int J;
  int M;
  int S;
  int lam;

  //  std::cin >> J;
  //  std::cin >> M;
  //  std::cin >> S;
  //  std::cin >> lam;


  //  1.78668
  //1.79886
  
  double s = pow(1.79886,2);//std::pow(1.8,2); //1.7745

  double sigma = d::m_rho * d::m_rho;
  //  std::cout << d::projection ( J, M, S, lam, s, sigma ) << std::endl;

    //  d::plot_1d_projection ( J, M, S, lam, M3pi_i, M3pi_f, nstep );

  //  d::plot_1d_intensity ( cos_th_i,
  //  			 cos_th_f,
  //  			 nstep );

  /*
  d::plot_2d_intensity ( M3pi_i,
			 M3pi_f,
			 cos_th_i,
			 cos_th_f,
			 nstep_x,
			 nstep_y );
  */
  return 0;
}
