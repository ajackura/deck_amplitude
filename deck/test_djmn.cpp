#include "rotational_functions.h"

namespace rf = rotational_functions;
int main()
{

  double J;
  double m_p;
  double m;

  //  std::cin >> J;
  //  std::cin >> m;
  //  std::cin >> m_p;

  //  rf::plot_wigner_d ( J, m_p, m );

  //  rf::print_clebsch_gordan ();

  //  rf::plot_legendre ( J );

  int l,ml;
  std::cin >> l;
  std::cin >> ml;
  rf::plot2d_spherical_harmonics ( l, ml );


  


  return 0;
}
