/****************************************************************************
   Name         :: pion_nucleon_scattering.h
   Author       :: Andrew W. Jackura, V. Mathieu
   Contact      :: ajackura@iu.edu
   Date         :: Feb 06, 2019
   
   Dependencies :: rotational_functions.h
   
   Description  :: Module 'namespace pion_nucleon_scattering',
                   contains functions to compute pion-nucleon, 
                   pi N --> pi N, scattering amplitudes and observables.
                   Details can be found in Ref. [1] and 
                   http://cgl.soic.indiana.edu/jpac/PiN.php.
                   The low-energy model is the SAID partial wave set [2],
                   while the high-energy model is a Regge parameterization
                   discussed in Ref. [1]. 

                   Some important kinematic variable definitions are
                   's = center-of-momentum energy-squared'
                   't = invariant momentum-transfer-squared'

                    pi (p1)   \  /  pi (p2)       s = ( p1 + p2 )^2
                               \/                 t = ( p2 - p4 )^2
                               /\
                     N (p2)   /  \   N (p4)

                   Included functions:
                   - legendre_polynomial ( n, x )
                   - legendre_derivative ( n, x )
                   - legendre_recursion  ( n, x )
                   - clebsch_gordan ( j1, m1, j2, m2, j, m )
                   - factorial10 (n)
                   - wigner_d_matrix ( j, m, n, beta )
                   - wigner_D_matrix ( j, m, n, alpha, beta, gamma )
                   - spherical_harmonics ( l, m, beta, alpha )
                   - cartesian_spherical_harmonics ( l, m, r )

   References   :: [1] V. Mathieu, I. V. Danilkin, C. Fernandez-Ramirez,
                       M. R. Pennington, D. Schott, A. P. Szczepaniak,
                       and G. Fox,
                       Phys. Rev. D 92, 074004 (2015)                          
                             
                   [2] R.L. Workman, R.A. Arndt, W.J. Briscoe,
                       M.W. Paris, I.I. Strakovsky,
                       Phys. Rev. C 86, 035202 (2012)   
****************************************************************************/
#ifndef PION_NUCLEON_SCATTERING_H_
#define PION_NUCLEON_SCATTERING_H_

#include<iostream>
#include<sstream>
#include<fstream>
#include<cstring>
#include<cmath>
#include<vector>
#include<complex>

#include "rotational_functions.h"

namespace rf = rotational_functions;

typedef std::complex<double> cd;
typedef std::vector<double> vd;
typedef std::vector<vd> v2d;
typedef std::vector<v2d> v3d;
typedef std::vector<cd> cvd;
typedef std::vector<cvd> cv2d;
typedef std::vector<cv2d> cv3d;

  /****************************************************************************
   Namespace :: pion_nucleon_scattering
   ***************************************************************************/
namespace pion_nucleon_scattering
{
  double m_N  = 0.9382720;    // nucleon mass
  double m_pi = 0.13957018;   // pion mass
  double hbc  = 0.19732;      // h-bar * c
  double pi = 3.1415926535897932384626433832795028841972;
  cd xr(1.0,0.0);
  cd xi(0.0,1.0);
  cd xz(0.0,0.0);
  vd regge_param ( { 0.4895696553113608,   // alpha rho
	             0.9425449473756354,
	             0.4895696553113608,   // alpha f
                     0.9425449473756354,
                     1.07473280868023,     // alpha Pomeron
	             0.4343313339478808,
	             0.16164225485369113,
                     5.00657330231329,     // C{0,1,2} rho
	            10.0982220619908,     
	             0.5730288256207643,   
                   128.86955863248116,     // D{0,1} rho
                     1.3840421677412782,
                    71.35385175012772,     // C{0,1} f
                     3.1790121902525574,
                    71.35385175012772,     // D{0,1} f
                     3.1790121902525574,
                    23.894196195009776,    // C{0,1} Pomeron
                     2.205005669899481,
                    23.894196195009776,    // D{0,1} Pomeron
                     2.205005669899481 } );


  vd plab;
  v3d RePW( 4, v2d( 10 ) );
  v3d ImPW( 4, v2d( 10 ) );
  cv3d PW( 4, cv2d( 10 ) );

  /****************************************************************************
   Function declarations
   ***************************************************************************/
  void total_cross_section ( double s,
			     vd &sigma_tot );
  
  void observables ( double s,
		     double t,
		     v2d &obs );

  void Full_ABC ( double s,
                  double t,
                  cv2d & ABC );

  void Regge_ABC ( double s,
		   double t,
		   cv2d &ABC );
  
  void SAID_ABC ( double s,
                  double t,
                  cv2d &ABC );

  void SAID_pion_nucleon_pw ( );

  double legendre ( int n,
		    double x );
  
  double legendre_derivative ( int n,
			       double x );

  cd cgamma ( cd z,
	      int OPT );

  void plot_legendre ( double xi,
		       double xf,
		       int nsteps,
		       int n );
 
  void plot_legendre_derivative ( double xi,
				  double xf,
				  int nsteps,
				  int n );

  void plot_gamma ( double xi,
                    double xf,
                    int nsteps );

  void plot_pi_N_cross ( double si,
                         double sf,
                         int nsteps );

  void plot_pi_N_ABC ( double si,
                       double sf,
                       double t,
                       int nsteps );

  
  /****************************************************************************
   Function total_cross_section ( s, sigma_tot )
   
   Computes the total cross section for pi N -> pi N

   ***************************************************************************/
  void total_cross_section ( double s,
			     vd &sigma_tot )
  {
    cv2d ABC(2, cvd (3,0) );

    double t    = 0.;
    double Elab = ( s - m_N * m_N - m_pi * m_pi ) / 2. / m_N;
    double Plab = sqrt ( Elab * Elab - m_pi * m_pi );
    
    Full_ABC ( s, t, ABC );
    sigma_tot[0] = 0.389352 * imag ( ABC[0][2] + ABC[1][2] ) / Plab;
    sigma_tot[1] = 0.389352 * imag ( ABC[0][2] - ABC[1][2] ) / Plab;
    
    return;
  }

  /***************************************************************************
   Function observables ( s, t, obs )
   ***************************************************************************/
  void observables ( double s,
		     double t,
		     v2d & obs )
  {
    cv2d ABC(2, cvd(3,0) );
    cd B(0,0);
    cd C(0,0);

    double Elab = ( s - m_N * m_N - m_pi * m_pi ) / 2.0 / m_N;
    double Plab = sqrt ( Elab * Elab - m_pi * m_pi );
    double sth  = pow( m_N * m_N + m_pi * m_pi, 2 );
    double spth = pow( m_N * m_N - m_pi * m_pi, 2 );
    double qcm  = sqrt ( ( s - sth ) * ( s - spth ) ) / 2.0 / sqrt( s );
    double cost = 1.0 + t / 2.0 / qcm / qcm;
    double sint = sqrt ( 1.0 - cost * cost );
    double fac  = 1.0 / pi / s * pow ( m_N / 4.0 / qcm, 2 );
    double tau  = t / 4.0 / m_N / m_N;

    Full_ABC ( s, t, ABC );

    // pi^- p --> pi^- p
    B = ABC[0][1] + ABC[1][1];
    C = ABC[0][2] + ABC[1][2];
    
    obs[0][0] = fac * ( ( 1.0 - tau ) * abs(C) * abs(C)
			- tau *( ( s * tau + Plab * Plab )/( 1.0 - tau ) ) 
			* abs(B) * abs(B) );
    obs[0][1] = -sint / ( 16.0 * pi * sqrt(s) ) * std::imag( C * std::conj(B) )
      / obs[0][0];
    obs[0][0] = obs[0][0] * 389.352; 

    // pi^+ p --> pi^+ p
    B = ABC[0][1] - ABC[1][1];
    C = ABC[0][2] - ABC[1][2];

    obs[1][0] = fac * ( ( 1.0 - tau ) * abs(C) * abs(C)
                        - tau *( ( s * tau + Plab * Plab )/( 1.0 - tau ) )
                        * abs(B) * abs(B) );
    obs[1][1] = -sint / ( 16.0 * pi * sqrt(s) ) * std::imag( C * std::conj(B) )
      / obs[1][0];
    obs[1][0] = obs[1][0] * 389.352;

    // pi^- p --> pi^0 n
    B = - sqrt( 2.0 ) * ABC[1][1];
    C = - sqrt( 2.0 ) * ABC[1][2];
    
    obs[2][0] = fac * ( ( 1.0 - tau ) * abs(C) * abs(C)
                        - tau *( ( s * tau + Plab * Plab )/( 1.0 - tau ) )
                        * abs(B) * abs(B) );
    obs[2][1] = -sint / ( 16.0 * pi * sqrt(s) ) * std::imag( C * std::conj(B) )
      / obs[2][0];
    obs[2][0] = obs[2][0] * 389.352;

    return;
  }



  /***************************************************************************

   **************************************************************************/
  void Full_ABC ( double s, 
		  double t,
		  cv2d &ABC )
  {
    cv2d ABCSaid(2, cvd(3,0) );
    cv2d ABCRegge(2, cvd(3,0) );

    double thresh1 = 1.5;
    double thresh2 = 2.1;
    double Elab = ( s - m_N * m_N - m_pi * m_pi ) / 2.0 / m_N;

    if ( Elab < thresh1 )
      {
	SAID_ABC( s, t, ABC );
	return;
      }
    if ( Elab > thresh2 )
      {
	Regge_ABC( s, t, ABC );
      }
    if ( Elab >= thresh1 && Elab <= thresh2 )
      {
	SAID_ABC( s, t, ABCSaid );
	Regge_ABC( s, t, ABCRegge );
	
	double eps = ( Elab - thresh1 ) / ( thresh2 - thresh1 );

	ABC[0][0] = (1 - eps) * ABCSaid[0][0] + eps * ABCRegge[0][0];
	ABC[0][1] = (1 - eps) * ABCSaid[0][1] + eps * ABCRegge[0][1];
	ABC[0][2] = (1 - eps) * ABCSaid[0][2] + eps * ABCRegge[0][2];

	ABC[1][0] = (1 - eps) * ABCSaid[1][0] + eps * ABCRegge[1][0];
	ABC[1][1] = (1 - eps) * ABCSaid[1][1] + eps * ABCRegge[1][1];
	ABC[1][2] = (1 - eps) * ABCSaid[1][2] + eps * ABCRegge[1][2];

	return;
      }
  }
  


  void Regge_ABC ( double s,
		   double t,
		   cv2d &ABC )
  {
    // the parameters:
    double a0rho = regge_param[0]; 
    double aprho = regge_param[1];
    double a0f = regge_param[2]; 
    double apf = regge_param[3];
    double a0P = regge_param[4]; 
    double apP = regge_param[5]; 
    double appP = regge_param[6];
    double C0rho = regge_param[7]; 
    double C1rho = regge_param[8]; 
    double C2rho = regge_param[9];
    double D0rho = regge_param[10]; 
    double D1rho = regge_param[11];
    double C0f = regge_param[12]; 
    double C1f = regge_param[13];
    double D0f = regge_param[14]; 
    double D1f = regge_param[15];
    double C0P = regge_param[16]; 
    double C1P = regge_param[17];
    double D0P = regge_param[18]; 
    double D1P = regge_param[19];


    double Elab = ( s - m_N * m_N - m_pi * m_pi ) / 2.0 / m_N;
    double som  = 2.0 * m_N * m_N + 2.0 * m_pi * m_pi;
    double u    = - s - t + som;
    double nu   = ( s - u ) / 4.0 / m_N;

      // Trajectories
    cd alpRho = a0rho + t * aprho;
    cd alpF = a0f + t * apf;
    cd alpPom = a0P + t * apP + t * t * appP;

    // Regge factors:
    cd R1rho = cgamma( 0.0 - alpRho, 0)/2.0 
      * ( 1.0-std::exp(-xi * pi *alpRho) ) * std::pow(nu,alpRho);
    cd R2rho = cgamma( 1.0 - alpRho, 0)/2.0 
      * ( 1.0-std::exp(-xi * pi * alpRho) ) * std::pow(nu,alpRho-1.0);
    cd R1f   = cgamma( 1.0 - alpF, 0)/2.0 
      * ( 1.0+std::exp(-xi * pi * alpF) ) * std::pow(nu,alpF);
    cd R2f   = cgamma( 1.0 - alpF, 0)/2.0 
      * ( 1.0+std::exp(-xi * pi * alpF) ) * std::pow(nu,alpF-1.0);
    cd R1Pom = cgamma( 1.0 - alpPom, 0)/2.0 
      * ( 1.0+std::exp(-xi * pi * alpPom) ) * std::pow(nu,alpPom);
    cd R2Pom = cgamma( 1.0 - alpPom, 0)/2.0 
      * ( 1.0+std::exp(-xi * pi * alpPom) ) * std::pow(nu,alpPom-1.0);

    cd Crho = - C0rho*( (1.0+C2rho)*std::exp(C1rho*t) - C2rho )* R1rho;
    cd Brho = + D0rho * std::exp(D1rho*t)* R2rho;
    cd Cf   = - C0f * std::exp(C1f*t) * R1f;
    cd Bf   = - D0f * std::exp(D1f*t) * R2f;
    cd CPom = - C0P * std::exp(C1P*t) * R1Pom;
    cd BPom = - D0P * std::exp(D1P*t) * R2Pom;

    // reconstruct A from B and C = A'
    cd Arho = Crho - (Elab + t/(4*m_N) )/(1 - t/(4*m_N*m_N))*Brho;
    cd Af   = Cf   - (Elab + t/(4*m_N) )/(1 - t/(4*m_N*m_N))*Bf;
    cd APom = CPom - (Elab + t/(4*m_N) )/(1 - t/(4*m_N*m_N))*BPom;

    // Scalar Amplitudes (A,B,C)
    ABC[0][0] = APom + Af;
    ABC[0][1] = BPom + Bf;
    ABC[0][2] = CPom + Cf;

    ABC[1][0] = Arho;
    ABC[1][1] = Brho;
    ABC[1][2] = Crho;


    return;
  }


  void SAID_ABC ( double s,
		  double t,
		  cv2d &ABC )
  {
    double w    = sqrt ( s );
    double E    = ( s + m_N * m_N - m_pi * m_pi ) / 2.0 / w;
    double sth  = pow( m_N  + m_pi, 2 );
    double spth = pow( m_N  - m_pi, 2 );
    double qcm  = sqrt ( ( s - sth ) * ( s - spth ) ) / 2.0 / sqrt( s );
    double z    = 1.0 + t / 2.0 / qcm / qcm;
    double Tlab = ( s - m_N * m_N - m_pi * m_pi ) / 2.0 / m_N - m_pi;
    double Elab = Tlab + m_pi;
    double Plab = sqrt ( Elab * Elab - m_pi * m_pi );

    //    std::cerr << w << " " << E << " " << qcm << " " << z << " " << 
    //      Tlab << " " << Elab << " " << Plab << std::endl;

    if ( Plab > 2.5 )
      {
	return;
      }
    cvd PW1m(10);
    cvd PW1p(10);
    cvd PW3m(10);
    cvd PW3p(10);
    
    int k = std::floor(Plab/0.025);
    double eps = ( Plab - plab[k] ) / 0.025;


    for ( int L = 0; L < 8; L++ )
      {
	PW1m[L] = (1.0-eps)*PW[0][L][k] + eps*PW[0][L][k+1];
	PW1p[L] = (1.0-eps)*PW[1][L][k] + eps*PW[1][L][k+1];
	PW3m[L] = (1.0-eps)*PW[2][L][k] + eps*PW[2][L][k+1];
	PW3p[L] = (1.0-eps)*PW[3][L][k] + eps*PW[3][L][k+1];
      }


    //    std::cerr << "PW " << PW1m[3] << " " << PW1p[3] << std::endl;
    
    cd fp(0,0);
    cd fm(0,0);
    cd gp(0,0);
    cd gm(0,0);
    cd f1p(0,0);
    cd f1m(0,0);
    cd f2p(0,0);
    cd f2m(0,0);
    cd f1(0,0);
    cd g1(0,0);
    cd f3(0,0);
    cd g3(0,0);

    double pol=0.0;
    double pol1d=0.0;

    for ( int L = 0; L < 8; L++ )
      {
	pol = rf::legendre_polynomial ( L, z );
	pol1d = rf::legendre_derivative ( L, z );
	double ll = (double) L;
	f1 = f1 + ( ( ll + 1.0 ) * PW1p[L] + ll * PW1m[L] ) * pol;
	f3 = f3 + ( ( ll + 1.0 ) * PW3p[L] + ll * PW3m[L] ) * pol;
	g1 = g1 + ( PW1p[L] - PW1m[L] ) * pol1d;
	g3 = g3 + ( PW3p[L] - PW3m[L] ) * pol1d;
      }

    f1 = f1 / qcm;
    f3 = f3 / qcm;
    g1 = g1 / qcm;
    g3 = g3 / qcm;



    fp = ( f1 + 2.0 * f3 ) / 3.0;
    gp = ( g1 + 2.0 * g3 ) / 3.0;
    fm = ( f1 - f3 ) / 3.0;
    gm = ( g1 - g3 ) / 3.0;

    f2p = - gp;
    f1p =   fp - z * f2p;
    f2m = - gm;
    f1m =   fm - z * f2m;

    ABC[0][0] = 4* pi * ( (w+m_N)/(E+m_N)*f1p - (w-m_N)/(E-m_N)*f2p  );
    ABC[1][0] = 4* pi * ( (w+m_N)/(E+m_N)*f1m - (w-m_N)/(E-m_N)*f2m  );

    ABC[0][1] = 4* pi * ( 1.0/(E+m_N)*f1p + 1.0/(E-m_N)*f2p  );
    ABC[1][1] = 4* pi * ( 1.0/(E+m_N)*f1m + 1.0/(E-m_N)*f2m  );

    ABC[0][2] =  ABC[0][0] + (Elab + t/(4*m_N) )/(1.0 - t/(4*m_N*m_N) )*ABC[0][1];
    ABC[1][2] =  ABC[1][0] + (Elab + t/(4*m_N) )/(1.0 - t/(4*m_N*m_N) )*ABC[1][1];

    return;
  }


  /*************************************************************************** 
   Function SAID_pion_nucleon_pw ( );

   

   The SAID WI08 solution is valid up to Plab = 2.5 GeV
   The SAID solution includes 8 waves for both parity
   and s-channel isospin - In the spectroscopic notation L(2I)(2J) 

        P11  D13  F15  G17   H19  I111  J113
   S11  P13  D15  F17  G19  H111  I113  J115
        P31  D33  F35  G37   H39  I311  J313
   S31  P33  D35  F37  G39  H311  I313  J315


   The SAID solution computed for the reaction
   pi N --> pi N
   for all isospin compibation
   
   return wval, the vector of W = Sqrt{s} values
   return double complex PW[5][10]
   PW[0][L] is I = 1/2; J = L-1/2
   PW[1][L] is I = 1/2; J = L+1/2
   PW[2][L] is I = 3/2; J = L-1/2
   PW[3][L] is I = 3/2; J = L+1/2
   
   'plab' is plab the momentum of the pion in the lab
   
   The last index is the index of W
   MAXIMUM L is 8
   The partial waves are dimensionsless
   ***************************************************************************/
  void SAID_pion_nucleon_pw ( ) 
  {
    std::string SAIDfileslist = "SAID_PW_list.txt";
    std::vector<std::string> SAIDfilename;

    std::ifstream inputFile(SAIDfileslist);

    if (inputFile)
      {
	std::string name;
	while ( inputFile >> name )
	  {
	    SAIDfilename.push_back(name);
	  }
      }

    //    vd plab; //( 101 );
    //    v3d RePW( 4, v2d( 10 ) ); // vd( 101 ) ) );
    //    v3d ImPW( 4, v2d( 10 ) ); // , vd( 101 ) ) );

    //    cv3d PW( 4, cv2d( 10 ) ); 

    std::ifstream SAIDinput(SAIDfilename[0]);
    {
      if (SAIDinput)
	{
	  double tmp;
	  double val;
	  while (SAIDinput >> val >> tmp >> tmp >> tmp >> tmp
		 >> tmp >> tmp >> tmp >> tmp )
	    {
	      plab.push_back(val / 1000.0); // convert to GeV
	    }
	}
    }
    

    for ( int j = 0; j <= 1; j++ )
      {
	std::ifstream SAIDinput(SAIDfilename[j]);

	if (SAIDinput)
	  {
	    double tmp;
	    double val1;
	    double val2;
	    while ( SAIDinput >> tmp >> tmp >> tmp >> tmp >> tmp 
		    >> val1 >> val2 >> tmp >> tmp )
	      {
		cd val(val1,val2);
		PW[2 * j + 1][0].push_back(val);
		PW[2 * j][0].push_back(xz);
		//		std::cout << "val2 = " << val2 << std::endl;
		RePW[2 * j + 1][0].push_back(val1);
		ImPW[2 * j + 1][0].push_back(val2);
		
	      }
	  }
      }
  

    for ( int L = 1; L <= 7; L++ )
      {
	for ( int j = 0; j <= 3; j++ )
	  { 
	    std::ifstream SAIDinput(SAIDfilename[4 * L + j - 2]);

	    if (SAIDinput)
	      {
		double tmp;
		double val1;
		double val2;
		while ( SAIDinput >> tmp >> tmp >> tmp >> tmp >> tmp
			>> val1 >> val2 >> tmp >> tmp )
		  {
		    cd val(val1,val2);
		    PW[j][L].push_back(val);
		    RePW[j][L].push_back(val1);
		    ImPW[j][L].push_back(val2);
		  }
	      }
	  }
      }
    return;
  }


  void plot_SAID_pw ( int j,
		      int L )
  {
    for (int i = 0; i < plab.size(); i++ )                            
      {                                                                
	std::cout << plab[i] << " " << std::real(PW[j][L][i]) << " " 
		  << std::imag(PW[j][L][i]) << std::endl;
      }
    return;
  }


  /****************************************************************************
   Function legendre ( n, x )
    
   Computes the Legendre polynomial P_n (x)
   using the Bonnet recursion formula
    
   (n+1) P_{n+1}(x) = (2n+1) x P_{n}(x) - n P_{n-1}(x) for n > 1,
    
   with P_0 = 1, P_1 = x
   ***************************************************************************/
  double legendre ( int n,
		    double x )
  {
    double r = 0.; // r = P_0
    double s = 1.; // s = P_1
    double t;
    int m;
    for ( m = 0; m < n; m++ )
      {
	t = r;
	r = s;
	s = ( 2. * m + 1. ) * x * r - m * t;
	s /= ( m + 1. );
      } 
    return s;
  } 

  /****************************************************************************
   Function legendre_derivative ( n, x )

   Computes the derivative of the Legendre polynomial P_n' (x),
   where ()' = d/dx, using the formula

   ( ( x^2 - 1 ) / n ) P_n'(x) = x P_n(x) - P_{n-1}(x) for n > 0,
  
   with P_0' = 0
   ***************************************************************************/
  double legendre_derivative ( int n,
			       double x )
  {
    if ( n == 0 )
      {
	return 0;
      }
    if ( n == 1 )
      {
	return 1.0;
      }
    if ( n == 2 )
      {
	return 3.0 * x;
      }
    if ( n == 3 )
      {
	return (-3 + 15*x*x)/2.0;
      }
    if ( n == 4 )
      {
	return (-60*x + 140*x*x*x)/8;
      }
    if ( n == 5 )
      {
	return (15 - 210*x*x + 315*x*x*x*x)/8;
      }
    if ( n == 6 )
      {
	return (210*x - 1260*x*x*x + 1386*pow(x,5))/16;
      }
    if ( n == 7 )
      {
	return (-35 + 945*x*x - 3465*pow(x,4) + 3003*pow(x,6))/16;
      }
    if ( n == 8 )
      {
	return (-2520*x + 27720*x*x*x - 72072*pow(x,5) + 51480*pow(x,7))/128;
      }
    if ( n == 9 )
      {
	return (315 - 13860*x*x + 90090*pow(x,4) - 180180*pow(x,6) + 109395*pow(x,8))/128;
      }
    //    else
    //   {
    //	double d = x * legendre ( n, x ) - legendre ( n-1, x );
    //	d = d * n / ( pow ( x, 2 ) - 1. );
    //	return d;
    //  }
  }
  /****************************************************************************
   Function plot_legendre ( xi, xf, nsteps, n )
     
   plotting routine for the Legendre polynomial of order 'n'
   between the bounds 'xi' and 'xf' in 'nsteps' equal steps
     
   ***************************************************************************/
  void plot_legendre ( double xi,
		       double xf,
		       int nsteps,
		       int n )
  {
    double step = ( xf - xi ) / ( (double) nsteps );
    double x    = xi;
    for ( int i = 1; i <= nsteps; i++ )
      {
	std::cout << x << " " << legendre(n,x) << std::endl;
	x = x + step;
      }
  }

  /****************************************************************************
   Function plot_legendre_derivative ( xi, xf, nsteps, n )

   plotting routine for the derivative of the Legendre polynomial
   of order 'n' between the bounds 'xi' and 'xf' in 'nsteps' equal steps   

   ***************************************************************************/
  void plot_legendre_derivative ( double xi,
				  double xf,
				  int nsteps,
				  int n )
  {
    double step = ( xf - xi ) / ( (double) nsteps );
    double x    = xi;
    for ( int i = 1; i <= nsteps; i++ )
      {
	std::cout << x << " " << legendre_derivative(n,x) << std::endl;
	x = x + step;
      }
  }


  /****************************************************************************
   Function cgamma ( z, OPT )

   Computes the complex gamma function for complex z

   ***************************************************************************/
  cd cgamma ( cd z, int OPT )
  {
    cd g;
    cd infini(1e308,0.0); // z0,z1
    double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
    double na,t,x1 = 1,y1,sr,si;
    int j,k;
    

    static vd a( { 8.333333333333333e-02,
	          -2.777777777777778e-03,
	           7.936507936507937e-04,
	          -5.952380952380952e-04,
	           8.417508417508418e-04,
	          -1.917526917526918e-03,
	           6.410256410256410e-03,
	          -2.955065359477124e-02,
	           1.796443723688307e-01,
	          -1.39243221690590 } );
    
    x = std::real(z);
    y = std::imag(z);

    if (x > 171)
      {
	return infini;
      }
    if ( ( y == 0.0 ) && ( x == (int) x ) && ( x <= 0.0 ) )
      {
	return infini;
      }
    if (x < 0.0) 
      {
	x1 =  x;
	y1 =  y;
	x  = -x;
	y  = -y;
      }
    x0 = x;
    if ( x <= 7.0 ) 
      {
	na = (int)( 7.0 - x );
	x0 = x+na;
      }
    q1 = std::sqrt( x0 * x0 + y * y );
    th = std::atan( y / x0 );
    gr = ( x0 - 0.5 )*std::log( q1 ) - th * y 
         - x0 + 0.5 * std::log( 2.0 * pi );
    gi = th * ( x0 - 0.5 ) + y * std::log( q1 ) - y;
    for ( k = 0; k < 10; k++ )
      {
	t = pow( q1, -1.0 -2.0 * k );
	gr += ( a[k] * t * std::cos( ( 2.0 * k + 1.0 ) * th ) );
	gi -= ( a[k] * t * std::sin( ( 2.0 * k + 1.0 ) * th ) );
      }
    if ( x <= 7.0 ) 
      {
	gr1 = 0.0;
	gi1 = 0.0;
	for ( j = 0; j < na; j++ )
	  {
	    gr1 += ( 0.5 * std::log( ( x+j ) * ( x+j ) + y * y ) );
	    gi1 += std::atan( y / ( x + j ) );
	  }
	gr -= gr1;
	gi -= gi1;
      }


    if ( x1 <= 0.0 ) 
      {
	q1  =  std::sqrt( x * x + y * y );
	th1 =  std::atan( y / x );
	sr  = -std::sin( pi * x ) * std::cosh( pi * y );
	si  = -std::cos( pi * x ) * std::sinh( pi * y );
	q2  =  std::sqrt( sr * sr + si * si );
	th2 =  std::atan( si / sr );
	if ( sr < 0.0 ) 
	  {
	    th2 += pi;
	  }
	gr =  std::log( pi / ( q1 * q2 ) ) - gr;
	gi = -th1 - th2 - gi;
	x  =  x1;
	y  =  y1;
      }

    if ( OPT == 0 ) 
      {
      g0 = std::exp( gr );
      gr = g0 * std::cos( gi );
      gi = g0 * std::sin( gi );
      }
    g = cd( gr, gi );
    return g;
  }

  /****************************************************************************
   Function plot_gamma ( xi, xf, nsteps )
   

   ***************************************************************************/
  void plot_gamma ( double xi,
                    double xf,
                    int nsteps )
  {
    double x = xi;
    double step = ( ( xf - xi ) / ( (double) nsteps ) );
    cd cx;
    cd g;
    int zero = 0;
    for ( int i = 1; i < nsteps; i++ )
      {
        cx = cd(x,0.0);
        g = cgamma(cx,zero);
	std::cout << x << " " << real(g) << " " << imag(g) << std::endl;
        x = x + step;
      }
  }

  

  void plot_pi_N_cross ( double Pi,
			 double Pf,
			 int nsteps )
    {
      double P = Pi;
      double step = ( ( Pf - Pi ) / ( (double) nsteps ) );
      double s;
      vd sigma(2,0.0);

      for ( int i = 1; i < nsteps; i++ )
	{
	  s = m_pi * m_pi + m_N * m_N + 
	    2.0 * m_N * sqrt( m_pi * m_pi +  P * P );
	  total_cross_section (s, sigma);
	  std::cout << P << " " << sigma[0] << " " << sigma[1] << std::endl;
	  P = P + step;
	    
	}
    }


  void plot_pi_N_ABC ( double si,
                       double sf,
                       double t,
                       int nsteps )
  {

    cv2d ABC(2,cvd(3));
    double s = si;
    double step = ( ( sf - si ) / ( (double) nsteps ) );

    for ( int i = 1; i < nsteps; i++ )
      {
        Full_ABC ( s, t, ABC );
	std::cout << s << " " << real(ABC[0][0]) << " " << imag(ABC[0][0])
                  << " " << real(ABC[0][1]) << " " << imag(ABC[0][1]) << " "
                  << real(ABC[0][2]) << " " << imag(ABC[0][2]) << " "
                  << real(ABC[1][0]) << " " << imag(ABC[1][0]) << " "
                  << real(ABC[1][1]) << " " << imag(ABC[1][1]) << " "
                  << real(ABC[1][2]) << " " << imag(ABC[1][2]) << std::endl;
        s = s + step;
      }
    return;
  }

}

#endif
