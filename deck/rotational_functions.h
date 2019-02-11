/****************************************************************************
   Name         :: rotational_functions.h
   Author       :: Andrew W. Jackura
   Contact      :: ajackura@iu.edu
   Date         :: Feb 06, 2019
   
   Dependencies :: none
   
   Description  :: Module 'namespace rotational_functions', which contains
                   various functions pertaining to rotations and 
                   angular momenta. 
   
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

   References   :: [1] Quantum Theory Of Angular Momemtum,
                       D.A. Varshalovich
  ****************************************************************************/
#ifndef ROTATIONAL_FUNCTIONS_H_
#define ROTATIONAL_FUNCTIONS_H_

#include <iostream>
#include <complex>
#include <algorithm>
#include <cmath>
#include <vector>

typedef std::complex<double> cd;
typedef std::vector<double> vd;

 /****************************************************************************
  Namespace :: rotational_functions
  ****************************************************************************/
namespace rotational_functions
{
 /**************************************************************************** 
   Global Variables
  ****************************************************************************/
  const double pi    = 3.1415926535897932384626433832795028841972;
  const double twopi = 2.0 * pi;
  const cd ii (0.0,1.0); // complex 'i'
  const cd ir (1.0,0.0); // complex '1'
  const cd iz (0.0,0.0); // complex '0'
  /***************************************************************************
   Function Prototypes                                                         
  ****************************************************************************/
  double legendre_polynomial ( int n,
			       double x );

  double clebsch_gordan ( double j1,
			  double m1,
			  double j2,
			  double m2,
			  double j,
			  double m );

  double factorial10 ( int n );

  double wigner_d_matrix ( double j,
			   double m,
			   double n,
			   double beta );

  cd wigner_D_matrix (double j,
		      double m,
		      double n,
		      double alpha,
		      double beta,
		      double gamma );

  cd spherical_harmonics ( int l, 
			   int m, 
			   double beta, 
			   double alpha );

  cd cartesian_spherical_harmonics ( vd &rvec,
				     int l,
				     int m );
  
  /****************************************************************************
   Function    :: legendre_polynomial 

   Input       :: int n
                  double x

   Output      :: double 

   Description :: Returns the Legendre polynomial of degree 'n'
                  given an argument 'x', P_{n}(x)
  ****************************************************************************/
  double legendre_polynomial ( int n,
		               double x )
  {
    double x2  = std::pow ( x, 2 ) ;
    double x3  = std::pow ( x, 3 ) ;
    double x4  = std::pow ( x, 4 ) ;
    double x5  = std::pow ( x, 5 ) ;
    double x6  = std::pow ( x, 6 ) ;
    double x7  = std::pow ( x, 7 ) ;
    double x8  = std::pow ( x, 8 ) ;
    double x9  = std::pow ( x, 9 ) ;
    double x10 = std::pow ( x, 10 ) ;
    double Pl;
    switch ( n )
      {
      case 0: 
	Pl = 1.0; 
	break;
      case 1: 
	Pl = x;
	break;
      case 2: 
	Pl = ( - 1.0 + 3.0 * x2 ) / 2.0;
	break;
      case 3: 
	Pl = ( - 3.0 * x + 5.0 * x3 ) / 2.0;
	break;
      case 4: 
	Pl = (   3.0 - 30.0 * x2 + 35.0 * x4 ) / 8.0;
	break;
      case 5: 
	Pl = (   15.0 * x - 70.0 * x3 + 63.0 * x5 ) / 8.0;
	break;
      case 6: 
	Pl = ( - 5.0 + 105.0 * x2 - 315.0 * x4 + 231.0 * x6 ) / 16.0;
	break;
      case 7: 
	Pl = ( - 35.0 * x + 315.0 * x3 
	       - 693.0 * x5 + 429.0 * x7 ) / 16.0;
	break;
      case 8: 
	Pl = (   35.0 - 1260 * x2 + 6930.0 * x4 - 12012.0 * x6
		 + 6435.0 * x8 ) / 128.0;
	break;
      case 9: 
	Pl = (   315.0 * x - 4620.0 * x3 + 18018.0 * x5
		 - 25740.0 * x7 + 12155.0 * x9 ) / 128.0;
	break;
      case 10: 
	Pl = ( - 63.0 + 3465.0 * x2 - 30030.0 * x4 + 90090.0 * x6
	       - 109395.0 * x8 + 46189.0 * x10 ) / 256.0;
	break;
      default: 
	std::cerr << "Error in 'rotational_functions' " << std::endl;
	std::cerr << "n > 10 in function 'legendre_polynomial' " << std::endl;
	exit(1);
      }
    return Pl;
  }
  /****************************************************************************
   Function    :: legendre_derivative

   Input       :: int n
                  double x

   Output      :: double 

   Description :: Returns the derivative of the Legendre polynomial 
                  of degree 'n' given an argument 'x', (d/dx) P_{n}(x)
  ****************************************************************************/
  double legendre_derivative ( int n,
                               double x )
  {
    double x2 = std::pow ( x, 2 ) ;
    double x3 = std::pow ( x, 3 ) ;
    double x4 = std::pow ( x, 4 ) ;
    double x5 = std::pow ( x, 5 ) ;
    double x6 = std::pow ( x, 6 ) ;
    double x7 = std::pow ( x, 7 ) ;
    double x8 = std::pow ( x, 8 ) ;
    double x9 = std::pow ( x, 9 ) ;
    
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
        return ( - 3.0 + 15.0 * x2 ) / 2.0;
      }
    if ( n == 4 )
      {
        return ( - 60.0 * x + 140.0 * x3 ) / 8.0;
      }
    if ( n == 5 )
      {
        return (   15.0 - 210.0 * x2 + 315.0 * x4 ) / 8.0;
      }
    if ( n == 6 )
      {
        return (   210.0 * x - 1260.0 * x3 + 1386.0 * x5 ) / 16.0;
      }
    if ( n == 7 )
      {
        return ( - 35.0 + 945.0 * x2 - 3465.0 * x4 + 3003.0 * x6 ) / 16.0;
      }
    if ( n == 8 )
      {
        return ( - 2520 * x + 27720.0 * x3 - 72072.0 * x5 
		 + 51480.0 * x7 ) / 128.0;
      }
    if ( n == 9 )
      {
        return (   315.0 - 13860.0 * x2 + 90090.0 * x4 
		 - 180180.0 * x6 + 109395.0 * x8 ) / 128.0;
      }
    if ( n == 10 )
      {
        return (   6930.0 * x - 120120.0 * x3 + 540540.0 * x5
                 - 875160.0 * x7 + 461890.0 * x9 ) / 256.0;
      }
  }
  /****************************************************************************
   Function    :: legendre_recursion

   Input       :: int n
                  double x

   Output      :: double 

   Description :: Returns the Legendre polynomial of degree 'n'
                  given an argument 'x', P_{n}(x),
                  via Bonnet’s recursion formula
  ****************************************************************************/
  double legendre_recursion ( int n,
                              double x )
  {
    double r, s, t;
    int m;
    r = 0; s = 1.;
    // Use recursion relation to generate p1 and p2
    for (m=0; m < n; m++ )
      {
        t = r; r = s;
        s = (2*m+1) * x * r - m * t;
        s /= (m+1);
      } // end of do loop
    return s;
  }
  /****************************************************************************
   Function    :: clebsch_gordan

   Input       :: double j1 
                  double m1
                  double j2
                  double m2
                  double j
                  double m

   Output      :: double 

   Description :: Returns the Clebsch-Gordan coefficients 
                  < j m | j1 m1 ; j2 m2 >, given individual
                  angular momenta and projection '(j1,m1)' and '(j2,m2)'
                  which are coupled to total angular momenta '(j,m)'
  ****************************************************************************/
  double clebsch_gordan( double j1,
			 double m1,
			 double j2,
			 double m2,
			 double j,
			 double m )
  {
    if ( ( abs(m1) > j1 ) || ( abs(m2) > j2 ) || ( abs(m) > j ) )
      {
	//	std::cerr << "input error: |m| larger than j" << std::endl;
	return 0.0;
      }

    double epsil = 1e-3;
    if ( abs( m1 + m2 - m ) > epsil )
      {
	return 0.0;
      }

    double jmin = std::abs ( j1 - j2 );
    double jmax = j1 + j2;
    if ( ( j < jmin ) || ( j > jmax ) )
      {
	return 0.0;
      }

    double a     = j1;
    double alpha = m1;
    double b     = j2;
    double beta  = m2;
    double c     = j;
    double gamma = m;

    double f1    = factorial10( (int)( a + b - c ) );
    double f2    = factorial10( (int)( a + c - b ) );
    double f3    = factorial10( (int)( b + c - a ) );
    double f4    = factorial10( (int)( a + b + c + 1 ) );
    double Delta = std::sqrt( f1 * f2 * f3 / f4 );

    double p1    = std::sqrt( 2.0 * c + 1.0 );
  
    double f5    = factorial10( (int)( a + alpha ) );
    double f6    = factorial10( (int)( a - alpha ) );
    double f7    = factorial10( (int)( b + beta ) );
    double f8    = factorial10( (int)( b - beta ) );
    double f9    = factorial10( (int)( c + gamma ) );
    double f10   = factorial10( (int)( c - gamma ) );
  
    double p2    = std::sqrt( f5 * f6 * f7 * f8 * f9 * f10 );

    int numin1   = (int)( b - c - alpha );
    int numin2   = (int)( -c + a + beta );
    int numin3   = 0;
    int numin    = std::max( { numin1, numin2, numin3 } );
  
    int numax1   = (int)( a - alpha );
    int numax2   = (int)( b + beta );
    int numax3   = (int)( a + b - c );
    int numax    = std::min( { numax1, numax2, numax3 } );

    int nuphase  = 0;
    double f11   = 0.0;
    double f12   = 0.0;
    double f13   = 0.0;
    double f14   = 0.0;
    double f15   = 0.0;
    double f16   = 0.0;
    double term  = 0.0;
    double denom = 0.0;
    double sumnu = 0.0;
    for ( int nu = numin; nu <= numax; nu++ )
      {
	nuphase  = std::pow( -1, nu );
	f11      = factorial10( (int)( a - alpha ) - nu );
	f12      = factorial10( (int)( c - b + alpha ) + nu );
	f13      = factorial10( (int)( b + beta ) - nu );
	f14      = factorial10( (int)( c - a - beta ) + nu );
	f15      = factorial10( nu );
	f16      = factorial10( (int)( a + b - c ) - nu );
	denom    = f11 * f12 * f13 * f14 * f15 * f16;
	term     = nuphase / denom;
	sumnu    = sumnu + term;
      }
    return Delta / std::sqrt( 10.0 ) * p1 * p2 * sumnu;
  }
  /****************************************************************************
   Function    :: factorial10

   Input       :: int n 

   Output      :: double 

   Description :: Returns the factorial of argument 'n', divided by 10^n
  ****************************************************************************/
  double factorial10 ( int n )
  {
    if ( n == 0 )
      {
	return 1.0;
      }
    if ( n != 0 )
      {
	double tmp = 1.0;
	double q   = 1.0;
	for ( int i = 1; i<= n; i++ )
	  {
	    tmp = tmp * q / 10.0;
	    q   = q + 1.0;
	  }
	return tmp;
      }
  }
  /****************************************************************************
   Function    :: wigner_d_matrix

   Input       :: double j
                  double m
                  double n
                  double beta

   Output      :: double 

   Description :: Returns the Wigner little-d-matrix  d^{j}_{mn}(beta),
                  where 'beta' is the polar angle
  ****************************************************************************/
  double wigner_d_matrix ( double j,
			   double m,
			   double n,
			   double beta )
  {
    if ( ( abs(m) > j ) || ( abs(n) > j ) )
      {
	//	std::cerr << "input error: |m| or |n| larger than j" << std::endl;
        return 0.0;
      }

    //    double cos2b    = std::sqrt ( ( 1.0 + z ) / 2 );
    //    double sin2b    = std::sqrt ( ( 1.0 - z ) / 2 );
    double cos2b    = std::cos ( beta / 2.0 );
    double sin2b    = std::sin ( beta / 2.0 );

    int itmin1      = 0;
    int itmin2      = (int) ( m - n );
    int itmin       = std::max ( itmin1, itmin2 );
  
    int itmax1      = (int) ( j + m );
    int itmax2      = (int) ( j - n );
    int itmax       = std::min ( itmax1, itmax2 );
  
    int ij1         = (int) ( j - m );
    int ij2         = (int) ( j + n );
    double sqrt_fac = std::sqrt ( factorial10 (itmax1) * 
				  factorial10 (ij1) * 
				  factorial10 (ij2) * 
				  factorial10 (itmax2) );

    int iphase      = 0;
    int ia          = 0;
    int ib          = 0;
    int ic          = 0;
    double sumt     = 0.0;
    double denom    = 0.0;
    double term     = 0.0;

    for ( int it = itmin; it <= itmax; it++ )
      {
	iphase   = std::pow( -1, it );
	ia       = itmax1 - it;
	ib       = itmax2 - it;
	ic       = it + (int) (n - m);
	denom    = factorial10 (ia) * factorial10 (ib) 
	         * factorial10 (it) * factorial10 (ic);
	term     = iphase * std::pow ( cos2b, ia + ib ) 
	         * std::pow ( sin2b, it + ic ) / denom;
	sumt     = sumt + term;
      }
    return sqrt_fac * sumt;
  }
  /****************************************************************************
   Function    :: wigner_D_matrix

   Input       :: double j
                  double m
                  double n
                  double alpha
                  double beta
                  double gamma

   Output      :: complex<double>

   Description :: Returns the Wigner D-matrix  D^{j}_{mn}(alpha,beta,gamma),
                  where '(alpha,beta,gamma)' are the Euler angles
  ****************************************************************************/
  cd wigner_D_matrix (double j,
		      double m,
		      double n,
		      double alpha,
		      double beta,
		      double gamma )
  {
    double djmn = wigner_d_matrix ( j, m, n, beta );
    return std::exp( - ii * m * alpha ) * djmn * std::exp ( - ii * n * gamma );
  } 
  /****************************************************************************
   Function    :: spherical_harmonics

   Input       :: int l
                  int m
                  double beta
                  double alpha

   Output      :: complex<double>

   Description :: Returns the spherical harmonics Y_{lm}(beta,alpha)
                  where 'beta' is the polar angle, 'alpha' is azimuthal angle,
                  'l' is the angular momentum, and 'm' is the projection
  ****************************************************************************/
 cd spherical_harmonics ( int l, 
			  int m, 
			  double beta, 
			  double alpha )
 {
   double L    = (double) l;
   double M    = (double) m;
   double zero = 0.0;
   cd Dlm = wigner_D_matrix ( L, M, zero, alpha, beta, zero );
   return std::sqrt( ( 2.0 * L + 1.0 ) / 4.0 / pi ) * std::conj( Dlm );
 }

 /****************************************************************************
   Function :: dot_product

   Returns the Euclidean scalar product of two vectors "A" and "B".
  ****************************************************************************/
 double dot_product ( vd &A,
		  vd &B )
 {
   if ( A.size() != B.size() )
     {
       exit(-1);
     }
   double dot = 0.0;
   for ( int i = 0; i < A.size(); i++ )
     {
       dot = dot + A[i] * B[i];
     }
   return dot;
 }

 /****************************************************************************
   Function :: magnitude

   Returns the magnitude of a given vector "A"
  ****************************************************************************/
 double magnitude ( vd &A )
 {
   return std::sqrt ( dot_product ( A, A ) );
 }
  /****************************************************************************
   Function    :: cartesian_spherical_harmonics

   Input       :: vector<double> rvec
                  int l
                  int m

   Output      :: complex<double>

   Description :: Returns the cartesian spherical harmonics r^l * Y_{lm}(r)
                  where 'r' is a 3-vector, r = (rx, ry, rz), 
                  'l' is the angular momentum, and 'm' is the projection
  ****************************************************************************/
  cd cartesian_spherical_harmonics ( vd &rvec,
				     int l,
				     int m )
  {
    if ( ( abs(m) > l ) )
      {
	std::cerr << "input error: |m| larger than l" << std::endl;
        return 0;
      }
    double r = magnitude ( rvec );
    double x = rvec[0];
    double y = rvec[1];
    double z = rvec[2];
    cd rYlm(0.0,0.0);

 /************
  Case l = 0
  ************/
    if ( l == 0 )
      {
	if ( m == 0 )
	  {
	    rYlm =  std::sqrt ( 1.0 / 4.0 / pi );
	  }
      }
 /************
  Case l = 1
  ************/
    else if ( l == 1 )
      {
	if ( m == -1 )
	  {
	    rYlm = (1.0/2.0) * std::sqrt ( 3.0 / ( 2.0 * pi ) ) * (x - ii * y);
	  }
	else if ( m == 0 )
	  {
	    rYlm = z * std::sqrt ( 3.0 / pi ) / 2.0;
	  }
	else if ( m == 1 )
	  {
	    rYlm = -(1.0/2.0) * std::sqrt ( 3.0 / (2.0 * pi) ) * (x + ii * y);
	  }
      }
 /************
  Case l = 2
  ************/
    else if ( l == 2 )
      {
	if ( m == -2 )
          {
            rYlm = (1.0/4.0) * std::sqrt ( 15.0 / (2.0 * pi) )
              * std::pow ( x - ii * y, 2 );
          }
	else if ( m == -1 )
          {
            rYlm = (1.0/2.0) * std::sqrt ( 15.0 / ( 2.0 * pi ) ) 
	      * z * ( x - ii * y );
          }
        else if ( m == 0 )
          {
            rYlm = (1.0/4.0) * std::sqrt ( 5.0 / pi ) 
	      * ( 2.0 * z * z - x * x - y * y );
          }
        else if ( m == 1 )
          {
            rYlm = -(1.0/2.0) * std::sqrt ( 15.0 / (2.0 * pi) ) 
	      * z * ( x + ii * y );
          }
	else if ( m == 2 )
          {
            rYlm = (1.0/4.0) * std::sqrt ( 15.0 / (2.0 * pi) )
              * std::pow ( x + ii * y, 2 );
          }
      }
 /************
  Case l = 3
  ************/
    else if ( l == 3 )
      {
	if ( m == -3 )
          {
            rYlm = (1.0 / 8.0) * std::sqrt ( 35.0 / pi )
              * std::pow ( x - ii * y, 3 );
          }
	else if ( m == -2 )
          {
            rYlm = (1.0/4.0) * std::sqrt ( 105.0 / (2.0 * pi) )
              * z * std::pow ( x - ii * y, 2 );
          }
        else if ( m == -1 )
          {
            rYlm = (1.0/8.0) * std::sqrt ( 21.0 / pi )
              * ( x - ii * y ) * ( 4.0 * z * z - x * x - y * y );
          }
        else if ( m == 0 )
          {
            rYlm = (1.0/4.0) * std::sqrt ( 7.0 / pi )
              * z * ( 2.0 * z * z - 3.0 * x * x - 3.0 * y * y );
          }
        else if ( m == 1 )
          {
            rYlm = -(1.0/8.0) * std::sqrt ( 21.0 / pi )
              * (x + ii * y) * ( 4.0 * z * z - x * x - y * y);
          }
        else if ( m == 2 )
          {
            rYlm = (1.0/4.0) * std::sqrt ( 105.0 / (2.0 * pi) )
              * z * std::pow ( x + ii * y, 2 );
          }
	else if ( m == 3 )
          {
            rYlm = -(1.0/8.0) * std::sqrt ( 35.0 / pi )
              * std::pow ( x + ii * y, 3 );
          }
      }
 /************
  Case l = 4
  ************/
    else if ( l == 4 )
      {
	if ( m == -4 )
          {
            rYlm = (3.0 / 16.0) * std::sqrt ( 35.0 / ( 2.0 * pi ) )
              * std::pow ( x - ii * y, 4 );
          }
        else if ( m == -3 )
          {
            rYlm = (3.0 / 8.0) * std::sqrt ( 35.0 / pi )
              * z * std::pow ( x - ii * y, 3 );
          }
        else if ( m == -2 )
          {
            rYlm = (3.0/8.0) * std::sqrt ( 5.0 / (2.0 * pi) )
              * std::pow ( x - ii * y, 2 ) * ( 7.0 * z * z - r * r );
          }
        else if ( m == -1 )
          {
            rYlm = (3.0/8.0) * std::sqrt ( 5.0 / pi )
              * z * ( x - ii * y ) * ( 7.0 * z * z - 3.0 * r * r );
          }
        else if ( m == 0 )
          {
            rYlm = (3.0/16.0) * std::sqrt ( 1.0 / pi )
              * ( 35.0 * std::pow (z, 4) - 30.0 * z * z * r * r 
		  + 3.0 * std::pow (r, 4) );
          }
        else if ( m == 1 )
          {
            rYlm = -(3.0/8.0) * std::sqrt ( 5.0 / pi )
              * z * ( x + ii * y ) * ( 7.0 * z * z - 3.0 * r * r );
          }
        else if ( m == 2 )
          {
            rYlm = (3.0/8.0) * std::sqrt ( 5.0 / (2.0 * pi) )
              * std::pow ( x + ii * y, 2 ) * ( 7.0 * z * z - r * r );
          }
        else if ( m == 3 )
          {
            rYlm = -(3.0/8.0) * std::sqrt ( 35.0 / pi )
              * z * std::pow ( x + ii * y, 3 );
          }
	else if ( m == 4 )
          {
            rYlm = (3.0/16.0) * std::sqrt ( 35.0 / ( 2.0 * pi ) )
              * std::pow ( x + ii * y, 4 );
          }
      }
    /************
     End of Cases
     ************/
    return rYlm;
  }



  void print_clebsch_gordan ()
  {

    double cg = 0.0;

    for ( int ij1 = 1; ij1 <= 2; ij1++ )
      {
	for ( int ij2 = 1; ij2 <= 2; ij2++ )
	  {
	    double j1    =   ( (double) ij1 ) / 2.0;
	    double j2    =   ( (double) ij2 ) / 2.0;
	    double m1min = - j1;
	    double m1max =   j1;
	    double m2min = - j2;
	    double m2max =   j2;
	    double jmin  =   std::abs ( j1 - j2 );
	    double jmax  =   j1 + j2;
	    double mmin  = - jmax;
	    double mmax  =   jmax;

	    double m1 = m1min;
	    std::cout << " j1 = " << j1 << ", j2 = " << j2 << std::endl;
	    for ( int im1 = 0; im1 <= ( 2 * (int) j1 + 1 ); im1++ )
	      {
		double m2 = m2min;
		for ( int im2 = 0; im2 <= ( 2 * (int) j2 + 1 ); im2++ )
		  {
		    double j = jmin;
		    for ( int ij = 0; ij <= 2; ij++ )
		      {
			double m = -j;
			for ( int im = 0; im < ( 2 * (int) j + 1 ); im++ )
			  {
			    cg = clebsch_gordan ( j1, m1, j2, m2, j, m );
			    if ( cg != 0.0 )
			      {
			    std::cout << " m1 = " << m1 << ", m2 = " << 
				      m2 << ", j =  " << j << ", m =  " 
				      << m << std::endl;
			    std::cout << " " << std::endl;
			    std::cout << cg << std::endl;
			      }
			    m = m + 1.0;
			  }
			j = j + 1.0;
			//m = -j;
		      }
		    m2 = m2 + 1.0;
		    j  = jmin;
		  }
		m1 = m1 + 1.0;
		m2 = m2min;
	      }
	    
	  }
      }
    

  }

  void plot_legendre ( int n )
  {
    int Npts    = 100;
    double xi   = -1.0;
    double xf   =  1.0;
    double step = ( xf - xi ) / ( (double) Npts );
    double x    = xi;
    for ( int i = 0; i <= Npts; i++ )
      {
	double leg   = legendre_polynomial ( n, x );
	double leg_d = legendre_derivative ( n, x );
	std::cout << x << " " << leg << " " << leg_d << std::endl;
	x = x + step;
      }
    return;
  }


  void plot_wigner_d ( double J,
		       double m_p,
		       double m )
  {
    double theta_i = 0.0;
    double theta_f = pi;
    int Npts = 1000;
    double step = ( theta_f - theta_i ) / ( (double) Npts );
    double wd = 0.0;
    double theta = theta_i;
    for ( int i = 0; i <= Npts; i++ )
      {
	wd = wigner_d_matrix ( J, m_p, m, theta );
	std::cout << theta << " " << wd << std::endl;
	theta = theta + step;
      }
    return;
  }

  void plot_spherical_harmonics ( int l,
				  int m )
  {
    int Npts       =   500;
    double theta_i = - pi;
    double theta_f =   pi;
    double phi     =   0.0;
    double step    =  ( theta_f - theta_i ) / ( (double) Npts );
    cd Ylm         =   iz;
    double theta   =   theta_i;
    for ( int i = 0; i <= Npts; i++ )
      {
	Ylm    = spherical_harmonics ( l, m, theta, phi );
	std::cout << theta << " " << std::real ( Ylm ) << std::endl;
	theta = theta + step;
      }
    return;
  }

  void plot2d_spherical_harmonics ( int l,
				    int m )
  {
    int Npts       =  150;
    double theta_i = -pi;
    double theta_f =  pi;
    double phi_i   =  0.0;
    double phi_f   =  2.0 * pi;
    double step_t  =  ( theta_f - theta_i ) / ( (double) Npts );
    double step_p  =  ( phi_f - phi_i ) / ( (double) Npts );
    double theta   =  theta_i;
    double phi     =  phi_i;
    cd Ylm         =  iz;
    double Ylm_sq  =  0.0;
    for ( int i = 0; i <= Npts; i++ )
      {
	for ( int j = 0; j <= Npts; j++ )
	  {
	    Ylm    = spherical_harmonics ( l, m, theta, phi );
	    Ylm_sq = std::real ( std::conj(Ylm) * Ylm );
	    std::cout << phi << " " << theta << " " 
		      << Ylm_sq << " " 
		      << std::real ( Ylm ) << " "
		      << std::imag ( Ylm ) << std::endl;
	    phi = phi + step_p;
	  }
	std::cout << " " << std::endl;
	theta = theta + step_t;
	phi   = phi_i;
      }
  }

  /****************************************************************************
  End of Namespace :: rotational_functions
  ****************************************************************************/
}
#endif
