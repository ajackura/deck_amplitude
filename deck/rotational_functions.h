#ifndef ROTATIONAL_FUNCTIONS_H_
#define ROTATIONAL_FUNCTIONS_H_

#include <iostream>
#include <complex>
#include <algorithm>
#include <cmath>
#include <vector>

typedef std::complex<double> cd;
typedef std::vector<double> vd;

namespace rotational_functions
{

  double legendre ( double x,
		    int n );

  double clebsch_gordan ( double j1,
			  double m1,
			  double j2,
			  double m2,
			  double j,
			  double m );

  double fac10 ( int n );

  double djmn ( double j,
		double m,
		double n,
		double z );



  /*

   */

  double legendre ( double x,
		    int n )
  {
 
    // This function computes the Legendre polynomial of degree N
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
  } // end of function Legendre 


  double clebsch_gordan( double j1,
			 double m1,
			 double j2,
			 double m2,
			 double j,
			 double m )
  {
    if ( ( abs(m1) > j1 ) || ( abs(m2) > j2 ) || ( abs(m) > j ) )
      {
	std::cout << "input error: |m| larger than j" << std::endl;
	return 0;
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

    double f1    = fac10( (int)( a + b - c ) );
    double f2    = fac10( (int)( a + c - b ) );
    double f3    = fac10( (int)( b + c - a ) );
    double f4    = fac10( (int)( a + b + c + 1 ) );
    double Delta = std::sqrt( f1 * f2 * f3 / f4 );

    double p1    = std::sqrt( 2.0 * c + 1.0 );
  
    double f5    = fac10( (int)( a + alpha ) );
    double f6    = fac10( (int)( a - alpha ) );
    double f7    = fac10( (int)( b + beta ) );
    double f8    = fac10( (int)( b - beta ) );
    double f9    = fac10( (int)( c + gamma ) );
    double f10   = fac10( (int)( c - gamma ) );
  
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
	f11      = fac10( (int)( a - alpha) - nu );
	f12      = fac10( (int)( c - b + alpha ) + nu );
	f13      = fac10( (int)( b + beta ) - nu );
	f14      = fac10( (int)( c - a - beta ) + nu );
	f15      = fac10( nu );
	f16      = fac10( (int)( a + b - c ) - nu );
	denom    = f11 * f12 * f13 * f14 * f15 * f16;
	term     = nuphase / denom;
	sumnu    = sumnu + term;
      }
    return Delta / std::sqrt( 10.0 ) * p1 * p2 * sumnu;
  }

  double fac10 ( int n )
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

  double djmn ( double j,
		double m,
		double n,
		double z )
  {
    double cos2b    = std::sqrt ( ( 1.0 + z ) / 2 );
    double sin2b    = std::sqrt ( ( 1.0 - z ) / 2 );

    int itmin1      = 0;
    int itmin2      = (int) ( m - n );
    int itmin       = std::max ( itmin1, itmin2 );
  
    int itmax1      = (int) ( j + m );
    int itmax2      = (int) ( j - n );
    int itmax       = std::min ( itmax1, itmax2 );
  
    int ij1         = (int) ( j - m );
    int ij2         = (int) ( j + n );
    double sqrt_fac = std::sqrt ( fac10 (itmax1) * 
				  fac10 (ij1) * 
				  fac10 (ij2) * 
				  fac10 (itmax2) );

    double sumt     = 0.0;
    int iphase      = 0;
    int ia          = 0;
    int ib          = 0;
    int ic          = 0;
    double denom    = 0.0;
    double term     = 0.0;

    for ( int it = itmin; it <= itmax; it++ )
      {
	iphase   = std::pow( -1, it );
	ia       = itmax1 - it;
	ib       = itmax2 - it;
	ic       = it + (int) (n - m);
	denom    = fac10 (ia) * fac10 (ib) * fac10 (it) * fac10 (ic);
	term     = iphase * std::pow ( cos2b, ia + ib ) 
	  * std::pow ( sin2b, it + ic ) / denom;
	sumt     = sumt + term;
      }
    return sqrt_fac * sumt;
  }


  
}

#endif
