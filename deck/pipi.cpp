
// Phase-shift and Inelasticities for PiPi scattering based off GKPRY parameterization.
//
// Dependencies: None
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "pipi.hpp"

//-----------------------------------------------------------------------------
const double pipi::conv = (M_PI / 180.);

//Masses
const double pipi::mPi = 0.13957018;
const double pipi::mK = 0.496;
const double pipi::mEta = 0.54753;

const double pipi::mRho = .77545;
const double pipi::mF2 = 1.2754;

//Thresholds for pi, eta, and K
const double pipi::sthPi = 4.*mPi*mPi;
const double pipi::sthK = 4.*mK*mK;
const double pipi::sthEta = 4.*mEta*mEta;

//Unit imaginary and real
const std::complex<double> pipi::xr(1., 0.);
const std::complex<double> pipi::xi(0., 1.);

//-----------------------------------------------------------------------------
//Constructor
pipi::pipi()
{
        std::cout << " pipi : Input valid Isospin. \n";
        std::cout << " pipi : Usage: 'pipi amp(1); \n";
        exit(1);
}
pipi::pipi(int i)
{
        if (i < 0 || i > 2)
        {
                std::cout << " pipi : Amplitude with Invalid Isospin Initialized. \n";
                exit(1);
        }
        qn_I = i;
}

//-----------------------------------------------------------------------------
void pipi::print_I()
{
        std::cerr << " Isospin = " << qn_I << "\n";
};

//-----------------------------------------------------------------------------
//Returns conformal variable given Mandelstam s and branching point s0.
double pipi::conformal(double s, double s0)
{
        double numerator = sqrt(s) - sqrt(s0 - s);
        double denominator = sqrt(s) + sqrt(s0 - s);
        return numerator/denominator;
}
//-----------------------------------------------------------------------------
// Elastic momentum above threshold sth, as a function of s.
double pipi::elastic_mom( double s, double sth)
{
        return sqrt(s - sth) / 2.;
}
//-----------------------------------------------------------------------------
// Produces the Phase-shift for given Partial Wave as function of s.
// Based on parameterizations from CFD Set:
// R. Garcia-Martin, R. Kaminski, J. R. Palaez, J. Ruiz de Elvira, F. J. Yndurain
// Phys. Rev. D 83, 074004 (2011).
//
// Ported from A. Jackura and N. Sherrill's fortran code
// ajackura@indiana.edu, nlsherri@indiana.edu
double pipi::phase_shift(int l, double s)
{
        //Error Checks :)
        int wave = (3*l - qn_I)/2;
        if (qn_I > 2 || qn_I < 0) {wave = 600;}
        if (l % 2 != qn_I % 2) {wave = 600;}
        if (l >= 4 && wave != 600)
        {std::cout << "Phase shift contributions of l > 3 are completely negligible. Quitting... \n";
         exit(1);}

        double cot_delta, delta, sh;
        double temp1, temp2, temp3;
        sh = pow(1.42, 2.);

        //Momenta
        double k, k2;
        k = elastic_mom(s, pipi::sthPi);
        k2 = elastic_mom(s, pipi::sthK);

//S2 wave (l = 0, iso = 2)
        if (wave == -1)
        {
                double B0, B1, z2, wl;
                double sl, sm;

                sl = pow(1.05, 2.);
                sm = pow(.850, 2.);

                wl = conformal(s, sl);

                z2 = .1435;
                B0 = -79.4;
                B1 = -63.0;

                temp1 = sqrt(s) / (2.*k);
                temp2 = mPi*mPi / (s - 2.*z2*z2);

                if ((s > sthPi) && (s <= sm)) //Low energy parameterization
                {
                        temp3 = B0 + B1 * wl;
                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else if ((s > sm) && ( s < sh ) ) //Intermediate energies
                {
                        double Bh0, Bh1, Bh2, wh, whsm, wlsm;
                        double temp4;

                        wh = conformal(s, sh);
                        whsm = conformal(sm, sh);
                        wlsm = conformal(sm, sl);

                        Bh2 = 32.;
                        Bh0 = B0 + B1 * wlsm;
                        temp4 = pow((sqrt(sm) + sqrt(sh - sm)) / (sqrt(sm) + sqrt(sl - sm)), 2.);
                        Bh1 = B1 * (sl / sh) * (sqrt(sh - sm) / sqrt(sl - sm)) * temp4;

                        temp3 = Bh0 + Bh1*(wh - whsm) + Bh2*pow(wh - whsm, 2.);

                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}

                return delta;
        }

//S0 wave (l = 0, iso = 0)
        if (wave == 0)
        {
                double sm = pow(0.85, 2.);
                double B0, B1, B2, B3, w;
                B0 = 7.14;
                B1 = -25.3;
                B2 = -33.2;
                B3 = -26.2;
                w = conformal(s, sthK);

                if ((s > sthPi) && (s <= sm)) //Low energy parameterization
                {
                        temp1 = sqrt(s) / (2.*k);
                        temp2 = mPi * mPi / (s - 0.5 * mPi * mPi);
                        temp3 = mPi / sqrt(s);

                        cot_delta = temp1 * temp2 * (temp3 + B0 + B1 * w + B2 * w * w + B3 *w*w*w);
                        delta = atan2(1., cot_delta);
                }
                else if (s < sh )
                {
                        double d0, C1, B, C2, D;
                        double k2, k2m;

                        d0 = 226.5 * pipi::conv; //converting to radians
                        C1 = -81. * pipi::conv;
                        B = 93.3 * pipi::conv;
                        C2 = 48.7 * pipi::conv;
                        D = -88.3 * pipi::conv;


                        k2m = elastic_mom(sthK, sm); //Switched inputs to avoid negative under radical

                        if ((s > sm) && (s < sthK)) // intermediate energy
                        {
                                double cot_delm, delm, delPm, km, wm;
                                double temp4, temp5, temp6;
                                k2 = elastic_mom(sthK, s);
                                temp4 = pow(1. - (k2 / k2m), 2.);
                                temp5 = k2 * (2. - (k2 / k2m)) / k2m;
                                temp6 = (k2m - k2) / pow(mK, 3.);


                                //Calculate delta(sm)
                                wm = conformal(sm, sthK);
                                km = elastic_mom(sm, sthPi);
                                temp1 = sqrt(sm) / (2. * km);
                                temp2 = mPi * mPi / (sm - 0.5 * mPi * mPi);
                                temp3 = mPi / sqrt(sm);
                                cot_delm = temp1 * temp2 * (temp3 + B0 + B1 * wm + B2 * wm * wm + B3 *wm*wm*wm);
                                delm = atan2(1., cot_delm);

                                //Derivative calculated in Mathematica (because im lazy)
                                delPm = 1.588503;

                                delta = d0 * temp4 + delm * temp5 + k2 * (k2m - k2) * (8.*delPm + C1 * temp6);
                        }

			if ((s > sthK) && (s < sh)) //above KK threshold
                        {
                                k2 = elastic_mom(s, sthK);
                                temp1 = (k2*k2) / (mK * mK);
                                delta = d0 + B*temp1 + C2 * temp1*temp1;

                                if ( s > sthEta )
                                {
                                        double k3 = elastic_mom(s, sthEta);
                                        temp2 = D * (k3*k3) / (mEta * mEta);
                                        delta += temp2;
                                }
			}
						
                }
		if ( s > sh )
		  {
		    delta = 0.0;
		  }
                return delta;
        }

//P1 wave (l = 1, iso = 1)
        else if (wave == 1)
        {
                double s0, w;

                s0 = pow(1.05, 2.);
                w = conformal(s, s0);

                if ((s > sthPi) && (s < sthK))
                {
                        double B0, B1;
                        B0 = 1.043;
                        B1 = .19;

                        temp1 = sqrt(s)*(mRho*mRho - s)/(2. * pow(k, 3.));
                        temp2 = 2.*pow(mPi, 3.)/(mRho*mRho*sqrt(s));

                        cot_delta = temp1*(temp2 + B0 + B1* w);
                        delta = atan2(1., cot_delta);
                }
                else if ((s >= sthK) && (s < sh))
                {
                        double lambda0, wK, KK;
                        double lambda1, lambda2;
                        lambda1 = 1.38;
                        lambda2 = -1.70;

                        //lambda0 = low energy (see above) at KK threshold
                        double B0, B1;
                        B0 = 1.043;
                        B1 = .19;
                        KK = elastic_mom(sthK, sthPi);
                        wK = conformal(sthK, s0);
                        temp1 = sqrt(sthK)*(mRho*mRho - s)/(2. * pow(KK, 3.));
                        temp2 = 2.*pow(mPi, 3.)/(mRho*mRho*sqrt(sthK));
                        lambda0 = temp1*(temp2 + B0 + B1* wK);

                        temp3 = (sqrt(s) / (2.*mK)) - 1.;

                        cot_delta = lambda0 + lambda1*temp3 + lambda2*temp3*temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}

                return delta;
        }

//D2 wave (l = 2, iso = 2)
        else if (wave == 2)
        {
	  if ((s > sthPi) && (s < sh))
                {
                        double w, s0, Del;
                        double B0, B1, B2;

                        s0 = pow(1.45,2.);
                        B0 = 4.1e3;
                        B1 = 8.6e3;
                        B2 = 25.5e3;
                        Del = .233;

                        w = conformal(s, s0);

                        temp1 = sqrt(2) / (2.* pow(k, 5.));
                        temp2 = (pow(mPi, 4.) * s) / (sthPi + 4.*Del*Del -s);
                        temp3 = B0 + B1 * w + B2 * w * w;

                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}
                return delta;

        }

//D0 wave (l = 2, iso = 0)
        else if (wave == 3)
        {
                double B0, B1;
                double s0 = pow(1.05, 2.);

                B0 = 12.40;
                B1 = 10.06;

                temp1 = sqrt(s) / (2.*pow(k, 5.));
                temp2 = (mF2*mF2 - s)*mPi*mPi;

                if ((s > sthPi) && (s <= sthK))
                {
                        double w = conformal(s, s0);
                        temp3 = B0 + B1* w;
                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else if ((s > sthK) && (s < sh))
                {
                        double wh, s02, B0h, B1h;

                        B1h = 43.2;
                        s02 = pow(1.45, 2.);
                        wh = conformal(s, s02);
                        B0h = 18.69;

                        temp3 = B0h + B1h * wh;
                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}
        }

//F1 wave (l = 3, iso = 1)
        else if (wave == 4)
        {
                if ((s > sthPi) && (s < sh))
                {
                        double B0, B1, lambda;
                        double s0, w;

                        s0 = pow(1.45, 2.);
                        w = conformal(s, s0);
                        B0 = 1.09e5;
                        B1 = 1.41e5;
                        lambda = 0.051e5;

                        temp1 = sqrt(s)*pow(mPi, 6.) / (2.*pow(k, 7.));
                        temp2 = 2.*lambda * mPi / sqrt(s);

                        cot_delta = temp1 * (temp2 + B0 + B1 * w);
                        delta = atan2(1., cot_delta);
                }
        }

//Unphysical partial waves (ie anything else)
        else
        {
                std::cout << "Invalid Phase Shift with l = " << l << " and isospin = " << qn_I << ". Quiting... \n";
                exit(1);

        }
}
//-----------------------------------------------------------------------------
// Produces the Inelasticity for given Partial Wave as a function of s.
// Based on CFD parameterizations from:
// R. Garcia-Martin, R. Kaminski, J. R. Palaez, J. Ruiz de Elvira, F. J. Yndurain
// Phys. Rev. D 83, 074004 (2011).
//
// Ported from A. Jackura and N. Sherrill's fortran code
// ajackura@indiana.edu, nlsherri@indiana.edu
double pipi::inelasticity(int l, double s)
{
        double eta;
        double sh = pow(1.42, 2.);

        int wave = (3*l - qn_I)/2;
        if (qn_I > 2 || qn_I < 0 || l < 0) {wave = 600;}
        if (l % 2 != qn_I % 2) {wave = 600;} //Check Bose symmetry
        if (l >= 3 && wave != 600) {return eta = 1.;} //Inelasticities for high partial waves are negligable

//S2 wave (l = 0, iso = 2)
        if (wave == -1)
        {
                double si = pow(1.05, 2.);

                if ((s > sthPi) && (s <= si))          //Low energy parameterization
                {
                        eta = 1;
                }
                else if ((s > si) && ( s <= sh))        //Intermediate energies
                {
                        double eps = .28;

                        eta = 1. - eps*pow(1 - (si/s), 1.5);
                }
                return eta;
        }

//S0 wave (l = 0, iso = 0)
        else if (wave == 0)
        {
                if ((s > sthPi) && (s <= sthK))
                {
                        eta = 1.;
                }
                else if ((s > sthK) && (s < sh))
                {
                        double ep1, ep2, ep3;
                        double k2;
                        double temp1, temp2, temp3;

                        ep1 = 4.9;
                        ep2 = -15.1;
                        ep3 = 4.7;
                        k2 = elastic_mom(s, sthK);

                        temp1 = -k2 / sqrt(s);
                        temp2 = ep1 + ep2 * (k2/sqrt(s)) + ep3 * (k2 * k2 / s);
                        temp3 = temp1 * temp2 * temp2;
                        if (s > sthEta)
                        {
                                double ep4, k3;
                                ep4 = 0.32;
                                k3 = elastic_mom(s, sthEta);

                                temp3 += -ep4 * k3 / sqrt(s);
                        }
                        eta = exp(temp3);
                }
                return eta;
        }

//P1 wave (l = 1, iso = 1)
        else if (wave == 1)
        {
                if ((s > sthPi) && (s < sthK))
                {
                        eta = 1.;
                }
                else if ((s > sthK) && (s < sh))
                {
                        double temp;
                        double ep1, ep2;
                        ep1 = 0.;
                        ep2 = .07;

                        temp = 1. - (sthK/s);
                        eta = 1 - ep1 * sqrt(temp) - ep2 * temp;
                }
                return eta;
        }

//D2 wave (l = 2, iso = 2)
        else if (wave == 2)
        {
                double si = pow(1.05, 2.);
                if ((s > sthPi) && (s <= si))
                {
                        eta = 1.;
                }
                else if ((s > si) && (s < sh))
                {
                        double ep = 0.;
                        eta = 1. - ep * pow(1 - (si/s),3.);
                }
                return eta;
        }

//D0 wave (l = 2, iso = 0)
        else if (wave == 3)
        {
                if ((s > sthPi) && (s <= sthK))
                {
                        eta = 1.;
                }
                if ((s > sthK) && (s < sh))
                {
                        double eps, r, temp1, temp2;
                        eps = 0.254;
                        r = 2.29;

                        double k2 = elastic_mom(s, sthK);

                        temp1 = (1. - sthK/s) / (1. - sthK/(mF2*mF2));
                        temp2 = 1. - k2 / elastic_mom(mF2*mF2, sthK);

                        eta = 1. - eps*pow(temp1, 2.5)*(1. + r*temp2);
                }
                return eta;
        }

//Unphysical partial waves (ie anything else)
        else
        {
                std::cout << "Invalid Inelasticity with l = " << l << " and isospin = " << qn_I << ". Quiting... \n";
                exit(1);
        }
}

//-----------------------------------------------------------------------------
std::complex<double> pipi::GKPRY_partial_wave(int l, double s)
{
        double delta, eta, k;
        std::complex<double> amp;

        k = elastic_mom(s, sthPi);
        delta = phase_shift(l, s);
        eta = inelasticity(l, s);

        amp = eta * exp(2.*xi*delta) - 1.;
        amp *= sqrt(s) / (4. * k * xi);

        return amp;
}
//-----------------------------------------------------------------------------
//Produces Isospin-definite Amplitude
std::complex<double> pipi::GKPRY_iso_amp(double s, double z)
{
        std::complex<double> amp, pw;
        double Pl;
        int l;

        for (int i = 0; i < 2; i++)
        {
                if (qn_I == 0 || qn_I == 2)
                {
                        l = 2*i; //only even spin
                }
                else if (qn_I == 1)
                {
                        l = 2*i + 1;
                }
                if (qn_I > 2 || qn_I < 0)
                {
                        std::cerr << " Invalid isospin! Quitting... \n";
                        exit(1);
                }


                Pl = pipi::legendre(l, z);
                pw = GKPRY_partial_wave(l, s);
                amp += (2.*double(l) + 1.) * Pl * pw;
        }
        amp *= 32. * M_PI;
        return amp;
}
//-----------------------------------------------------------------------------
// Legendre Polynomials P_l(x) (up to l = 5)
double pipi::legendre(int l, double x)
{
        double PL;
        switch (l)
        {
        case 0: PL = 1.; break;
        case 1: PL = x; break;
        case 2: PL = .5 * (3.*pow(x, 2.) - 1.); break;
        case 3: PL = .5 * (5.*pow(x,3.) - 3.*x); break;
        case 4: PL = (35.*pow(x,4.) - 30.*pow(x,2.) + 3.)/8.; break;
        case 5: PL = (63.*pow(x,5.) - 70.*pow(x,3.) + 15.*x)/8.; break;
        default: std::cerr << "Legendre Polynomials don't go that high!!! \n";
                exit(1);
        }
        return PL;
}
