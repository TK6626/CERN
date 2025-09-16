#include "computations.h"
#include "TMath.h"

float calc_inv_mass(const array_4 &p) {
	return p[0]*p[0] - p[1]*p[1] - p[2]*p[2] - p[3]*p[3];
}

array_4 add_4_vec(const array_4 &p1, const array_4 &p2){
	array_4 p;
	for (int i=0;i<4;++i){
		p[i] = p1[i]+p2[i];
		};
	return p;
	}

/**
 * @brief Gaussian fit function with offset.
 *
 * This function defines a Gaussian curve with an added constant offset.
 * It is typically used in ROOT for fitting histograms or graphs.
 *
 * @param x Pointer to the independent variable (usually a single-element array).
 * @param params Array of parameters:
 *		- params[0]: Mean (center of the peak)
 *		- params[1]: Standard deviation (width of the peak)
 *		- params[2]: Amplitude (height of the peak)
 *		- params[3]: Constant offset (baseline)
 *
 * @return Value of the Gaussian function at point x[0].
 *
 *
 * Note that the FWHM is not standard deviation but 2 * sqrt(2 ln(2)) * std (~2.355)
 */
 Double_t fit_gaussian(Double_t *x, Double_t *params)
{  
	Double_t mu = params[0];
	Double_t sig = params[1];
	Double_t A = params[2];
	Double_t C = params[3];


	Double_t fit =   C + A*  TMath::Exp(-TMath::Power((x[0] - mu) ,2) / (2 * TMath::Power(sig, 2)));
	return fit;
};

/**
 * @brief Exponential of a quartic polynomial with constant offset.
 *
 * This function defines: f(x) = C + exp(a·x⁴ + b·x³ + c·x² + d·x + e)
 * It's useful for fitting sharp or asymmetric peaks with heavy tails.
 *
 * @param x Pointer to the independent variable (usually a single-element array).
 * @param params Array of parameters:
 *		- params[0]: a (coefficient of x⁴)
 *		- params[1]: b (coefficient of x³)
 *		- params[2]: c (coefficient of x²)
 *		- params[3]: d (coefficient of x)
 *		- params[4]: e (constant term in exponent)
 *		- params[5]: C (constant offset added after the exponential)
 *
 * @return Value of the exponential-quartic function at point x[0].
 */
Double_t fit_exp_quartic(Double_t *x, Double_t *params)
{
	Double_t a = params[0];
	Double_t b = params[1];
	Double_t c = params[2];
	Double_t d = params[3];
	Double_t A = params[4];
	Double_t C = params[5];

	Double_t xx = x[0];
	Double_t exponent = a*TMath::Power(xx, 4) + b*TMath::Power(xx, 3) +
						c*TMath::Power(xx, 2) + d*xx;

	Double_t fit = C + A *TMath::Exp(exponent);
	return fit;
};



/** 
 * @breif Breit-Wigner fit with offset
 *
 *
 * Defines a relativistic breit wigner fit with (constant)
 * background offset 
 * As based on wiki formula - so sue me. Will find a text book
 * or something that is more official
 *
 */
Double_t fit_breit_wigner(Double_t *E, Double_t *par)
{	
	// Params
	Double_t m = par[0]; 	// mass of resonance - MeV/c^2
	Double_t G = par[1]; 	// FWHM of resonance - MeV/c^2
	Double_t A = par[2];	// Scaling factor - Unitless
	Double_t C = par[3]; 	// Background/ constant offset - MeV

	Double_t f = (m*m * G*G) / (TMath::Power(E[0]*E[0] - m*m, 2) + m*m * G*G);
	return A*f + C;
}

Double_t fit_relativistic_breit_wigner(Double_t *E, Double_t *par)
{
	// Params
	Double_t m = par[0]; 	// mass of resonance - MeV/c^2
	Double_t G = par[1]; 	// FWHM of resonance - MeV/c^2
	Double_t A = par[2];	// Scaling factor - Unitless
	Double_t C = par[3]; 	// /Constant offset - MeV


	// True constant k given in the numerator
	Double_t g = TMath::Power(m*m * (m*m + G*G), 0.5);
	Double_t k = (TMath::Power(2, 1.5) * m * G * g) / (TMath::Pi() * TMath::Power(m*m + g, 2));
	Double_t f = (k) / (TMath::Power(E[0]*E[0] - m*m, 2) + m*m * G*G);

	// Simplified version such that at m=E f = 1	

	return A*f + C;
};

Double_t fit_background(Double_t *E, Double_t *par)
{
	Double_t A_bg = par[0];
	Double_t B	= par[1];
	Double_t C	= par[2];
	Double_t D	= par[3];

	Double_t E0 = E[0];

	if (E0 > B) {
		return A_bg * TMath::Power(E0 - B, C) * TMath::Exp(D * E0);
	} else {
		return 0;
	}
}


Double_t fit_breit_wigner_bg(Double_t *E, Double_t *par)
{
	// First 4 parameters for Breit-Wigner
	Double_t bw = fit_relativistic_breit_wigner(E, par);

	// Next 4 parameters for background: shift pointer by 4
	Double_t bg = fit_background(E, &par[4]);

	return bw + bg;
}

Double_t fit_double_gaussian(Double_t *x, Double_t *par) {
	
	// use already existing gaussian
	Double_t gaus_1 = fit_gaussian(x, par);
	Double_t gaus_2 = fit_gaussian(x, &par[4]);

	return gaus_1 + gaus_2;
}

/** fills in a histogram taking from the invariant mass of the 4 pion events while also applying the rho mass cut to try and single out rho-rho events
 */
void fillHistFromP4(const RVecLor& p4, TH1* hist, Float_t m_rho, Float_t mass_bound) {
	auto m0 = p4[0].M() * 1e3;
	auto m1 = p4[1].M() * 1e3;
	auto m2 = p4[2].M() * 1e3;
	auto m3 = p4[3].M() * 1e3;

	bool cond1 = TMath::Abs(m0 - m_rho) < mass_bound && TMath::Abs(m1 - m_rho) < mass_bound;
	bool cond2 = TMath::Abs(m2 - m_rho) < mass_bound && TMath::Abs(m3 - m_rho) < mass_bound;

	if (cond1) hist->Fill((p4[0] + p4[1]).M() * 1e3);
	if (cond2) hist->Fill((p4[2] + p4[3]).M() * 1e3);
}


RooGenericPdf* ShiftedArgusPdf(const char* name, const char* title,
    RooRealVar& x_var, RooRealVar& x_var0,
    RooRealVar& c, RooRealVar& p, RooRealVar& shift)
{
    TString expr =
        "((@0-@4) > 0 && (@0-@4) < @1) ? "
        "((@0-@4) * pow(1 - pow((@0-@4)/@1, 2), @3) * exp(@2*(1 - pow((@0-@4)/@1, 2)))) "
        ": 1e-30";

	return new RooGenericPdf(name, title, expr, RooArgList(x_var, x_var0, c, p, shift));
}

/*
RooGenericPdf* BiGaussianPdf(const char* name, const char* title,
        RooRealVar& x, RooRealVar& mu1, RooRealVar& sigma1,
        RooRealVar& mu2, RooRealVar& sigma2, RooRealVar& frac)
{
	TString expr = "@5*exp(-0.5*pow((@0-@1)/@2,2))+(1-@5)*exp(-0.5*pow((@0-@3)/@4,2))";
*/

RooGenericPdf* CubicExpPdf(const char* name, const char* title,
        RooRealVar& x, RooRealVar& a, RooRealVar& b,
        RooRealVar& c)
{
    TString expr = "exp(@1*pow(@0 - @2,3) + @2*pow(@0 - @2,2) + @3*(@0 - @2))";
	return new RooGenericPdf(name, title, expr, RooArgList(x, a, b, c));
}

RooGenericPdf* BivariateGaussianPdf(const char* name, const char* title,
    RooRealVar& x,   RooRealVar& mux,  RooRealVar& sigx,
    RooRealVar& y,   RooRealVar& muy,  RooRealVar& sigy,
    RooRealVar& rho)
{
    TString formula =
        "exp(-0.5/(1-@6*@6)*(((@0-@1)*(@0-@1))/(@2*@2) + "
        "((@3-@4)*(@3-@4))/(@5*@5) - 2*@6*(@0-@1)*(@3-@4)/(@2*@5)))";

    return new RooGenericPdf(name, title, formula,
        RooArgList(x, mux, sigx, y, muy, sigy, rho));
}


RooGenericPdf* ThresholdBackgroundPdf(
    const char* name, const char* title,
    RooRealVar& x_var, RooRealVar& x_var0,
    RooRealVar& p, RooRealVar& a, RooRealVar& b, RooRealVar& c)
{
	TString expr = "(@0 > @1) ?  pow((@0 - @1), @2) * exp(@3*pow(@0 -@1,3) + @4*pow(@0 -@1,2) + @5*(@0 - @1)): 0.0";
	return new RooGenericPdf(name, title, expr, RooArgList(x_var, x_var0, p, a, b, c));

}

RooGenericPdf* LogThresholdBackgroundPdf(
    const char* name, const char* title,
    RooRealVar& x_var, RooRealVar& x_var0,
    RooRealVar& p, RooRealVar& a, RooRealVar& b, RooRealVar& c, RooRealVar& C)
{
    TString expr = "(@0 > @1) ? @2 * ln(@0 - @1) + @3*pow(@0,3) + @4*pow(@0,2) + @5*@0 + @6: 0.0";
    return new RooGenericPdf(name, title, expr, RooArgList(x_var, x_var0, p, a, b, c));
}

