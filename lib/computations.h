#pragma once
#include "custom_definitions.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"


float calc_inv_mass(const array_4 &p);
array_4 add_4_vec(const array_4 &p1, const array_4 &p2);
Double_t fit_gaussian(Double_t *x, Double_t *params);
Double_t fit_breit_wigner(Double_t *E, Double_t *par);
Double_t fit_relativistic_breit_wigner(Double_t *E, Double_t *par);
Double_t fit_exp_quartic(Double_t *x, Double_t *params);
Double_t fit_breit_wigner_bg(Double_t *E, Double_t *par);
Double_t fit_double_gaussian(Double_t *x, Double_t *par);
void fillHistFromP4(const RVecLor& p4, TH1* hist, Float_t m_rho, Float_t mass_bound);


RooGenericPdf* ShiftedArgusPdf(const char* name, const char* title,
    RooRealVar& x_var, RooRealVar& x_var0,
    RooRealVar& c, RooRealVar& p, RooRealVar& shift);

RooAddPdf* BiGaussianPdf(const char* name, const char* title,
        RooRealVar& x, RooRealVar& mu1, RooRealVar& sigma1,
        RooRealVar& mu2, RooRealVar& sigma2, RooRealVar& frac);

RooGenericPdf* ThresholdBackgroundPdf(
    const char* name, const char* title,
    RooRealVar& x_var, RooRealVar& x_var0,
    RooRealVar& p, RooRealVar& a, RooRealVar& b, RooRealVar& c);

RooGenericPdf* LogThresholdBackgroundPdf(
    const char* name, const char* title,
    RooRealVar& x_var, RooRealVar& x_var0,
    RooRealVar& p, RooRealVar& a, RooRealVar& b, RooRealVar& c, RooRealVar& C);

RooGenericPdf* CubicExpPdf(const char* name, const char* title,
        RooRealVar& x, RooRealVar& a, RooRealVar& b,
        RooRealVar& c);

















