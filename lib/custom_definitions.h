#pragma once

#include "TMath.h"
#include <array>
#include <stdio.h>
#include <iostream>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "TH1.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"


// common constants
const Float_t m_p = 938.27; // MeV/c^2
const Float_t m_rho = 766.5; // MeV/c^2
const Float_t m_kaon_0 = 497.611; // MeV/c^2
const Float_t m_kaon_char = 493.677; // MeV/c^2
const Float_t m_kaon_star = 891.66; // MeV/c^2
const Float_t m_phi = 1019.460; // MeV/c^2
const Float_t pi = TMath::Pi();
const Float_t m_pi_0 = 135.98; // MeV/c^2
const Float_t m_pi_char = 139.57; // MeV/c^2


// common ROOT vectors
using RVecF = ROOT::VecOps::RVec<Float_t>;
using RVecI = ROOT::VecOps::RVec<Int_t>;
using RVecD = ROOT::VecOps::RVec<Double_t>;
using RVecStr = ROOT::VecOps::RVec<TString>;
using RVecChar = ROOT::VecOps::RVec<const char*>;
using RVecLor = ROOT::VecOps::RVec<TLorentzVector>;
using LorCyl = ROOT::Math::PtEtaPhiMVector;
using RVecLorCyl = ROOT::VecOps::RVec<LorCyl>;
using array_4 = std::array<Float_t, 4>;
using array_5 = std::array<Float_t, 5>;
using RVecF_4 = ROOT::VecOps::RVec<array_4>;
using RVecF_5 = ROOT::VecOps::RVec<array_5>;

//Other definitions
using RDF = ROOT::RDataFrame;
using RN = ROOT::RDF::RNode;


// for canvas saving
struct Drawable {
    TObject* obj;
    TString drawOpt;
};
using RVecDraw = ROOT::VecOps::RVec<Drawable>;
