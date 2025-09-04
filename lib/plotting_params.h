// plotting_style.h
#pragma once

#include "custom_definitions.h"
#include "TObject.h"
#include "Rtypes.h"

inline const RVecI color_scheme = {kP10Cyan, kP10Ash, kP10Green, kP10Orange, kP10Brown, kP10Violet, kP10Gray, kP10Red, kP10Yellow, kP10Blue};

void SetPlotStyle();


// Forward declarations
class TH1;
class TH2;
class THStack;
class TLine;

enum LineOrientation { kVertical, kHorizontal, kDiagonal };
TLine* DrawLine(Float_t x1, Float_t y1, Float_t x2, Float_t y2,
                Color_t color = kBlack, Int_t width = 2, Int_t style = kDashed);

TLine* DrawLine(Float_t val1, Float_t val2, TObject* obj,
                LineOrientation orientation = kVertical,
                Color_t color = kRed, Int_t width = 2, Int_t style = kDashed,
                Float_t offsetFactor = 1.00);

