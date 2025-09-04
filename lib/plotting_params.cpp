
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TList.h"
#include "TIterator.h"
#include "THStack.h"

#include "plotting_params.h"

void SetPlotStyle() {
	//gROOT->SetStyle("Plain");
	//gStyle->SetTitleSize(0.05, "T");
    //gStyle->SetTextFont(4);
    //gStyle->SetTitleFont(4, "");      // Global title
    //gStyle->SetLabelFont(4, "XYZ");   // Axis labels
    //gStyle->SetTitleFont(4, "XYZ");   // Axis titles
    //gStyle->SetStatFont(4);
    //gStyle->SetLegendFont(4);

	gStyle->SetLegendFillColor(kWhite);    // Set background color to white
	gStyle->SetLegendFillStyle(1001);      // Solid fill
	gStyle->SetLegendBorderSize(1);        // Optional: 0 = no border, 1 = thin border
    //gStyle->SetOptStat(0);  // Optional: turn off statistics box
    //gStyle->SetPalette(1);  // Optional: set default color palette
}


// Simple line between two points
TLine* DrawLine(Float_t x1, Float_t y1, Float_t x2, Float_t y2,
                Color_t color, Int_t width, Int_t style)
{
    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(color);
    line->SetLineWidth(width);
    line->SetLineStyle(style);
    return line;
}

// Line relative to ROOT object
TLine* DrawLine(Float_t val1, Float_t val2, TObject* obj,
                LineOrientation orientation,
                Color_t color, Int_t width, Int_t style,
                Float_t offsetFactor)
{
    if (!obj) return nullptr;

    Double_t xmin = 0, xmax = 1, ymin = 0, ymax = 1;

    // Determine axis limits
    if (auto h2 = dynamic_cast<TH2*>(obj)) {
        xmin = h2->GetXaxis()->GetXmin();
        xmax = h2->GetXaxis()->GetXmax();
        ymin = h2->GetYaxis()->GetXmin();
        ymax = h2->GetYaxis()->GetXmax();
    } 
    else if (auto h1 = dynamic_cast<TH1*>(obj)) {
        xmin = h1->GetXaxis()->GetXmin();
        xmax = h1->GetXaxis()->GetXmax();
        ymin = h1->GetMinimum();
        ymax = h1->GetMaximum();
    } 
    else if (auto stack = dynamic_cast<THStack*>(obj)) {
        xmin = 1e20; xmax = -1e20;
        ymin = 1e20; ymax = -1e20;

        TList* hlist = stack->GetHists();
        TIter next(hlist);
        TH1* h;

        while ((h = (TH1*)next())) {
            if (h->GetEntries() == 0) continue;
            if (h->GetXaxis()->GetXmin() < xmin) xmin = h->GetXaxis()->GetXmin();
            if (h->GetXaxis()->GetXmax() > xmax) xmax = h->GetXaxis()->GetXmax();
            if (h->GetMinimum() < ymin) ymin = h->GetMinimum();
            if (h->GetMaximum() > ymax) ymax = h->GetMaximum();
        }
    } 

    // Set coordinates based on orientation
    Double_t x1, x2, y1, y2;
    switch (orientation) {
        case kVertical:
            x1 = x2 = val1;
            y1 = ymin;
            y2 = ymax * offsetFactor;
            break;
        case kHorizontal:
            y1 = y2 = val1;
            x1 = xmin;
            x2 = xmax * offsetFactor;
            break;
        case kDiagonal:
            x1 = xmin; y1 = val1;
            x2 = xmax; y2 = val2;
            break;
    }

    TLine* line = new TLine(x1, y1, x2, y2);
    line->SetLineColor(color);
    line->SetLineWidth(width);
    line->SetLineStyle(style);
    return line;
}

