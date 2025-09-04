#pragma once

#include "custom_definitions.h"
#include <functional>

class TCanvas;
class TObject;

void SaveCanvas(TCanvas* c,const RVecDraw& drawObjects,const char* fileName,const char* saveType);
void SaveCanvas(TCanvas* c, TObject* obj, const char* drawOpt, const char* fileName, const char* SaveType);
void SaveOnlyCanvas(TCanvas* c, const char* fileName, const char* SaveType);
void SaveRooFitObjects(TCanvas* c,
						const ROOT::VecOps::RVec<TObject*>& objects,
                        const TString& fileName,
                        RooPlot* frame = nullptr,
                        const TString& canvasName = "canvas");
