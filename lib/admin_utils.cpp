#include "admin_utils.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TObject.h"

void SaveCanvas(TCanvas* c, TObject* obj, const char* drawOpt, const char* fileName, const char* SaveType) {
    c->cd();	
    obj->Draw(drawOpt);
    c->Update();
    TFile f(fileName, SaveType);
	obj->Write();
    c->Write();
    c->Clear();
}

void SaveOnlyCanvas(TCanvas* c, const char* fileName, const char* SaveType) {
    c->cd();	
	c->Update();
    TFile f(fileName, SaveType);
    c->Write();
    c->Clear();
}

void SaveCanvas(TCanvas* c,
	const RVecDraw& drawObjects,
	const char* fileName,
	const char* saveType) {
    
	c->cd();
    Bool_t first = true;
    
	for (const auto& d : drawObjects) {
        if (!d.obj) continue;

        TString opt = d.drawOpt;
        if (!first) opt += " SAME";
        d.obj->Draw(opt);
        first = false; }
    
	c->Update();

	gPad->Update(); 
	gPad->Modified(); 
    TFile f(fileName, saveType);
    c->Write();
    
	for (const auto& d : drawObjects) {
        gPad->Update(); 
		gPad->Modified();
		d.obj->Write(); }
    f.Close();
    c->Clear();
};
void SaveRooFitObjects(TCanvas* c,
						const ROOT::VecOps::RVec<TObject*>& objects,
                        const TString& fileName,
                        RooPlot* frame,
                        const TString& canvasName)
	{

    if (!c) {
        std::cerr << "Canvas pointer is null!" << std::endl;
        return;
    }

    c->cd();
    c->Update();

    TFile f(fileName, "RECREATE");
    if (f.IsZombie()) {
        std::cerr << "Cannot open file: " << fileName << std::endl;
        return;
    }

    c->Write(canvasName, TObject::kOverwrite);

    for (auto obj : objects) {
        if (!obj) continue;

        // Ensure the object has a name
        if (!obj->GetName() || strlen(obj->GetName()) == 0) {
            TNamed* named = dynamic_cast<TNamed*>(obj);
            if (named) {
                TString autoName = TString(named->ClassName()) + "_obj";
                named->SetName(autoName);
            }
        }
        obj->Write(obj->GetName(), TObject::kOverwrite);
    }

    if (frame) frame->Write(frame->GetName(), TObject::kOverwrite);

    f.Close();

    //c->Clear();
}

