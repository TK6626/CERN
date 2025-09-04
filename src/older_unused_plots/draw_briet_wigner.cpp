#include "../lib/computations.h"
#include "TCanvas.h"
#include "TF1.h"




int main() {
// Wrap in a TF1
    
	TCanvas* canvas = new TCanvas("canvas");
	TF1 *bw = new TF1("bw", fit_breit_wigner, 0, 2000, 4); // Energy range [0, 2000] MeV

    // Set parameter names and values
    bw->SetParNames("Mass", "Width", "Amplitude", "Background");
    bw->SetParameters(1000, 100, 1e12, 0); // Example values

    // Styling (optional)
    bw->SetLineColor(kBlue + 1);
    bw->SetLineWidth(2);
	bw->SetNpx(200);

    // Draw the function
    bw->SetTitle("Breit-Wigner Resonance;Energy [MeV];f(E)");
    bw->Draw();

	TFile* outFile = new TFile("media/root_files/Younes_NTuples_TOTEM20/glueball_mass_reconstruction/draw_breit_wigner.root", "RECREATE");
    canvas->Update();
    canvas->Write();
    outFile->Close();
    canvas->Clear();

    return 0;
}

