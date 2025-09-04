#include <typeinfo>
#include <stdio.h>
#include <iostream>
#include <vector>
#include "TMath.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "ROOT/RVec.hxx"

#include "../../lib/computations.h"
#include "../../lib/cut_branch.h"



// command
// !g++ src/branch_plots/dedx_p_charge.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/branch_plots/dedx_p_charge; ./bin/branch_plots/dedx_p_charge



/*
	calculates the specific energy loss against the rigidity 

	Cut1: cut data based on n (4) tracks
	Cut2: cut based on the the value of the transverse momentum
	
	Plots dE/dx against p 
*/
int main()
{  

	const char* file_name_20 = "data/TOTEM20.root";
	const char* file_name_40 = "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM40.root"; 
  const char* tree_name = "tree";	

	//const char* file_name_20= "data/glueball_mass_reconstruction/data_20.root";
  	//const char* file_name_40= "data/glueball_mass_reconstruction/data_40.root";
    //const char* tree_name = "tree1";
	
	// check the identifier type of the entries to know 
	// how to use RVec
	// TFile *f = TFile::Open(file_name);
	// TTree *t = (TTree*)f->Get(tree_name);
	// t->Print();

    // Load in the Tree

	ROOT::EnableImplicitMT();
	RDF df_20_df(tree_name, file_name_20);
	RDF df_40_df(tree_name, file_name_40);

	// cut on track number 
    Int_t track_num = 4;
    VectorConstraint<Int_t> track_constraint{
        .vector_condition = [track_num](Int_t ntrk) { return ntrk == track_num; },
        .op = Operator::all
    };
    RN df_20 = cut_branch(df_20_df, "ntrk", track_constraint, true);
    RN df_40 = cut_branch(df_40_df, "ntrk", track_constraint, true);

    // cut on transverse momentum
	float_t low_pt_threshold = 0.8;
    VectorConstraint<float_t> momentum_constraint{
        .vector_condition = [low_pt_threshold](float_t pt) { return pt < low_pt_threshold; },
        .op = Operator::all
    };
    df_20 = cut_branch(df_20, "alltrk_pt", momentum_constraint, true);
    df_40 = cut_branch(df_40, "alltrk_pt", momentum_constraint, true);


    // cut on track charge
	Int_t num_pos = 2, num_neg = 2;
    VectorConstraint<Int_t> charge_constraint_pos{
        .vector_condition = [](Int_t q) { return q == +1; },
        .op = Operator::equal,
        .element_count = num_pos
    };
    VectorConstraint<Int_t> charge_constraint_neg{
        .vector_condition = [](Int_t q) { return q == -1; },
        .op = Operator::equal,
        .element_count = num_neg
    };
    df_20 = cut_branch(df_20, "trk_q", charge_constraint_pos, true);
    df_20 = cut_branch(df_20, "trk_q", charge_constraint_neg, true);
    df_40 = cut_branch(df_40, "trk_q", charge_constraint_pos, true);
    df_40 = cut_branch(df_40, "trk_q", charge_constraint_neg, true);
	









	//Want to plot  momemetum against de/dx

    TCanvas *c1 = new TCanvas("c1", "trk_dE/dx vs p");

    // Create the 2D histogram
    TH2F *histo_20 = new TH2F("histo_20", "trk_dedx vs p;trk_p (MeV/c);dE/dx (MeV/m)",
                               1000, 0, 2.5,  // x-axis bins and range
                               1000, 0, 10); // y-axis bins and range
    TH2F *histo_40 = new TH2F("histo_40", "trk_dedx vs p;trk_p (MeV/c);dE/dx (MeV/m)",
                               100, 0,10,  // x-axis bins and range
                               1000, 0, 1e-6); // y-axis bins and range

    // Fill the histogram from RVecs
    df_20.Foreach(
        [histo_20] (const RVecF &p, const RVecF &dedx) 
        {int n = p.size();
		for (int i = 0; i < n; ++i)
    	histo_20->Fill(p[i], dedx[i]);
		},
	{"trk_p", "trk_dedx"}
	);
	
	histo_20->Draw("COLZ");
	c1->Update();	
	TFile* f20 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/branch_plots/dedx_p_20_charge.root", "RECREATE");
    c1->Write();
    f20->Close();	
	c1->Clear();

    df_40.Foreach(
        [histo_40] (const RVecF &p, const RVecF &dedx) 
        {int n = p.size();
		for (int i = 0; i < n; ++i)
    	histo_40->Fill(p[i], dedx[i]);
		},
	{"trk_p", "trk_dedx"}
);
    TFile* f40 = new TFile("media/root_files/Younes_NTuples_TOTEM20_and_40/branch_plots/dedx_p_40_charge.root", "RECREATE");
    histo_40->Draw("COLZ");
	c1->Update();	
    histo_40->Write();
 

//plot_dedx->Draw("COLZ");
//c1->Update();
//c1->SaveAs("media/photos/momentum_against_dedx.pdf");

	return 0;
}
