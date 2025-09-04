#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"

#include "../../lib/custom_definitions.h"
#include "../../lib/cut_branch.h"
#include "../../lib/apply_cuts.cpp"
#include "../../lib/plotting_params.h"
#include "../../lib/admin_utils.h"

// !./bashing/run_file.sh phi_phi_reconstruction/generic_data


/** make cuts on data to 
 * try and remove extraneous kaon backroudn
 * due arising from secodary interactions
 */

void AddTrkSNR(RN &df);

int main() {

	// relevent constnats
	const Float_t kaon_mass = m_kaon_char * 1e-3;

    ROOT::EnableImplicitMT();
    std::vector<std::string> SNR_uncut_20 = {
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM21.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM22.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM23.root"
    };
    std::vector<std::string> SNR_uncut_40 = {
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM40.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM41.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM42.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM43.root"
    };
   
	 ROOT::VecOps::RVec<std::vector<std::string>> file_bunch = {SNR_uncut_20, SNR_uncut_40};
    RVecStr file_name = {"uncut_SNR_20.root", "uncut_SNR_40.root"};

	 Int_t i = 0;
	for (const std::vector<std::string> files : file_bunch) {
	


	RDF df_df("tree", files);

	std::cout << df_df.Describe() << std::endl;
	// for the  ΦΦ -> K+K-K+K- we require 4 tracks, and precisely 2+ve charges, 2-ve charges  -apply that here
	RN df = CutApplier(df_df)
		.apply_ntrk_cut(4) // there are precisely 4 tracks	
		.apply_charge_cut() // precisely 2+ve 2-ve
		.result();


	// because usefull define the SNR (signicance) branches here
	df.Define("trk_dz_snr", [] (Int_t ntrk, RVecF z, RVecF err)
			{
				RVecF result(ntrk);
        		for (Int_t i = 0; i < ntrk; ++i) {
           			 result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
				}
				return result;
			}, {"ntrk", "trk_dz", "trk_dzerr"})
			
	 	// Define trk_dxy snr for use later		
		.Define("trk_dxy_snr", [] (Int_t ntrk, RVecF z, RVecF err)
			{
				RVecF result(ntrk);
        		for (Int_t i = 0; i < ntrk; ++i) {
           			 result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
				}
				return result;
			}, {"ntrk", "trk_dxy", "trk_dxyerr"})
	

		// Define four momentum just assume that the tracks at the moment are indeed kaons
		.Define("kaon_four_momentum",
        [kaon_mass](const RVecF &trk_pt, const RVecF &trk_eta, const RVecF &trk_phi, RVecI& isK) {
			RVecLorCyl kaon_4p {4};
            for (Int_t i = 0; i < 4; ++i) {
			//	if (isK[i]){
                	LorCyl p_kaon{trk_pt[i], trk_eta[i], trk_phi[i], kaon_mass};
                	kaon_4p[i] = p_kaon;
			//	}
            }
				return kaon_4p;

        	},
        	{"trk_pt", "trk_eta", "trk_phi", "trk_isK"}
    		)
		.Define("phi_four_momentum",[](const RVecI &trk_q, const RVecLorCyl &kaon_4p) {
		
			std::array<Int_t, 2> kaon_neg, kaon_pos;
            Int_t pos_count = 0, neg_count = 0;
           
			// check for charge matching (this is redundant as it has been vetoed on but make sure
			for (Int_t i = 0; i < 4; ++i) {
                if (trk_q[i] == +1 && pos_count < 2) {
                    kaon_pos[pos_count++] = i;
                } 
				else if (trk_q[i] == -1 && neg_count < 2) {
                    kaon_neg[neg_count++] = i;
                }
            }

            if (!(pos_count == 2 && neg_count == 2)) {
				std::cout << " something has gone wrong" << std::endl;
                return RVecLorCyl {};
            }	

			// create possible parings for di phi- production 
            LorCyl P_00 = kaon_4p[kaon_pos[0]] + kaon_4p[kaon_neg[0]];
            LorCyl P_11 = kaon_4p[kaon_pos[1]] + kaon_4p[kaon_neg[1]];

			LorCyl P_01 = kaon_4p[kaon_pos[0]] + kaon_4p[kaon_neg[1]];
            LorCyl P_10 = kaon_4p[kaon_pos[1]] + kaon_4p[kaon_neg[0]];
            
			return RVecLorCyl{P_00, P_11, P_10, P_01}; // each pari (00,11) (01,10) are pairings that are consitent with charge
			
			},{"trk_q", "kaon_four_momentum"}
		)
		.Snapshot("tree", TString::Format("data/phi_phi_reconstruction/" + file_name[i]));
	++i;	
	}
	
	return 0;
}
