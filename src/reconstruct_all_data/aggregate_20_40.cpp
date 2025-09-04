#include "ROOT/RDataFrame.hxx"

// include headers containing useful functions and definitions
#include "../../lib/computations.h"
#include "../../lib/cut_branch.h"
#include "../../lib/common_cuts.h"

// command
// !g++ src/reconstruct_all_data/aggregate_20_40.cpp lib/*.cpp -Llib -ldict `root-config --cflags --libs` -o bin/reconstruct_all_data/aggregate_20_40; ./bin/reconstruct_all_data/aggregate_20_40

int main() {

    // Load in all 20 and 40 combined type events, get column names, enable multithreading
    ROOT::EnableImplicitMT();

    std::vector<std::string> files_20 = {
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM20.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM21.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM22.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM23.root"
    };

    std::vector<std::string> files_40 = {
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM40.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM41.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM42.root",
        "/eos/cms/store/group/phys_diffraction/CMSTotemLowPU2018/YounesNtuples/TOTEM43.root"
    };
    RDF df_20("tree", files_20);
    RDF df_40("tree", files_40);


	// because usefull define the SNR (signicance) branches here
	
	// while it could be useful the extra 500mb just arent worht the extra wait
	df_20.Define("trk_dz_snr", [] (Int_t ntrk, RVecF z, RVecF err)
			{
				RVecF result(ntrk);
        		for (Int_t i = 0; i < ntrk; ++i) {
           			 result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
				}
				return result;
			}, {"ntrk", "trk_dz", "trk_dzerr"})
			
		.Define("trk_dxy_snr", [] (Int_t ntrk, RVecF z, RVecF err)
			{
				RVecF result(ntrk);
        		for (Int_t i = 0; i < ntrk; ++i) {
           			 result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
				}
				return result;
			}, {"ntrk", "trk_dxy", "trk_dxyerr"})
		.Snapshot("tree", "data/glueball_mass_reconstruction/data_20_uncut.root");
	


	df_40.Define("trk_dz_snr", [] (Int_t ntrk, RVecF z, RVecF err)
			{
				RVecF result(ntrk);
        		for (Int_t i = 0; i < ntrk; ++i) {
           			 result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
				}
				return result;
			}, {"ntrk", "trk_dz", "trk_dzerr"})
	.Define("trk_dxy_snr", [] (Int_t ntrk, RVecF z, RVecF err)
			{
				RVecF result(ntrk);
        		for (Int_t i = 0; i < ntrk; ++i) {
           			 result[i] = (err[i] != 0.0f) ? z[i] / err[i] : 0.0f;
				}
				return result;
			}, {"ntrk", "trk_dxy", "trk_dxyerr"})
	.Snapshot("tree", "data/glueball_mass_reconstruction/data_40_uncut.root");

}
