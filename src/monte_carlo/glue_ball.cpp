#include <stdio.h>
#include <array>

#include "TRandom3.h"

#include "../../lib/computations.h"
#include "../../lib/common_cuts.h"


/**
 * Simple Monte Carlo Simluation of the vector (J=2)  Glueball decay 
 * using most common decay channels
 *
 *	F -> ρ ρ	-> π+ π− π+ π−				(1)
 *	F -> ϕ ϕ 	-> K+ K− K+ K−				(0.84)
 *	F -> K∗ K∗ 	-> K+ π− K− π+				(0.32)
 *  F -> ω ω 	-> π+ π− π0 ->π+ π− γγ 		(0.11)
 *
 */

void glue_ball()
{
	// Initiate random number generator, and number of decays
	TRandom3 rng(0);
	const Int_t N = 100000000;
	const Int_t bn = 4; // number of branches we are considering
	
	// Set up Decay Counters
	std::array<Int_t, bn> Counter = {0,0,0,0};

	// Set up relative decay rates
	std::array<Float_t, bn> relative_branch_ratios  = {1, 0.84, 0.32, 0.11};
	
	// calculate total relative probability	
	Double_t total = 0.0;
    for (int i = 0; i < bn; ++i) {
        total += relative_branch_ratios[i];
    }

	// normalise to sum to 1 (or here include other posibilities of decays
	std::array<Float_t, bn> branch;
	
	for (int i = 0; i < bn; ++i) {
		branch[i] = relative_branch_ratios[i] / total; 
	}
	
	// create cumulative distribution function 
	std::array<Float_t, bn> cumulative;
    cumulative[0] = branch[0];
    for (int i = 1; i < bn; ++i) {
        cumulative[i] = cumulative[i - 1] + branch[i];
    }


	for (Int_t i = 0; i < N; ++i)
	{
        Double_t r = rng.Uniform();

        for (int j = 0; j < bn; ++j)
		{
            if (r < cumulative[j])
			{
                Counter[j]++;
                break;
            }
        }
	}


	std::cout << "Kaon decay MC (" << N << " events):\n";
    for (int i = 0; i < bn; ++i) {
        std::cout << "  Branch " << i << ": " << Counter[i] << " decays ("
                  << 100.0 * Counter[i] / N << "%)\n";
    }


}

