#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom3.h"

const double me = 0.511;       // Electron mass [MeV]
const double E0 = 1.0;         // Endpoint energy [MeV]
const int N_EVENTS = 100000;   // Number of simulated decays

// Beta spectrum: N(E) ∝ p E (E0 - E)^2
double beta_spectrum(double E) {
    if (E <= 0 || E >= E0) return 0;
    double p = std::sqrt(E * E + 2 * E * me);
    return p * E * std::pow(E0 - E, 2);
}

double sample_energy(TRandom3 &rng, double ymax) {
    double E, y;
    do {
        E = rng.Uniform(0, E0);
        y = rng.Uniform(0, ymax);
    } while (y > beta_spectrum(E));
    return E;
}

void example_beta_decay() {
    TRandom3 rng(0);  // ROOT random number generator

    // Estimate ymax for rejection sampling (not exact)
    double ymax = beta_spectrum(E0 / 2);

    TH1D *h = new TH1D("h", "Beta Decay Electron Spectrum;Electron Energy (MeV);Counts", 100, 0, E0);

    for (int i = 0; i < N_EVENTS; ++i) {
        double E = sample_energy(rng, ymax);
        h->Fill(E);
    }

    TCanvas *c = new TCanvas("c", "Beta Decay Spectrum", 800, 600);
    h->Draw();
    c->SaveAs("beta_decay_spectrum.png");
}

