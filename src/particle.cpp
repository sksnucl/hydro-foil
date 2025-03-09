#include "particle.h"

Particle::Particle() {
    mass = 0.0;
    width = 0.0;
    gspin = 1;
    baryon = 0;
    strange = 0;
    charm = 0;
    bottom = 0;
    gisospin = 1;
    charge = 0;
    Nd = 0;

    EdN_d3p_primary.resize(0, std::vector<std::vector<double>>(0, std::vector<double>(0, 0.0)));
    EdN_d3p_feeddown.resize(0, std::vector<std::vector<double>>(0, std::vector<double>(0, 0.0)));
    Pv_primary.resize(0, std::vector<std::vector<std::vector<double>>>(0, std::vector<std::vector<double>>(0, std::vector<double>(0, 0.0))));
    Ps_primary.resize(0, std::vector<std::vector<std::vector<double>>>(0, std::vector<std::vector<double>>(0, std::vector<double>(0, 0.0))));
    Pv_feeddown.resize(0, std::vector<std::vector<std::vector<double>>>(0, std::vector<std::vector<double>>(0, std::vector<double>(0, 0.0))));
    Ps_feeddown.resize(0, std::vector<std::vector<std::vector<double>>>(0, std::vector<std::vector<double>>(0, std::vector<double>(0, 0.0))));
}

Particle::Particle(int pt_bins, int phi_bins, int y_bins) {
    EdN_d3p_primary.resize(pt_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>(y_bins, 0.0)));
    EdN_d3p_feeddown.resize(pt_bins, std::vector<std::vector<double>>(phi_bins, std::vector<double>(y_bins, 0.0)));
    Pv_primary.resize(pt_bins, std::vector<std::vector<std::vector<double>>>(phi_bins, std::vector<std::vector<double>>(y_bins, std::vector<double>(4, 0.0))));
    Ps_primary.resize(pt_bins, std::vector<std::vector<std::vector<double>>>(phi_bins, std::vector<std::vector<double>>(y_bins, std::vector<double>(4, 0.0))));
    Pv_feeddown.resize(pt_bins, std::vector<std::vector<std::vector<double>>>(phi_bins, std::vector<std::vector<double>>(y_bins, std::vector<double>(3, 0.0))));
    Ps_feeddown.resize(pt_bins, std::vector<std::vector<std::vector<double>>>(phi_bins, std::vector<std::vector<double>>(y_bins, std::vector<double>(3, 0.0))));
}
