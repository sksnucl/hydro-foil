#ifndef INTEGRALS_H
#define INTEGRALS_H

#include <vector>
#include <map>
#include "surface.h"
#include "utils.h"
#include "particle.h"

double momentumCM(double s, double m1, double m2);

std::vector<double> ResonanceThreeMomentum(const double MR, const double E, const std::vector<double>& p, const double Es, const std::vector<double>& ps);

double Jacobian(const double Mm, const double E, const std::vector<double>& p, const double Mp, const double Es, const std::vector<double>& ps);

std::vector<double> RestFrameSpinVector(const double mass, const std::vector<double>& p, const std::vector<double>& S) ;

void sum_over_surface(const std::vector<element> &freeze_out_sup, const Particle& particle, const std::array<double, 4> p, double& dndp, double P_vorticity[4], double P_shear[4]);

double aux_exact_polarization(double spin, double pu, double T, double mutot, double abs_theta);

void sum_over_surface_exact(const std::vector<element> &freeze_out_sup, const Particle& particle, const std::array<double, 4> p, double& dndp, double P_vorticity[4], double P_shear[4]);

void compute_primary(const std::vector<double>& pT_vec, const std::vector<double>& phi_vec, const std::vector<double>& y_vec, std::map<int, Particle>& particles, std::vector<element>& freeze_out_sup, const bool exact, const bool decay, int pdg_id);

void compute_polarization_feeddown(const std::vector<double>& pT_vec, const std::vector<double>& phi_vec, const std::vector<double>& y_vec, Particle& primary, std::map<int, Particle>& particles);

#endif
