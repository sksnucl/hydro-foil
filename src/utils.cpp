#include<tuple>
#include<iostream>
#include <algorithm> 
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <cmath>
#include "utils.h"

int levi(int i, int j, int k, int l){
// Levi-Civita symbols
// i,j,k,l = 0...3 i.e. upper indices
 if((i==j)||(i==k)||(i==l)||(j==k)||(j==l)||(k==l)) return 0;
 else return ( (i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12 );
}

int statistics(double spin) {
    // returns 1 if Fermi, -1 if Bose, 0 if invalid
    double twice_spin = 2.0 * spin;
    double nearest = std::round(twice_spin);
    if (std::abs(twice_spin - nearest) > 1e-6) return 0; // not valid spin

    int n = static_cast<int>(nearest); // 2*spin as int
    return (n % 2 == 1) ? 1 : -1; // odd → Fermion (1), even → Boson (-1)
}

//double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
//    if (a.size() == 3 && b.size() == 3) {
//        // Standard Euclidean dot product for 3D vectors
//        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//    } 
//    else if (a.size() == 4 && b.size() == 4) {
//        // Minkowski metric dot product for 4-vectors: (+,-,-,-)
//        return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
//    } 
//    else {
//        std::cerr << "Error: dotProduct requires two vectors of size 3 or 4." << std::endl;
//        exit(1);
//    }
//}

double dotProduct(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    // Standard Euclidean dot product for 3D vectors
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

int g(int mu, int nu){
    if(mu > 3 || mu < 0 || nu > 3 || nu < 0){
        std::cout<<"error with the indices of the metric"<<std::endl;
        exit(1);
    }
    if(mu==nu){
        if(mu==0){
            return 1;
        }
        else{
            return -1;
        }
    }
    return 0;
}


std::vector<NonZeroLevi> getNonZeroLevi() {
    std::vector<NonZeroLevi> nonzerolevi;
    
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rh = 0; rh < 4; rh++) {
                for (int sg = 0; sg < 4; sg++) {
                    double value = levi(mu, nu, rh, sg); 
                    if (value != 0) {
                        NonZeroLevi element;
                        element.mu = mu;
                        element.nu = nu;
                        element.rh = rh;
                        element.sg = sg;
                        element.value = value;
                        nonzerolevi.push_back(element);
                    }
                }
            }
        }
    }

    return nonzerolevi;
}

std::vector<double> linspace(double min, double max, int size){
    //returns a vector of doubles from min to max (included), equally spaced
    std::vector<double> result;
    double interval = (max-min)/(size-1);
    double tmp_min = min;
    for( int i=0;i<size;i++){
        result.push_back(tmp_min);
        tmp_min+=interval;
    }
    return result;
}

trilinearInterpolator::trilinearInterpolator(const std::vector<double>& pG_, const std::vector<double>& phG_, const std::vector<double>& yG_)
    : pT_grid(pG_), phi_grid(phG_), y_grid(yG_) {
    if (!std::is_sorted(pT_grid.begin(), pT_grid.end())) {
        std::sort(pT_grid.begin(), pT_grid.end());
    }
    if (!std::is_sorted(phi_grid.begin(), phi_grid.end())) {
        std::sort(phi_grid.begin(), phi_grid.end());
    }
    if (!std::is_sorted(y_grid.begin(), y_grid.end())) {
        std::sort(y_grid.begin(), y_grid.end());
    }
    npT = pT_grid.size();
    nphi = phi_grid.size();
    ny = y_grid.size();
    pTmin = pT_grid.front();
    pTmax = pT_grid.back();
    phimin = phi_grid.front();
    phimax = phi_grid.back();
    ymin = y_grid.front();
    ymax = y_grid.back();
    dpT = (pTmax - pTmin) / (npT - 1);
    dphi = (phimax - phimin) / (nphi - 1);
    dy = (ymax - ymin) / (ny - 1);
}

double trilinearInterpolator::trilinearInterpolation(const double pT, const double phi, const double y_rap, const std::vector<std::vector<std::vector<double>>>& P_primary) {
    size_t i_pT = (pT - pTmin) / dpT;
    size_t i_phi = (phi - phimin) / dphi;
    size_t i_y = (y_rap - ymin) / dy;

    i_pT = std::max(size_t(0), std::min(i_pT, pT_grid.size() - 2));
    i_phi = std::max(size_t(0), std::min(i_phi, phi_grid.size() - 2));
    i_y = std::max(size_t(0), std::min(i_y, y_grid.size() - 2));

    double pT_frac = (pT - pT_grid[i_pT]) / dpT;
    double phi_frac = (phi - phi_grid[i_phi]) / dphi;
    double y_frac = (y_rap - y_grid[i_y]) / dy;
    
    pT_frac = std::max(0.0, std::min(pT_frac, 1.0));
    phi_frac = std::max(0.0, std::min(phi_frac, 1.0));
    y_frac = std::max(0.0, std::min(y_frac, 1.0));

    double c00 = P_primary[i_pT][i_phi][i_y] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi][i_y] * pT_frac;
    double c01 = P_primary[i_pT][i_phi + 1][i_y] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi + 1][i_y] * pT_frac;
    double c10 = P_primary[i_pT][i_phi][i_y + 1] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi][i_y + 1] * pT_frac;
    double c11 = P_primary[i_pT][i_phi + 1][i_y + 1] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi + 1][i_y + 1] * pT_frac;

    double c0 = c00 * (1 - phi_frac) + c01 * phi_frac;
    double c1 = c10 * (1 - phi_frac) + c11 * phi_frac;

    return c0 * (1 - y_frac) + c1 * y_frac;
}

std::vector<double> trilinearInterpolator::trilinearInterpolation(
    const double pT, const double phi, const double y_rap, 
    const std::vector<std::vector<std::vector<std::vector<double>>>>& P_primary
) {
    size_t i_pT = (pT - pTmin) / dpT;
    size_t i_phi = (phi - phimin) / dphi;
    size_t i_y = (y_rap - ymin) / dy;

    i_pT = std::max(size_t(0), std::min(i_pT, P_primary.size() - 2));
    i_phi = std::max(size_t(0), std::min(i_phi, P_primary[0].size() - 2));
    i_y = std::max(size_t(0), std::min(i_y, P_primary[0][0].size() - 2));

    double pT_frac = (pT - pT_grid[i_pT]) / dpT;
    double phi_frac = (phi - phi_grid[i_phi]) / dphi;
    double y_frac = (y_rap - y_grid[i_y]) / dy;

    pT_frac = std::max(0.0, std::min(pT_frac, 1.0));
    phi_frac = std::max(0.0, std::min(phi_frac, 1.0));
    y_frac = std::max(0.0, std::min(y_frac, 1.0));

    // Determine number of components dynamically
    size_t num_components = P_primary[0][0][0].size();
    std::vector<double> result(num_components, 0.0);

    for (size_t comp = 0; comp < num_components; ++comp) {
        double c00 = P_primary[i_pT][i_phi][i_y][comp] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi][i_y][comp] * pT_frac;
        double c01 = P_primary[i_pT][i_phi + 1][i_y][comp] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi + 1][i_y][comp] * pT_frac;
        double c10 = P_primary[i_pT][i_phi][i_y + 1][comp] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi][i_y + 1][comp] * pT_frac;
        double c11 = P_primary[i_pT][i_phi + 1][i_y + 1][comp] * (1 - pT_frac) + P_primary[i_pT + 1][i_phi + 1][i_y + 1][comp] * pT_frac;

        double c0 = c00 * (1 - phi_frac) + c01 * phi_frac;
        double c1 = c10 * (1 - phi_frac) + c11 * phi_frac;

        result[comp] = c0 * (1 - y_frac) + c1 * y_frac;
    }

    return result;
}
