#ifndef UTILS_H
#define UTILS_H

#include <array>
#include <cmath>
#include <vector>

const std::array<int,4> gmumu = {1,-1,-1,-1};
const std::array<int,4> t_vector = {1,0,0,0};
int g(int mu, int nu);

const double Gevtofm = 5.067728853;
const double hbarC = 1. / 5.067728853; //=0.197 Gev*fm
const double PI = std::acos(-1);

int levi(int i, int j, int k, int l);
int statistics(double spin);
std::vector<double> linspace(double min, double max, int size);
std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> file_to_4_column(std::string filename, int f_column);

double dotProduct(const std::array<double, 3>& a, const std::array<double, 3>& b);

struct NonZeroLevi {
    int mu, nu, rh, sg;
    double value;
};

std::vector<NonZeroLevi> getNonZeroLevi();

class trilinearInterpolator {
private:
    int npT, nphi, ny;
    double pTmin, phimin, ymin, pTmax, phimax, ymax, dpT, dphi, dy;
    std::vector<double> pT_grid, phi_grid, y_grid;

public:
    trilinearInterpolator(const std::vector<double>& pG_, const std::vector<double>& phG_, const std::vector<double>& yG_);
    double trilinearInterpolation(const double pT, const double phi, const double y_rap, const std::vector<std::vector<std::vector<double>>>& P_primary);
    std::vector<double> trilinearInterpolation(const double pT, const double phi, const double y_rap, const std::vector<std::vector<std::vector<std::vector<double>>>>& P_primary);
};

#endif
