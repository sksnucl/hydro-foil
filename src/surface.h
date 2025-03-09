#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
#include <iostream>
#include <array>

//meant to read a the beta.dat file from the output of vhlle (polarization branches)

struct element {
    double tau=0, x=0, y=0, eta=0, ed=0, p=0;
    std::array<double,4> u={1,0,0,0};
    std::array<double,4> dsigma={0};
    double T=0, mub=0, muq=0, mus=0; //GeV
    std::array<std::array<double,4>,4> dbeta={0};
};

void read_hypersurface(std::string filename, std::vector<element> &hypersurface);

#endif
