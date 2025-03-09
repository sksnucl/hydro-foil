#ifndef PDG_READER_H
#define PDG_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "particle.h"

void read_pdg_file(const std::string &filename, std::map<int, Particle> &particles, int pt_bins, int phi_bins, int y_bins);

#endif // PDG_READER_H
