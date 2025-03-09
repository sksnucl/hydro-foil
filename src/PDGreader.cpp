#include "PDGreader.h"

void read_pdg_file(const std::string &filename, std::map<int, Particle> &particles, int pt_bins, int phi_bins, int y_bins) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        Particle p(pt_bins, phi_bins, y_bins);
        
        iss >> p.pdg_id >> p.name >> p.mass >> p.width >> p.gspin >> p.baryon 
            >> p.strange >> p.charm >> p.bottom >> p.gisospin >> p.charge >> p.Nd;
        
        for (int i = 0; i < p.Nd; ++i) {
            DecayChannel dc;
            std::getline(file, line);
            std::istringstream iss_decay(line);
            int res_id;
            iss_decay >> res_id >> dc.num_daughters >> dc.branching_ratio;

            for (int j = 0; j < dc.num_daughters; ++j) {
                int daughter_id;
                iss_decay >> daughter_id;
                dc.daughters.push_back(daughter_id);
            }
                
            p.decay_channels.push_back(dc);
                
            for (int j = dc.num_daughters; j < 5; ++j) {
                int dummy;
                iss_decay >> dummy;
            }
        }

        particles[p.pdg_id] = p;
    }
    
    //Add anti-baryons too. Not needed for this work

    file.close();
}
