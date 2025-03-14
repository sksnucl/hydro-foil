#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

struct DecayChannel {
    int num_daughters;
    double branching_ratio;
    std::vector<int> daughters;
};

class Particle {
private:
    int pt_bins, phi_bins, y_bins;
public:
    int pdg_id, gspin, baryon, strange, charm, bottom, gisospin, charge, Nd;
    double mass, width;
    std::string name;
    
    std::vector<std::vector<std::vector<double>>> EdN_d3p_primary;
    std::vector<std::vector<std::vector<std::vector<double>>>> Pv_primary, Ps_primary;
    std::unordered_map<int, std::vector<std::vector<std::vector<double>>>> EdN_d3p_feeddown;
    std::unordered_map<int, std::vector<std::vector<std::vector<std::vector<double>>>>> Pv_feeddown, Ps_feeddown;
    
    std::vector<DecayChannel> decay_channels;
    
    Particle();

    Particle(int pt_bins, int phi_bins, int y_bins);
    void allocate_feeddown(int mother_pdg);
    
    double get_mass() const {return mass;}
    double get_width() const {return width;}
    int get_gspin() const {return gspin;}
    int get_baryon() const {return baryon;}
    int get_strange() const {return strange;}
    int get_charm() const {return charm;}
    int get_bottom() const {return bottom;}
    int get_gisospin() const {return gisospin;}
    int get_charge() const {return charge;}
    int get_Nd() const {return Nd;}
    
    const std::vector<std::vector<std::vector<double>>>& get_dndp_primary() const { return EdN_d3p_primary; }
    const std::vector<std::vector<std::vector<std::vector<double>>>>& get_Pv_primary() const { return Pv_primary; }
    const std::vector<std::vector<std::vector<std::vector<double>>>>& get_Ps_primary() const { return Ps_primary; }
};

#endif // PARTICLE_H
