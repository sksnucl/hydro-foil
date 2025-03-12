#include "integrals.h"
#include "particle.h"
#include "PDGreader.h"
#include "utils.h"
#include "surface.h"
#include <filesystem>

#ifdef OPEN_MP
    #include<omp.h>
#endif

//TO TURN PARALLEL INTEGRATION OFF(ON) COMMENT(UN-COMMENT) THE OPEN_MP_FLAG LINE OF THE MAKEFILE

int main(int argc, char** argv){

#ifdef OPEN_MP
    int NTHREADS = omp_get_max_threads();
    omp_set_num_threads(NTHREADS);
#endif

    bool decay = false;
    bool exact = false;
    bool midrapidity = false;
    std::string surface_file;
    std::string name_file_primary;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-f" && i + 1 < argc) {
            name_file_primary = argv[i + 1];
            i++;
        } else if (arg == "-s" && i + 1 < argc) {
            surface_file = argv[i + 1];
            i++;
        } else if (arg == "-d") {
            decay = true;
        } else if (arg == "-e") {
            exact = true;
        } else if (arg == "-m") {
            midrapidity = true;
        }
    }
    
    if (surface_file.empty()) {
        std::cout << "Error: You must provide the hypersurface file using -s flag." << std::endl;
        exit(1);
    }
    
    if (name_file_primary.empty()) {
        if (decay){
            name_file_primary = "./primary_feeddown";
        }else{
            name_file_primary = "./primary";
        }
    }
    std::cout << "Data will be saved in the file: " << name_file_primary << std::endl;
   
    int size_pt = 31;
    int size_phi = 21;
    int size_y;
    const double pT_min = 0.0, pT_max = 6.0;
    const double phi_min = 0.0, phi_max = 2*PI;
    const double y_min = -5.0, y_max = 5.0;
    
    std::vector<double> pT = linspace(pT_min,pT_max,size_pt);
    std::vector<double> phi =  linspace(phi_min,phi_max,size_phi);
    std::vector<double> y_rap;
    
    if (midrapidity) {
        size_y = 1;
        y_rap = {0.0};  // Midrapidity, a single value
    } else {
        size_y = 21;
        y_rap = linspace(y_min, y_max, size_y);
    }
    
    std::map<int, Particle> particles;
    //read_pdg_file("pdg-urqmd_v3.3+.dat", particles, size_pt, size_phi, size_y);
    read_pdg_file("./pdg_database/pdg-this_project.dat", particles, size_pt, size_phi, size_y);

    std::vector<element> hypersup = {};
    read_hypersurface(surface_file, hypersup);
    
    compute_primary(pT, phi, y_rap, particles, hypersup, exact, decay, 3122);
    
    if(decay){
        compute_polarization_feeddown(pT, phi, y_rap, particles[3122], particles);
    }
    
    std::ofstream fout(name_file_primary);
    if (!fout.is_open()) {
        std::cerr << "Error: Could not open file " << name_file_primary << std::endl;
        return 1;
    }
    
    for (size_t i = 0; i < size_pt; ++i) {
        for (size_t j = 0; j < size_phi; ++j) {
            for (size_t k = 0; k < size_y; ++k) {
                fout << "   " << pT[i] << "   " << phi[j] << "   " << y_rap[k] << "   " << particles[3122].EdN_d3p_primary[i][j][k];
                for(int mu=0; mu<4; mu++)
                    fout << "   " << particles[3122].Pv_primary[i][j][k][mu];
                for(int mu=0; mu<4; mu++)
                    fout << "   " << particles[3122].Ps_primary[i][j][k][mu];
                if(decay){
                    fout << " " << particles[3122].EdN_d3p_feeddown[3212][i][j][k];
                    for(int mu=0; mu<3; mu++)
                        fout << "   " << particles[3122].Pv_feeddown[3212][i][j][k][mu];
                    for(int mu=0; mu<3; mu++)
                        fout << "   " << particles[3122].Ps_feeddown[3212][i][j][k][mu];
                    fout << " " << particles[3122].EdN_d3p_feeddown[3224][i][j][k];
                    for(int mu=0; mu<3; mu++)
                        fout << "   " << particles[3122].Pv_feeddown[3224][i][j][k][mu];
                    for(int mu=0; mu<3; mu++)
                        fout << "   " << particles[3122].Ps_feeddown[3224][i][j][k][mu];
                }
                fout << std::endl;
            }
        }
    }

    fout.close();
    
    std::cout << "The calculation is done!" << std::endl;
    return 0;
}
