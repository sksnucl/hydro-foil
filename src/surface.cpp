#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>

#include "surface.h"
#include "utils.h"


using namespace std;

void read_hypersurface(string filename, vector<element> &hypersurface){
    ifstream input_file(filename);
    if(!input_file.is_open()){
        cout << "Failed to open " << filename << endl;
        exit(1);
    }
    cout << "Reading hypersurface from " << filename << endl;

    string line;
    while (getline(input_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // skip empty or comment lines
        }

        std::istringstream iss(line);
        element cell;
        double dummy, vx, vy, vz;
        double dtemp[4];
        double dmu[4][4];
        
        iss >> cell.tau >> cell.x >> cell.y >> cell.eta;
        
        for(int mu=0;mu<4;mu++){
            iss >> cell.dsigma[mu];
        }
        
        for(int mu=0;mu<4;mu++){
            iss >> cell.u[mu];
        }
        
        iss >> cell.T >> cell.mub >> cell.muq >> cell.mus;

        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                iss >> cell.dbeta[i][j];
            }
        }

        iss >> dummy;
                
        hypersurface.push_back(cell);
    }
    cout << "Reading successful!" << endl;
}
