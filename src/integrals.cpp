#include <fstream>
#include <iostream>

#ifdef OPEN_MP
    #include<omp.h>
#endif

#include "utils.h"
#include "integrals.h"

const std::vector<NonZeroLevi> nonZeroLeviElements = getNonZeroLevi();

double momentumCM(double s, double m1, double m2) {
    if (s < (m1 + m2) * (m1 + m2)) {
        std::cerr << "Error: sqrt of negative number. Ensure s is large enough." << std::endl;
        exit(1);
    }
    
    double term = (s - (m1 + m2) * (m1 + m2)) * (s - (m1 - m2) * (m1 - m2));
    return sqrt(term) / (2.0 * sqrt(s));
}

std::array<double, 3> ResonanceThreeMomentum(const double MR, const double E, const double Es, const std::array<double, 3>& p, const std::array<double, 3>& ps) {
    std::array<double, 3> pdiff, PR;

    for (size_t mu = 0; mu < 3; ++mu) {
        pdiff[mu] = p[mu] - ps[mu];
    }

    const double pnorm = dotProduct(pdiff, pdiff);
    const double denom = (Es + E) * (Es + E) - pnorm;
    if (denom == 0) {
        std::cerr << "Error: Division by zero in ResonanceThreeMomentum." << std::endl;
        exit(1);
    }

    for (size_t mu = 0; mu < 3; ++mu) {
        PR[mu] = 2 * MR * (Es + E) * pdiff[mu] / denom;
    }

    return PR;
}

double Jacobian(const double Mm, const double E, const std::array<double,3>& p, const double Mp, const double Es, const std::array<double,3>& ps) {
    double pps = dotProduct(p,ps);
    double num = pow(Mm,3)*(Es+E)*(Es+E)*((Es+E)*(Es+E)- (E*Es + pps + Mp*Mp));
    double den = Es*pow(E*Es + pps + Mp*Mp,3);
    return num/den;
}

std::array<double,3> RestFrameSpinVector(const double mass, const std::array<double,3>& p, const std::array<double,3>& S) {
    std::array<double, 3> Sstar = {0.,0.,0.};
    double en = sqrt(dotProduct(p,p)+mass*mass);
    double pS = dotProduct(p,S);
    
    for (size_t mu = 0; mu < 3; ++mu) {
        Sstar[mu] = S[mu] - p[mu]*pS/(en*(en+mass));
    }
    return Sstar;
}

double aux_exact_polarization(double spin, double pu, double T, double mutot, double abs_theta) {
    double num = 0;
    double den = 1e-20;
    int fermi_or_bose = statistics(spin);
    
    if (fermi_or_bose == 0) {
        std::cerr << "Error: Unknown statistics (fermi_or_bose == 0)!" << std::endl;
        exit(1);
    }

    double exp_factor = (pu - mutot) / T;
    
    for (double k = -spin; k <= spin; ++k) {
        double denominator = exp(exp_factor - k * abs_theta) + fermi_or_bose;
        num += k / denominator;
        den += 1 / denominator;
    }

    if (den == 0) {
        std::cerr << "Denominator is zero in aux_exact_polarization!" << std::endl;
        exit(1);
    }

    return num / den;
}

void sum_over_surface(const std::vector<element> &freeze_out_sup, const Particle& particle, const std::array<double, 4> p, double& dndp, double P_vorticity[4], double P_shear[4]){
    dndp = 0.0;
    for (int i = 0; i < 4; ++i) {
        P_vorticity[i] = 0.0;
        P_shear[i] = 0.0;
    }
    const std::array<double, 4> p_ = {p[0], -p[1], -p[2], -p[3]}; // p_ is covariant
    const double mass = particle.get_mass();

    #ifdef OPEN_MP
    #pragma omp parallel for reduction(+:P_vorticity[:4], P_shear[:4], dndp)
    #endif
    for(size_t i = 0; i < freeze_out_sup.size(); ++i){ //loop over the FO hypersurface
        const element& cell = freeze_out_sup[i];
        double pdSigma = p[0] * cell.dsigma[0] + p[1] * cell.dsigma[1] + p[2] * cell.dsigma[2] + p[3] * cell.dsigma[3]; // GeV fm3
        double pu = p[0] * cell.u[0] - p[1] * cell.u[1] - p[2] * cell.u[2] - p[3] * cell.u[3]; // GeV

        const double mutot = cell.mub*particle.get_baryon() + cell.muq*particle.get_charge() + cell.mus*particle.get_strange();
	const double nf = 1.0 / (exp( (pu - mutot) / cell.T) + 1.0);
	
        dndp += pdSigma * nf ;
        //variables to the sum involving Levi Civita tensor. The indices in all expressions have been renamed so that Levi Civita tensor has indices in the order 
        
        double S1_thw[4] = {0.,0.,0.,0.};
        double S1_ths[4] = {0.,0.,0.,0.};
        
        for (const auto& element : nonZeroLeviElements) {
            int mu = element.mu;
            int nu = element.nu;
            int rh = element.rh;
            int sg = element.sg;
            double Levi = element.value; // Can be only 1 or -1. 0 has been discarded earlier to avoid unnecessary loops
            
            S1_thw[mu] += Levi * p_[sg] * cell.dbeta[nu][rh]; // epsilon^{mu nu rho sigma} p_{sigma} d_nu beta_rho. Eq(40) fm^{-1}
            
            if (nu == 0) {
                for (int ta = 0; ta < 4; ta++) {
                    S1_ths[mu] += -Levi * p_[sg] * p[ta] * (cell.dbeta[rh][ta] + cell.dbeta[ta][rh]);// epsilon^{mu nu rho sigma} t_nu p_{sigma} p^{tau} xi_{rho tau}. Eq(40) GeV fm^{-1}
		}
            }
        }
        
        for (int mu = 0; mu < 4; ++mu) {
            P_vorticity[mu] += pdSigma * nf * (1. - nf) * S1_thw[mu]; // GeV fm^{2}
            P_shear[mu] += pdSigma * nf * (1. - nf) * S1_ths[mu] / p[0]; // GeV fm^{2}
        }
    }
    
    for (int mu = 0; mu < 4; ++mu) {
        P_vorticity[mu] = P_vorticity[mu]/ (8.0 * mass) *hbarC; // GeV fm^3
        P_shear[mu] = P_shear[mu]/ (8.0 * mass) *hbarC; // GeV fm^3
    }
}

void sum_over_surface_exact(const std::vector<element> &freeze_out_sup, const Particle& particle, const std::array<double, 4> p, double& dndp, double P_vorticity[4], double P_shear[4]){
    dndp = 0.0;
    for (int i = 0; i < 4; ++i) {
        P_vorticity[i] = 0.0;
        P_shear[i] = 0.0;
    }
    double dndp_=0., P_vorticity_[4]={0.,0.,0.,0.}, P_shear_[4]={0.,0.,0.,};
    const std::array<double, 4> p_ = {p[0], -p[1], -p[2], -p[3]}; // p_ is covariant
    const double spin = (particle.get_gspin()-1.0)/2.0;
    const int fermi_or_bose = statistics(spin);
    const double phase_space = ((2*spin +1)/pow( 2*PI, 3.0)); //dimensionless
    const double mass = particle.get_mass();

    #ifdef OPEN_MP
    #pragma omp parallel for reduction(+:P_vorticity[:4], P_shear[:4], dndp)
    #endif
    for(size_t i = 0; i < freeze_out_sup.size(); i++){ //loop over the FO hypersurface
        const element& cell = freeze_out_sup[i];
        double pdSigma = p[0] * cell.dsigma[0] + p[1] * cell.dsigma[1] + p[2] * cell.dsigma[2] + p[3] * cell.dsigma[3]; // GeV fm3
        double pu = p[0] * cell.u[0] - p[1] * cell.u[1] - p[2] * cell.u[2] - p[3] * cell.u[3]; // GeV
        const double mutot = cell.mub*particle.get_baryon() + cell.muq*particle.get_charge() + cell.mus*particle.get_strange();
        const double dist = 1.0 / (exp( (pu - mutot) / cell.T) + fermi_or_bose);
	
        dndp += phase_space * pdSigma * dist ; // GeV fm^3
        
        //variables to the sum involving Levi Civita tensor. The indices in all expressions have been renamed so that Levi Civita tensor has indices in the order 
        
        std::array<double,4> theta_vector = {0,0,0,0};
        double S1_ths[4] = {0.,0.,0.,0.};
        
        for (const auto& element : nonZeroLeviElements) {
            int mu = element.mu;
            int nu = element.nu;
            int rh = element.rh;
            int sg = element.sg;
            double Levi = element.value; // Can be only 1 or -1. 0 has been discarded earlier to avoid unnecessary loops
            
            theta_vector[mu] += Levi*p_[sg]*cell.dbeta[nu][rh]*hbarC/(2.0*mass); //dimensionless
            
            if (nu == 0) {
                for (int ta = 0; ta < 4; ta++) {
                    S1_ths[mu] += -(spin/3.0) * (spin + 1) * Levi * (p_[sg]/mass) * (p[ta]/p[0]) * (cell.dbeta[rh][ta] + cell.dbeta[ta][rh]) * hbarC / 2.0;// epsilon^{mu nu rho sigma} t_nu p_{sigma} p^{tau} xi_{rho tau}. Eq(40) (dimensionless)
		}
            }
        }
        
        double theta_sq = theta_vector[0]*theta_vector[0] - theta_vector[1]*theta_vector[1] - theta_vector[2]*theta_vector[2] - theta_vector[3]*theta_vector[3];
        
        for (int mu = 0; mu < 4; ++mu) {
            P_vorticity[mu] += phase_space * pdSigma * dist * (theta_vector[mu]/sqrt(-theta_sq)) * aux_exact_polarization(spin, pu, cell.T, mutot, sqrt(-theta_sq));  // GeV fm^3
            P_shear[mu] += phase_space * pdSigma * dist * (1. - fermi_or_bose*dist) * S1_ths[mu]; // GeV fm^3
        }
    }
}

void compute_primary(const std::vector<double>& pT_vec, const std::vector<double>& phi_vec, const std::vector<double>& y_vec, std::map<int, Particle>& particles, std::vector<element>& freeze_out_sup, const bool exact, const bool decay, int pdg_id) {

    if(decay){
    for (auto& [p_id, particle] : particles) {
        if (p_id == 22){
            std::cout << "Skipping calculation for photon." << std::endl;
            continue;
        }
        for (size_t i = 0; i < pT_vec.size(); ++i) {
            for (size_t j = 0; j < phi_vec.size(); ++j) {
                for (size_t k = 0; k < y_vec.size(); ++k) {
                    
                    double pT = pT_vec[i];
                    double phi = phi_vec[j];
                    double y = y_vec[k];

                    double P_vorticity[4] = {0., 0., 0., 0.};
                    double P_shear[4] = {0., 0., 0., 0.};
                    double dndp = 0.;

                    const double mass = particle.get_mass();  // GeV
                    double mT = sqrt(mass * mass + pT * pT); // GeV
                    
                    // Define 4-momentum. p is contravariant
                    std::array<double, 4> p = {mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)};
                    
                    if(exact){
                        sum_over_surface_exact(freeze_out_sup, particle, p, dndp, P_vorticity, P_shear);
                    }else{
                        sum_over_surface(freeze_out_sup, particle, p, dndp, P_vorticity, P_shear);
                    }

                    particle.EdN_d3p_primary[i][j][k] = dndp;
                    for (int mu = 0; mu < 4; ++mu) {
                        particle.Pv_primary[i][j][k][mu] = P_vorticity[mu]; // GeV fm^3
                    }
                    for (int mu = 0; mu < 4; ++mu) {
                        particle.Ps_primary[i][j][k][mu] = P_shear[mu]; // GeV fm^3
                    }
                }
            }
        }
    }
    }else{
        Particle& particle = particles[pdg_id];
        for (size_t i = 0; i < pT_vec.size(); ++i) {
            for (size_t j = 0; j < phi_vec.size(); ++j) {
                for (size_t k = 0; k < y_vec.size(); ++k) {
                    
                    double pT = pT_vec[i];
                    double phi = phi_vec[j];
                    double y = y_vec[k];

                    double P_vorticity[4] = {0., 0., 0., 0.};
                    double P_shear[4] = {0., 0., 0., 0.};
                    double dndp = 0.;

                    const double mass = particle.get_mass();  // GeV
                    double mT = sqrt(mass * mass + pT * pT); // GeV
                    
                    // Define 4-momentum. p is contravariant
                    std::array<double, 4> p = {mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)};
                    
                    if(exact){
                        sum_over_surface_exact(freeze_out_sup, particle, p, dndp, P_vorticity, P_shear);
                    }else{
                        sum_over_surface(freeze_out_sup, particle, p, dndp, P_vorticity, P_shear);
                    }

                    particle.EdN_d3p_primary[i][j][k] = dndp;
                    for (int mu = 0; mu < 4; ++mu) {
                        particle.Pv_primary[i][j][k][mu] = P_vorticity[mu]; // GeV fm^3
                    }
                    for (int mu = 0; mu < 4; ++mu) {
                        particle.Ps_primary[i][j][k][mu] = P_shear[mu]; // GeV fm^3
                    }
                }
            }
        }    
    }
}

void compute_polarization_feeddown(const std::vector<double>& pT_vec, const std::vector<double>& phi_vec, const std::vector<double>& y_vec, Particle& primary, std::map<int, Particle>& particles) {
    
    trilinearInterpolator interpolator(pT_vec, phi_vec, y_vec);
    
    const int pdg_id = primary.pdg_id;
    const int size_th = 21;
    const int size_phi = 41;
    
    std::vector<double> thstar = linspace(0.0,PI,size_th);
    std::vector<double> phistar =  linspace(0.0,2*PI,size_phi);
    
    for (const auto& [mother_id, mother] : particles) {
        if (mother.get_mass() <= primary.get_mass()) {
            continue;
        }
        
        const auto& dNdP_mother_primary = mother.get_dndp_primary();
        const auto& Sv_mother_primary = mother.get_Pv_primary();
        const auto& Ss_mother_primary = mother.get_Ps_primary();

        for (const auto& decay : mother.decay_channels) {
            if (std::find(decay.daughters.begin(), decay.daughters.end(), primary.pdg_id) == decay.daughters.end()) {
                continue;
            }
            
            if (decay.num_daughters > 2){
                std::cout << "Skipping calculation for 3 or more body decay of mother particle " << mother_id << std::endl;
                continue;
            }
            
            int other_daughters[1];
            for (int daughter : decay.daughters) {
                if (daughter != primary.pdg_id) {
                    other_daughters[0] = daughter;
                }
            }
            
            const Particle& d2 = particles.at(other_daughters[0]);
    
            for (size_t i = 0; i < pT_vec.size(); ++i) {
                double pT = pT_vec[i];
                double mT = sqrt(primary.mass * primary.mass + pT * pT); // GeV
                for (size_t j = 0; j < phi_vec.size(); ++j) {
                    double phi = phi_vec[j];
                    for (size_t k = 0; k < y_vec.size(); ++k) {
                        double y = y_vec[k];
                        // particle 4-momentum in lab frame
                        std::array<double, 4> p = {mT * cosh(y), pT * cos(phi), pT * sin(phi), mT * sinh(y)};
                        // extract spatial part of 4-momentum above
                        std::array<double, 3> p3 = {p[1],p[2],p[3]};
                        // Integrate over solid angle = sin(th) dth dphi. Hence two loops
                        double num_eq29_vo[3] = {0., 0., 0.};
                        double num_eq29_sh[3] = {0., 0., 0.};
                        double den_eq29 = 1e-20;
                        #ifdef OPEN_MP
                        #pragma omp parallel for reduction(+:num_eq29_vo[:3], num_eq29_sh[:3], den_eq29) collapse(2)
                        #endif
                        for(size_t ith = 0; ith < size_th; ith++){
                            for (size_t iph = 0; iph < size_phi; ++iph) {
                                //particle 3-momentum and energy in mother's rest frame
                                std::array<double, 3> pstar;
                                double pstar_norm = momentumCM(mother.mass*mother.mass, primary.mass, d2.mass); // magnitude particle 3-momentum
                                double Estar = sqrt(pstar_norm*pstar_norm + primary.mass*primary.mass);//particle energy
                                pstar[0] = pstar_norm*sin(thstar[ith])*cos(phistar[iph]); //p^x
                                pstar[1] = pstar_norm*sin(thstar[ith])*sin(phistar[iph]); //p^y
                                pstar[2] = pstar_norm*cos(thstar[ith]); //p^z
                                //mother's energy, momentum, pT, phi, y in lab frame
                                std::array<double, 3> P_mother;
                                P_mother = ResonanceThreeMomentum(mother.mass, p[0], Estar, p3, pstar);
                                double E_mother = sqrt(dotProduct(P_mother,P_mother) + mother.mass*mother.mass);
                                double mother_pT = sqrt(P_mother[0]*P_mother[0]+P_mother[1]*P_mother[1]);
                                double mother_phi = atan2(P_mother[1],P_mother[0]);
                                if(mother_phi < 0) mother_phi += 2*PI;
                                double mother_y = atanh(P_mother[2]/E_mother);
                                
                                //interpolation of quantities in lab frame
                                // EdN/d^3p
                                double dNdP_mother = interpolator.trilinearInterpolation(mother_pT, mother_phi, mother_y, dNdP_mother_primary);
                                //polarization
                                std::vector<double> Sv_mother_lab(4,0.0), Ss_mother_lab(4,0.0);
                                Sv_mother_lab = interpolator.trilinearInterpolation(mother_pT, mother_phi, mother_y, Sv_mother_primary);
                                Ss_mother_lab = interpolator.trilinearInterpolation(mother_pT, mother_phi, mother_y, Ss_mother_primary);

                                // Remember: the primary computation yield numerator (surface integ of vorticity or shear) and denominator (yield) of spin vector separately. 
                                // Mean spin vector is numerator/denominator and is done in plotting script. Here, we must divide as the mean goes into integrand.
                                for (size_t mu = 0; mu<4; ++mu){
                                    if(dNdP_mother != 0){
                                        Sv_mother_lab[mu] = Sv_mother_lab[mu]/dNdP_mother;
					Ss_mother_lab[mu] = Ss_mother_lab[mu]/dNdP_mother;
				    }else{
					Sv_mother_lab[mu] = 0.0;
					Ss_mother_lab[mu] = 0.0;
				    }
				}
                                
                                //mother rest frame polarization
                                std::array<double, 3> S3_mother_lab = {Sv_mother_lab[1],Sv_mother_lab[2],Sv_mother_lab[3]};
                                std::array<double, 3> Sv_mother_RF = {0.,0.,0.}, Ss_mother_RF = {0.,0.,0.};
                                Sv_mother_RF = RestFrameSpinVector(mother.mass, P_mother, S3_mother_lab);
                                
                                S3_mother_lab[0] = Ss_mother_lab[1];
                                S3_mother_lab[1] = Ss_mother_lab[2];
                                S3_mother_lab[2] = Ss_mother_lab[3];
                                Ss_mother_RF = RestFrameSpinVector(mother.mass, P_mother, S3_mother_lab);
                                
                                //solid angle integrand
				double jac = Jacobian(mother.mass, p[0], p3, primary.mass, Estar, pstar);
                                double spv = dotProduct(Sv_mother_RF,pstar);
                                double sps = dotProduct(Ss_mother_RF,pstar);

                                den_eq29 += sin(thstar[ith])*jac*(dNdP_mother/E_mother);
                                switch(mother.pdg_id){
                                    case 3212: //Sigma0 to Lambda and photon
                                    for(int mu=0;mu<3;mu++){
                                        num_eq29_vo[mu] += -sin(thstar[ith])*(dNdP_mother/E_mother)*jac*spv*pstar[mu]/(pstar_norm*pstar_norm);
                                        num_eq29_sh[mu] += -sin(thstar[ith])*(dNdP_mother/E_mother)*jac*sps*pstar[mu]/(pstar_norm*pstar_norm);
                                    }
                                    break;
                                    case 3224:
                                    for(int mu=0;mu<3;mu++){
                                        num_eq29_vo[mu] += sin(thstar[ith])*(dNdP_mother/E_mother)*jac*0.4*(Sv_mother_RF[mu]-0.5*spv*pstar[mu]/(pstar_norm*pstar_norm));
                                        num_eq29_sh[mu] += sin(thstar[ith])*(dNdP_mother/E_mother)*jac*0.4*(Ss_mother_RF[mu]-0.5*sps*pstar[mu]/(pstar_norm*pstar_norm));
                                    }
                                    break;
                                    default:
                                        std::cout << "The current implementation considers feeddown from 3212 and 3224 only" << std::endl;
                                        exit(1);
                                }
                            }
                        }
                        
                        primary.EdN_d3p_feeddown[mother_id][i][j][k] = den_eq29;
                        for (int mu = 0; mu < 3; mu++) {
                            primary.Pv_feeddown[mother_id][i][j][k][mu] = num_eq29_vo[mu];
                        }
                        for (int mu = 0; mu < 3; mu++) {
                            primary.Ps_feeddown[mother_id][i][j][k][mu] = num_eq29_sh[mu];
                        }
                    }
                }
            }
        }
    }
}
