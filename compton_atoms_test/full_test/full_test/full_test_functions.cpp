//
//  torus_test_functions.cpp
//  torus_test_convergence
//
//  Created by Yash Gursahani on 6/10/24.
//

#include "full_test_functions.hpp"
#include <iostream>
#include <cmath>

#include "boost/multi_array.hpp" // included so that we can read KN file into a table
#include "boost/qvm.hpp" // deals with vectors and matrices for rotation when scattering
#include "boost/random/mersenne_twister.hpp" // included so that we can generate random integers properly. Uses Matsumoto & Nishimura 1998 prescription
#include "boost/random/uniform_int_distribution.hpp" // with this, we can generate a random number in a predetermined range using the Mersenne Twister ^

// Global Variable Declarations (Constant Scalars, Vectors, Matrices)

const double rme_keV = 511.04; // electron rest-mass energy in keV
const double sigma_T = 6.652e-25; // Thomson cross section in cm^2
const double pi = 3.14159265358979323846; // pi
const double Emin = 0.2; // minimum photon energy in keV
const double Emax = 200.; // maximum photon energy in keV
const double pc_cm = 3.086e18; // 1 parsec in cm

const boost::qvm::vec<double, 3> k_hat = {0.0, 0.0, 1.0}; // k-hat unit vector for photon direction rotations

const boost::qvm::mat<double, 3, 3> I = boost::qvm::identity_mat<double, 3>(); // identity matrix for use in rotation

boost::random::mt19937 rng(time(nullptr));

// Atomic Data. Order: C, O, Ne, Mg, Si, S, Ar, Ca, Cr, Fe, Ni. Cross sections from Verner+ 1995. Abundances from Lodders+ 2009. Line energies from LBNL X-ray Data Booklet. Yields from Meddouh+ 2023.
const double E_th[11] = {0.2910E+03, 0.5380E+03, 0.8701E+03, 0.1311E+04, 0.1846E+04, 0.2477E+04, 0.3203E+04, 0.4043E+04, 0.5996E+04, 0.7124E+04, 0.8348E+04};
const double E_0[11] = {0.8655E+02, 0.1774E+03, 0.3144E+03, 0.2711E+03, 0.5322E+03, 0.8114E+03, 0.1135E+04, 0.6947E+03, 0.1035E+04, 0.8044E+03, 0.7366E+03};
const double sigma_0[11] = {0.7421E+02, 0.3237E+02, 0.1664E+02, 0.3561E+02, 0.1184E+02, 0.6649E+01, 0.4280E+01, 0.1586E+02, 0.1025E+02, 0.2055E+02, 0.2836E+02};
const double y_a[11] = {0.5498E+02, 0.3812E+03, 0.2042E+06, 0.2374E+02, 0.2580E+03, 0.3734E+04, 0.3285E+08, 0.2563E+02, 0.3343E+02, 0.3633E+02, 0.3622E+02};
const double P[11] = {0.1503E+01, 0.1083E+01, 0.8450E+00, 0.1952E+01, 0.1102E+01, 0.8646E+00, 0.7631E+00, 0.1966E+01, 0.1822E+01, 0.2118E+01, 0.2316E+01};
const double y_w[11] = {0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00, 0.0000E+00};

const double abund[11] = {2.45E-04, 5.37E-04, 1.12E-04, 3.47E-05, 3.31E-05, 1.38E-05, 3.16E-06, 2.14E-06, 4.37E-07, 2.82E-05, 1.70E-06};
const double line_E[11] = {0.277, 0.5249, 0.8486, 1.2536, 1.740, 2.3078, 2.9577, 3.6917, 5.4147, 6.4038, 7.4782};
const double yield[11] = {0.00244, 0.00594, 0.01243, 0.02333, 0.04025, 0.06477, 0.09817, 0.14110, 0.25344, 0.31925, 0.38789};

// Functions

//double draw_energy(double rand) {
//    double range = Emax - Emin;
//    return rand*range + Emin;
//}

double draw_energy(double rand) {
    double logE = rand*3.0 - 0.6989700043;
    return pow(10.0, logE);
}

double weight_energy(double E, double pl_ind) {
    return pow((E/Emin), -pl_ind);
}

int interact_medium(double n, double sigma, std::vector<double> cross_secs, double ds, double rand1, double rand2)
{
    // /////////////// dT CALCULATION SECTION ////////////////////////////
    double dtau = n*sigma*ds; // start with the optical depth from hydrogen (Compton scattering) alone
    
    double dtau_array[13]; // array with cumulative dtau intervals ranging from 0 to the total dtau. index 1 = H (Compton scattering) and indices 2-12 are the elements in the atomic data arrays
    dtau_array[0] = 0;
    dtau_array[1] = dtau; // first element of the array is dtau from H
    
    for (int i = 0; i < sizeof(sigma_0)/sizeof(sigma_0[0]); i++) {
        dtau = dtau + n*abund[i]*cross_secs[i]*ds; // add cross sections from all elements to the total optical depth
        dtau_array[i+2] = dtau;
    }
    
    // /////////////// DETERMINING PHYSICAL PROCESS ////////////////////////////
    double prob = 1.0 - exp(-dtau); // the probability that an interaction will occur
        
    if (rand1 <= prob) {
        // interaction!
        for (int i = 0; i < sizeof(dtau_array)/sizeof(dtau_array[0]) - 1; i++) {
            if (rand2 >= dtau_array[i]/dtau && rand2 <= dtau_array[i+1]/dtau) {
                return i-1; // returns -1 for Compton otherwise 0-10 for all elements
            }
        }
    }
    
    return -2; // no interaction :(
}

boost::qvm::vec<double, 3> choose_direction(double rand_x, double rand_y, double rand_z) {
    double x = -1.0 + (rand_x * 2.0); // this puts the x direction between -1 and 1
    double y = -1.0 + (rand_y * 2.0); // this puts the y direction between -1 and 1
    double z = -1.0 + (rand_z * 2.0); // this puts the z direction between -1 and 1
    
    boost::qvm::vec<double, 3> u = {x, y, z};
    return u/boost::qvm::mag(u);
}

double choose_sigma(double E)
{
    if (E < 10.0) {return sigma_T;} // Thomson sigma for E < 10 keV, implemented by GF91
    
    else {
        double x = E/rme_keV;
        double f = 1 + 2*x;
        return 0.75*sigma_T*(((1+x)/pow(x, 3.0))*((2*x*(1+x)/f) - log(f)) + log(f)/(2*x) - (1+3*x)/pow(f, 2.0)); // full K-N sigma for E >= 10 keV
    }
}

double choose_rand_angle(double rand) {
    return 2.0*pi*rand;
}

std::vector<double> calc_sigmas(double E) {
    double E_eV = E*1000.; // put photon energy in eV to match atomic data
    
    std::vector<double> cross_secs; // will store every cross section in this vector before returning it
    
    for (int i = 0; i < sizeof(sigma_0)/sizeof(sigma_0[0]); i++) {
        
        if (E_eV < E_th[i]) {
            cross_secs.push_back(0.0); // zero cross section if out of energy range
        }
        
        else {
            double y = E_eV/E_0[i];
            double Q = 5.5 - 0.5*P[i];
            
            double A = pow(y-1., 2.0) + pow(y_w[i], 2.0);
            double B = pow(y, -Q);
            double C = pow(1 + sqrt(y/y_a[i]), -P[i]);
            
            cross_secs.push_back(sigma_0[i]*A*B*C*1e-18); // calculate cross section and convert from Mb to cm^2
        }
        
    }
    
    return cross_secs;
}


double new_energy(double E, double theta) {
    double x = E/rme_keV;
    return E/(1 + x*(1 - cos(theta))); // Compton formula for energy shift
}

double fe_energy(double rand1, double rand2) {
    if (rand1 < 0.04) {return 7.056;} // Fe Kb line
    else {
        if (rand2 < 0.666) {return 6.404;} // Fe Ka1 line
        else {return 6.391;} // Fe Ka2 line
    }
}

boost::qvm::vec<double, 3> rotate(boost::qvm::vec<double, 3> vel, double theta, double phi) {
    boost::qvm::vec<double, 3> w_hat = boost::qvm::cross(k_hat, vel); // k x v
    double w1 = X(w_hat);
    double w2 = Y(w_hat);
    double w3 = Z(w_hat);
    
    double c = boost::qvm::dot(k_hat, vel); // k . v
    
    boost::qvm::vec<double, 3> scat_vel = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)}; // scattering relative to k-hat
    
    boost::qvm::mat<double, 3, 3> V = {0, -w3, w2, w3, 0, -w1, -w2, w1, 0}; // V matrix from Rodrigues Rotation Formula
    
    boost::qvm::mat<double, 3, 3> R = I + V + (1.0/(1.0 + c))*(V*V); // Rotation matrix from k-hat to pre-scattering direction
    
    boost::qvm::vec<double, 3> new_vel = R*scat_vel; // scat_vel rotated by R into new_vel
    return new_vel;
}

bool escape_1(boost::qvm::vec<double, 3> vel_cart, double angle) {
    double v_xy = sqrt(pow(X(vel_cart), 2.0) + pow(Y(vel_cart), 2.0)); // horizontal component in the plane
    double v_z = Z(vel_cart); // vertical component (in place or otherwise)
    double test_angle = abs(atan(v_xy/v_z));
    
    if (test_angle < angle) {return true;} // if initial direction allows photon to leave, return true.
    else {return false;} // else, return false.
}

bool escape_2(boost::qvm::vec<double, 3> pos_cart, boost::qvm::vec<double, 3> vel_cart, double R, double r, double half_oa, double ds) {
    double x = X(pos_cart);
    double y = Y(pos_cart);
    double z = Z(pos_cart);
    
    double vx = X(vel_cart);
    double vy = Y(vel_cart);
    double vz = Z(vel_cart);
    
    double R_c = R*sin(half_oa);
    double z_0 = R*cos(half_oa);
    
    if (vz == 0) {return false;} // photon cannot escape if it is not moving vertically
    
    if (vz < 0) {z_0 = z_0*-1.0;} // want to test for lower boundary if photon is moving downward
    
    double delta_t = (z_0 - z)/vz;
    double x_0 = x + vx*delta_t;
    double y_0 = y + vy*delta_t;
    
    if (pow(x_0, 2.0) + pow(y_0, 2.0) < pow(R_c, 2.0)) {return true;} // escape if the photon will go through the circular cap at the top/bottom of the hollow region!
    
    else {return false;}
}

bool escape_3(boost::qvm::vec<double, 3> pos_cart, double R) {
    double x = X(pos_cart);
    double y = Y(pos_cart);
    double z = Z(pos_cart);
    
    double r_pos = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
    
    if (r_pos <= R) {return false;} // if still in the large sphere, return false
    else {return true;} // if out of the large sphere, return true -- will move to next photon
}

bool check_medium(boost::qvm::vec<double, 3> pos_cart, double r, double f) {
//    double x = X(pos_cart);
//    double y = Y(pos_cart); replace with norm in r_pos below. See if this speeds things up.
    double z = Z(pos_cart);
    
    double r_pos = boost::qvm::mag(pos_cart); // current r coordinate of the photon in spherical coordinates
    double theta = acos(abs(z)/r_pos); // current theta angle in speherical coordinates (w.r.t. +z axis) of the photon
    double psi = (pi/2.0) - theta; // the angle we care about when calculating whether the photon is in the medium or not
    
    double d = r*f*pow(((sin(psi)*sqrt(1 - pow(f, 2.0))) + (cos(psi)*f)), -1.0);
    
    if (r_pos < r) {return false;} // not in the medium if inside the inner radius; inside the bicone-like region in the middle of the torus
    else if (r_pos >= r && r_pos >= d) {return true;} // in the medium if beyond the threshold value d calculated for this geometry
    else {return false;} // not in the medium, but inside the bicone-like region in the middle of the torus
}

int int_energy(double E) {
    int M = 2000;
    
    double A = (1./M)*log(Emax/Emin);
    double i = (1.0/A)*log(E/Emin);
    //std::cout << E << " " << i << std::endl;
    int ind = std::round(i);
    return ind;
}

size_t sizet_theta(double rand) {
    double rand_1000 = rand*1000.0;
    size_t ind = std::round(rand_1000);
    return ind;
}

std::vector<double> propagate_photon(boost::multi_array<double, 2> table, double PL_ind, double R, double r, double half_oa, double f, double ds, double n) {
    boost::random::uniform_int_distribution<> dist_theta(0,1000); // used when we want to generate a random index for our scattering angle lookup table
    
    double phot_E_rand = (double)rng()/(double)rng.max(); // random number to choose the initial photon energy
    double phot_E = draw_energy(phot_E_rand); // initial photon energy    
    double weight = weight_energy(phot_E, PL_ind); // store weight of this photon
    double phot_E_initial = phot_E; // store initial energy
    
    double scat_counter = 0.0; // how many scatterings have occurred for this photon
    
    std::vector<double> cross_sections = calc_sigmas(phot_E); // all atomic cross sections for metals, will be updated every time a the photon changes energy (Compton scattering or line emission)
    
    double rand_x = (double)rng()/(double)rng.max(); // random number for choosing x direction
    double rand_y = (double)rng()/(double)rng.max(); // random number for choosing y direction
    double rand_z = (double)rng()/(double)rng.max(); // random number for choosing z direction
    
    boost::qvm::vec<double, 3> phot_pos = {0.0, 0.0, 0.0}; // cm, position of the photon
    boost::qvm::vec<double, 3> phot_vel = choose_direction(rand_x, rand_y, rand_z); // unitless direction vector for photon
    
    if (escape_1(phot_vel, half_oa)) {
        double pos_X = X(phot_pos);
        double pos_Y = Y(phot_pos);
        double pos_Z = Z(phot_pos);
        
        double vel_X = X(phot_vel);
        double vel_Y = Y(phot_vel);
        double vel_Z = Z(phot_vel);
        
        double num_scatters = scat_counter;
        double phot_E_final = phot_E;
        
        std::vector<double> phot_properties {pos_X, pos_Y, pos_Z, vel_X, vel_Y, vel_Z, num_scatters, phot_E_initial, phot_E_final, weight};
        return phot_properties;
        
        // immediately returns properties and exits once escape condition #1 is met: instant escape from the source
    }
    
    while (escape_3(phot_pos, R) == false) // Photon propagation and atomic physics loop
    {
        
        phot_pos = phot_pos + (ds*phot_vel); // propagate photon by one step length in current direction
        
        if (check_medium(phot_pos, r, f)) {
            
            double sigma = choose_sigma(phot_E);
            
            double interact_rand = (double)rng()/(double)rng.max(); // random number which determines whether we interact with medium or not
            double process_rand = (double)rng()/(double)rng.max(); // random number which determines which process the photon undergoes if it interacts with medium
            
            int process = interact_medium(n, sigma, cross_sections, ds, interact_rand, process_rand); // do we interact with medium or not?
            
            if (process == -1) // if process == -1, we scatter!
            {
                scat_counter = scat_counter + 1.0;
                
                double phi_rand = (double)rng()/(double)rng.max(); // random number with which we draw a new phi
                
                int k = int_energy(phot_E);
                int l = dist_theta(rng);
                
                if (k < 0) {k = 0;} // since we are scattering a lot of low-E photons, it's possible they will get scattered to just below Emin. KN table cannot handle that, so we fix the lowest possible E as Emin.
                
                double theta = table[k][l]; // choose new theta
                double phi = choose_rand_angle(phi_rand); // choose new phi
                
                phot_vel = rotate(phot_vel, theta, phi); // rotate our photon's direction
                phot_E = new_energy(phot_E, theta); // update photon energy due to Compton scattering
                cross_sections = calc_sigmas(phot_E); // update all the cross sections for metals
            }
            
            if (process >= 0) // if process >= 0, fluorescence or absorption!
            {
                double atomic_rand = (double)rng()/(double)rng.max(); // random number with which we decide on the atomic process
                
                if (atomic_rand < yield[process]) // fluorescence!
                {
                    if (process == 9) { // special case for iron branching into 3 possible lines
                        double rand_Ka_Kb = (double)rng()/(double)rng.max();
                        double rand_Ka1_Ka2 = (double)rng()/(double)rng.max();
                        phot_E = fe_energy(rand_Ka_Kb, rand_Ka1_Ka2);
                    }
                    
                    else {phot_E = line_E[process];} // fluorescent line energy
                    
                    double theta_rand = (double)rng()/(double)rng.max();
                    double phi_rand = (double)rng()/(double)rng.max();
                    
                    double theta = choose_rand_angle(theta_rand); // randomly choose new theta
                    double phi = choose_rand_angle(phi_rand); // randomly choose new phi
                    
                    phot_vel = rotate(phot_vel, theta, phi); // rotate our photon's direction
                    cross_sections = calc_sigmas(phot_E); // update all the cross sections for metals
                }
                
                else { // absorption!
                    phot_E = 0.0; // a flag that tells us the photon has been absorbed
                    break; // break out of propagation loop and return values
                }
            }
        }
        
        else if (escape_2(phot_pos, phot_vel, R, r, half_oa, ds)) {break;} // escape condition #2 has been satisfied
   
    }

    double pos_X = X(phot_pos);
    double pos_Y = Y(phot_pos);
    double pos_Z = Z(phot_pos);
    
    double vel_X = X(phot_vel);
    double vel_Y = Y(phot_vel);
    double vel_Z = Z(phot_vel);
    
    double num_scatters = scat_counter;
    double phot_E_final = phot_E;
    
    std::vector<double> phot_properties {pos_X, pos_Y, pos_Z, vel_X, vel_Y, vel_Z, num_scatters, phot_E_initial, phot_E_final, weight};
    return phot_properties;
    // either escape condition #2 (escape from cavity) or condition #3 (photon out of geometry entirely) has been triggered, OR a photon has been absorbed and the properties are returned
}
