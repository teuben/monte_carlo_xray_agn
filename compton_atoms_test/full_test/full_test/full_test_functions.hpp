//
//  torus_test_functions.hpp
//  torus_test_convergence
//
//  Created by Yash Gursahani on 6/10/24.
//

#ifndef full_test_functions_hpp
#define full_test_functions_hpp

#include <stdio.h>
#include <iostream>

#include "boost/multi_array.hpp" // included so that we can read KN file into a table
#include "boost/qvm.hpp" // Boost QVM library for vectors and matrices
#include "boost/random/mersenne_twister.hpp" // included so that we can generate random integers properly. Uses Matsumoto & Nishimura 1998 prescription

// Global Variable Declarations (Constants)
extern const double rme_keV; // electron rest-mass energy in keV
extern const double sigma_T; // Thomson cross section in cm^2
extern const double pi; // pi
extern const double Emin; // minimum photon energy in keV
extern const double Emax; // maximum photon energy in keV
extern const double pc_cm; // 1 parsec in centimeters

// Atomic Data. Order: C, O, Ne, Mg, Si, S Ar, Ca, Cr, Fe, Ni
extern const double E_th[11]; // eV
extern const double E_max[11]; // eV
extern const double E_0[11]; // eV
extern const double sigma_0[11];
extern const double y_a[11];
extern const double P[11];
extern const double y_w[11];
extern const double y_0[11];
extern const double y_1[11];

extern const double abund[11];
extern const double line_E[11]; // keV
extern const double yield[11];

// Functions

double draw_energy(double rand);
// rand: a randomly generated double between 0 and 1
// returns: an initial photon energy between 0.2 and 200 keV, drawn from a uniform distribution

double weight_energy(double E, double pl_ind);
// E: the initial photon energy (keV)
// pl_ind: the instrinsic AGN power-law index (in terms of energy)
// returns: the weighting (normalized to 1) of this photon when it goes into the spectrum later, according to the power-law

int interact_medium(double n, double sigma, std::vector<double> cross_secs, double ds, double rand1, double rand2);
// n: the number density of scatterers (electrons) in the medium
// sigma: the total cross section of the Compton scatterer (H) at the current photon energy
// cross_secs: the vector defining the cross sections of all metals at the current photon energy
// ds: the step length
// rand1: a randomly generated double between 0 and 1 to determine if an interaction occurs
// rand2: a randomly generated double between 0 and 1 that will only be used if an interaction occurs, to determine which interaction it is (Compton Scattering, fluorescent line emission, or Auger ejection/absorption)
// returns: -2 for no interaction, -1 for Compton scattering, [0-10] for any atomic process corresponding to the index in the atomic data arrays

boost::qvm::vec<double, 3> choose_direction(double rand_x, double rand_y, double rand_z);
// rand_x, rand_y, rand_z: randomly generated doubles between 0 and 1
// returns: a unit vector that initializes the direction for the photon

double choose_sigma(double E);
// E: the current photon energy (keV)
// returns: the TOTAL scattering cross section for the given energy regime (0.2-10 keV Thomson, > 10 keV full KN formula)

std::vector<double> calc_sigmas(double E);
// E: the current photon energy (keV)
// returns: a vector containing all the cross sections for neutral metals at a given energy (cm^2)

double choose_rand_angle(double rand);
// rand: a randomly generated double between 0 and 1
// returns: a new RANDOM angle for the photon, selected isotropically

double fe_energy(double rand1, double rand2);
// rand1: a randomly generated double between 0 and 1 that branches between Fe Ka and Kb
// rand2: a randomly generated double between 0 and 1 that branches between Fe Ka1 and Ka2
// returns: Fe fluorescent energy for ONE of either Fe Ka1, Fe Ka2, or Fe Kb

double new_energy(double E, double theta);
// E: the current photon energy (keV)
// theta: the scattering angle
// returns: the new photon energy based on Compton Scattering

int int_energy(double E);
// E: the current photon energy (keV)
// returns: an integer between 0 and 2000, which indicates which energy to use in our theta lookup table

size_t sizet_theta(double rand);
// rand: a randomly generated double between 0 and 1
// returns: an integer between 0 and 1000, which indicates which theta to choose in our theta lookup table. Should be used together with int_energy function above

bool escape_1(boost::qvm::vec<double, 3> vel_cart, double angle);
// vel_cart: the initial velocity (direction) vector of the photon in Cartesian (x,y,z) coordinates. The position here is (0,0,0)
// angle: the half-opening angle of the spherical-toroidal geometry (radians)
// returns: True if the photon will immediately escape the torus, False if we need to let the MC code run. This will be used as escape condition #1 (see Escape Condition Notes)

bool escape_2(boost::qvm::vec<double, 3> pos_cart, boost::qvm::vec<double, 3> vel_cart, double R, double r, double half_oa, double ds);
// pos_cart: the current position vector of the photon in Cartesian (x,y,z) coordinates
// vel_cart: the current velocity (direction) vector of the photon in Cartesian (x,y,z) coordinates
// R: the outer radius of the spherical-toroidal geometry; the outer sphere's full radius
// r: the inner radius of the spherical-toroidal geometry
// half_oa: the half-opening angle of the torus (radians)
// ds: step length being used in the rest of the code
// returns: True if the photon will leave the torus without hitting the medium, False if it will intersect the medium again. This will be used as escape condition #2 (see Escape Condition Notes)

bool escape_3(boost::qvm::vec<double, 3> pos_cart, double R);
// pos_cart: the current position vector of the photon in Cartesian (x,y,z) coordinates
// R: the outer radius of the spherical-toroidal geometry; the outer sphere's full radius
// returns: True if photon has left the large sphere, False if photon is still inside it. This will be used as escape condition #3 (see Escape Condition Notes)

bool check_medium(boost::qvm::vec<double, 3> pos_cart, double r, double f);
// pos_cart: the current position vector of the photon in Cartesian (x,y,z) coordinates
// r: the inner radius of the spherical-toroidal geometry
// f: a quantity encapsulating the geometry (R, half_oa, etc). Check geometry notes for a derivation
// returns: True if photon is in the medium, False if photon is out of the medium. This will be used to determine whether scattering or absorption occur

bool scatter(double rand, double dtau);
// rand: a randomly generated double between 0 and 1
// dtau: the differential optical depth for the step length ds
// returns: True if a scattering happens, False if no scattering happens

boost::qvm::vec<double, 3> rotate(boost::qvm::vec<double, 3> vel, double theta, double phi);
// vel: the pre-scattering direction vector (unitless) of the photon
// theta: the Compton scattering angle
// phi: the azimuthal scattering angle
// returns: the post-scattering direction vector (unitless) of the photon

std::vector<double> propagate_photon(boost::multi_array<double, 2> table, double PL_ind, double R, double r, double half_oa, double f, double ds, double n);
// random: random numbers using the Mersenne Twister through Boost libraries
// table: the KN angle lookup table, implemented as a Boost multi array
// PL_ind: the instrinsic AGN power-law index (in terms of energy)
// R: the outer radius of the spherical-toroidal geometry; the outer sphere's full radius
// r: the inner radius of the spherical-toroidal geometry
// half_oa: the half-opening angle of the torus (radians)
// f: a quantity encapsulating the geometry (R, half_oa, etc). Check geometry notes for a derivation
// ds: step length
// n: the number density of scatterers (electrons) in the medium
// returns: a vector with the following information about the photon after its exit from the spherical torus - x, y, z, vx, vy, vz, N_scat, E_initial, E_final, weight

#endif /* full_test_functions_hpp */
