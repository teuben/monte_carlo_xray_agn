//
//  main.cpp
//  full_test
//
//  Created by Yash Gursahani.

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "/usr/local/include/omp.h" // included for multithreading in the photon loop
#include "full_test_functions.hpp"
#include "boost/qvm.hpp"
#include "boost/multi_array.hpp" // included so that we can read KN file into a table
#include "boost/random/mersenne_twister.hpp" // included so that we can generate random integers properly. Uses Matsumoto & Nishimura 1998 prescription
#include "boost/random/uniform_int_distribution.hpp" // with this, we can generate a random number in a predetermined range using the Mersenne Twister ^

using namespace std::chrono;

int main() {

    // ///////////////////////////////////////// EDIT INPUT PARAMETERS HERE! ////////////////////////////////////////////////
    
    double n = 1.82e6; // cm^-3, scatterer (electron) number density
    double PL_ind = 2.0; // intrinsic source power-law index
    int N_phot = 1000000; // the number of photons we want to shoot out
    int phot_notif = 100000; // the number of photons after which you want to receive a notification that the program is running
    int phot_store = 1000; // the number of photons after which you want to store the data in the output file
    
    double r = 0.5 * pc_cm; // cm, the inner radius of the spherical-toroidal geometry - Li & Shen 2023
    double R = 1.5 * r; // cm, the outer radius of the spherical-toroidal geometry (radius of the full "sphere") - Li & Shen 2023
    double half_oa = pi/6.0; // radians, the half-opening angle of the torus
    double ds = 7.71e14; // cm, step length
    
    // Declaring some variables to reduce check_medium runtime
    double phi = (pi/2.0) - half_oa;
    double D = sqrt(pow(R, 2.0) + pow(r, 2.0) - 2*r*R*cos(phi));
    double f = R*sin(phi)/D; 
    
    bool writeout = false; // a switch for writing the results of the run to a file for analysis
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Reading in data from theta lookup table
    size_t E_size = 2001;
    size_t f_size = 1001;
    
    typedef boost::multi_array<double, 2> array_type;
    typedef array_type::index index;
    array_type table(boost::extents[E_size][f_size]);
    
    size_t i = 0; // table row index (energies)
    size_t j = 0; // table column index (fractions)
    
    std::ifstream KN;
    KN.open("/Users/yashgursahani/Desktop/UMD_Research/high_z_agn/Compton Scattering Test/kn_table/kn_table/kn_table.txt");
    if (KN.is_open()) {
        for (std::string line; std::getline(KN, line, ' ');) { //
            if (j < f_size) { // moving fraction by fraction (column by column)
                table[i][j] = std::stod(line);
                j = j + 1;
            }
                
            else if (j == f_size && i < E_size - 1) { // when we want to move to the next energy (row)
                i = i + 1;
                table[i][0] = std::stod(line);
                j = 1;
            }
                
            else {break;}
                
        }
            
    }

    KN.close(); 
    
/////////////////////////////////////////////// Start of the MC code ///////////////////////////////////////////////
    
// All output arrays
    // double pos_X[phot_store];
    // double pos_Y[phot_store];
    // double pos_Z[phot_store];

    // double vel_X[phot_store];
    // double vel_Y[phot_store];
    // double vel_Z[phot_store];

    // double num_scatters[phot_store]; // how many scatterings each photon has been through
    // double phot_E_initial[phot_store]; // the initial energy of each photon
    // double phot_E_final[phot_store]; // final energy of each photon
//    double weight[N_phot]; // the weight of each photon, determined by its initial energy

// Read the output directory from an environment variable in the CMakeLists.txt file
    const char* output_dir = std::getenv("OUTPUT_DIR");
    if (!output_dir) {
        output_dir = ".";  // Default to current directory if the environment variable is not set
    }
    
    // write out to file
    std::ofstream output; // create output file for our table
    if (writeout) {
        output.open(std::string(output_dir) + "/full_test_run_log6_phot_O3.txt"); // open the file and name it
        output << "x y z vx vy vz N_scat E_initial E_final" << "\n"; // all column names at the top of the file
    }
    
    // checking runtimes - start clock
    auto start = high_resolution_clock::now();
    
    // Start of the loop

    int writeout_counter = 0; // counter for when to write out to file
    std::stringstream ss; // to accumulate data in memory

    for (int m = 0; m < N_phot; m++) {
            
        if (m % phot_notif == 0) {
            std::cout << "Still running.. " << m/N_phot * 100 << "%" << std::endl; // update every phot_notif photons
        }

        if (writeout_counter == phot_store) {
            writeout_counter = 0; // reset the counter

            if (writeout) { // write last phot_store photons to file
            //     for (int k = 0; k < phot_store; k++) {
            //         ss << pos_X[k] << " " << pos_Y[k] << " " << pos_Z[k] << " " 
            //         << vel_X[k] << " " << vel_Y[k] << " " << vel_Z[k] << " "
            //         << num_scatters[k] << " " << phot_E_initial[k] << " " 
            //         << phot_E_final[k] << "\n";
            //     }

                // Write the accumulated stringstream data to the output file at once
                output << ss.str();
                ss.str("");   // Clear the stringstream for the next chunk
                ss.clear();   // Reset any flags (like EOF or fail)
            }
        }

        else {
            writeout_counter = writeout_counter + 1; // increment the counter
        }
            
        std::vector<double> phot_info = propagate_photon(table, PL_ind, R, r, half_oa, f, ds, n);
            
        // pos_X[writeout_counter] = phot_info[0];
        // pos_Y[writeout_counter] = phot_info[1];
        // pos_Z[writeout_counter] = phot_info[2];

        // vel_X[writeout_counter] = phot_info[3];
        // vel_Y[writeout_counter] = phot_info[4];
        // vel_Z[writeout_counter] = phot_info[5];

        // num_scatters[writeout_counter] = phot_info[6];
        // phot_E_initial[writeout_counter] = phot_info[7];
        // phot_E_final[writeout_counter] = phot_info[8];

        double pos_X = phot_info[0];
        double pos_Y = phot_info[1];
        double pos_Z = phot_info[2];

        double vel_X = phot_info[3];
        double vel_Y = phot_info[4];
        double vel_Z = phot_info[5];

        double num_scatters = phot_info[6];
        double phot_E_initial = phot_info[7];
        double phot_E_final = phot_info[8];

        if (writeout) { // store photon info in memory
            ss << pos_X << " " << pos_Y << " " << pos_Z << " " 
            << vel_X << " " << vel_Y << " " << vel_Z << " "
            << num_scatters << " " << phot_E_initial << " " 
            << phot_E_final << "\n";
        }


        
    }
    
    // checking runtimes - stop the clock and calculate runtime
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop-start);
    std::cout << "Time for this run: " << duration.count() << std::endl;
    
    return 0;
}
