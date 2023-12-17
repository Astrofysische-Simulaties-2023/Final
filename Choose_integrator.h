#ifndef choose_int
#define choose_int

//Functions used to choose the integrator, change filename etc.
//Author: Robrecht Keijzer
//Date: 24/11/2023

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include "Integrators.h"
#include "ButcherTableau.h"
using namespace std;



int choose_integrator() {
    std::string integr_num = "Bad input";    //By first storing the integrator as a string, we are better protected from bad input types
    while (integr_num == "Bad input") {
        cout << "==============================\n";
        cout << "1) Euler\n2) Fourth Order Runge Kutta\n3) Velocity Verlet\n4) Forest Ruth\n5) Implicit RK3\n6) Fifth Order Runge Kutta\n";
        cout << "==============================\n";
        cout << "Choose an integrator (number from 1-6): ";
        cin >> integr_num;
        
        if (integr_num != "1" && integr_num != "2" && integr_num != "3" && integr_num != "4" && integr_num != "5" && integr_num != "6") {
            integr_num = "Bad input";
            cout << "\nThis is not a valid option, try again.\n";
        }
    }
    return stoi(integr_num);    //Convert string to integer
}




//Using "typedef" to rename the integrator function to the generic integrator_func.
typedef int (*IntegratorFunction)(vector<Body>& bodies, double dt); 

//The integrator function should return two outputs: the name (to name the text file), and the integrator function.
//They can most easily be combined in a struct.
struct Integrator_struct {
    std::string name;
    IntegratorFunction integrator_func;
};


string make_filename(int numberbodies, string integr_name, double dt) {     //filename shows number of bodies, integrator name and timestep (e.g. 1e-2)
    std::stringstream timestep_stream;
    timestep_stream << std::fixed << std::setprecision(1) << log(dt)/log(10);
    std::string timestep_string = timestep_stream.str();
    string filename = to_string(numberbodies) + "body_" + integr_name + "_dt_1e" + timestep_string + ".txt";
    return filename;
}



//Helper function for the SRK3
int SRK3_function(vector<Body>& bodies, double dt) {
    return update_implicit_integrator(SRK3(), RK3(), bodies, dt);
}

//Helper function for the RK6
int RK5_function(vector<Body>& bodies, double dt) {
    return update_integrator(RK5(), bodies, dt);
}


Integrator_struct integrator(int integr_num) {
    Integrator_struct integr;
    if (integr_num == 1){
        integr.name = "Euler";
        integr.integrator_func = update_pos_vel_acc_EU;
    } else if (integr_num == 2){
        integr.name = "RK4";
        integr.integrator_func = update_pos_vel_acc_RK4;
    } else if (integr_num == 3){
        integr.name = "Verlet";
        integr.integrator_func = update_pos_vel_acc_VV;
    } else if (integr_num == 4){
        integr.name = "FoRu";
        integr.integrator_func = update_pos_vel_acc_FR;
    } else if (integr_num == 5){
        integr.name = "implicit_RK3";
        integr.integrator_func = SRK3_function;
    } else if (integr_num == 6) {
        integr.name = "RK5";
        integr.integrator_func = RK5_function;
    }
    return integr;
}



#endif