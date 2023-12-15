#ifndef classBody
#define classBody

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Vec.h"


using namespace std;

// double G = 6.67e-11; // m^3 s^-2 kg^-1
//double G = 2.6644e-3; // AE/M_earth * (km/s)^2, AE = 149.597.870.700 m, M_earth = 5.972e24 kg
//double G = 3.98589e5; // km/M_earth * (km/s)^2
double G = 8.887e-10; // AE/M_earth * (AE/day)^2, day = 86400 s

class Body {
public:
    double m,t,Epot,Ekin;
    Vec r,v,a;

    Body(double mass, Vec pos, Vec vel)
        : t(0),m(mass), r(pos), v(vel), a(0.0,0.0,0.0),Epot(0),Ekin(0) {
    }

    // Function to calculate and update accelerations between two bodies
    void calculate_acc(Body& body_2) {
        Vec dr = body_2.r-r;
        double R = dr.norm();
        
        // Avoid division by zero (softening parameter)
        if (R > 1e-10) {
            double acc_body_1 = G*body_2.m/(R*R);
            double acc_body_2 = G*m/(R*R);            
            
            // Update accelerations
            a += acc_body_1 * dr / R;
            body_2.a -= acc_body_2 * dr / R;
        }
    }

    void calculate_Epot(Body& body_2) {
        Vec dr = body_2.r-r;
        double R = dr.norm();

        // Avoid division by zero (softening parameter)
        if (R > 1e-10) {
            Epot -= G*m*body_2.m/R;
        }
    }
};

void update_acc(vector<Body>& bodies) {
    for (int i = 0; i < bodies.size(); i++) {
        for (int j = i + 1; j < bodies.size(); j++) {
            bodies[i].calculate_acc(bodies[j]);
        }
    }
};


double update_E(vector<Body>& bodies) {
    double Etot = 0.0;
    for (int i = 0; i < bodies.size(); i++) {
        for (int j = i + 1; j < bodies.size(); j++) {
            bodies[i].calculate_Epot(bodies[j]);
        }
    }
    for (int i = 0; i < bodies.size(); i++) {
        Etot += bodies[i].Epot;
        bodies[i].Epot = 0.0;

        bodies[i].Ekin = bodies[i].m * bodies[i].v.norm2()/2; //mv^2/2 = E_kin
        Etot += bodies[i].Ekin;
        bodies[i].Ekin = 0.0;
    }
    return Etot;
}


double calcR(vector<Body>& bodies, int numspacecrafts, int g) {
    double closestDistance = 100;

    for (int spacecraft = 0; spacecraft < numspacecrafts; ++spacecraft) {
        // Calculate distance from the current spacecraft to body g
        Vec dr = bodies[g].r - bodies[spacecraft].r;
        double distance = dr.norm();
        closestDistance = min(closestDistance, distance);
    }
    return closestDistance;
}



double update_dt(vector<Body>& bodies, vector<Body>& ref, double timestep, double min_dt, double power){
    double dt = timestep;
    
    vector<Body> helpref = bodies;
    update_acc(helpref);

    for (int i = 0; i < bodies.size(); i++) {
        double n = timestep * pow(ref[i].a.norm() / helpref[i].a.norm(), power); // acc (bodies at time =0) / acc bodies (at time = t)
        if (n<dt){
            dt = max(n,min_dt);
        }
    }
    
    return dt;
}


double update_dt2(vector<Body>& bodies, double timestep, double min_dt, double E, double Emax){

    double dt = 2*timestep;
    double Estep = 1;

    vector<Body> ref = bodies;
    update_acc(ref);

    while (Estep > Emax && dt >= min_dt) {
        vector<Body> helpref = ref;

        dt = dt*0.5;

        for (int i = 0; i < ref.size(); i++) {
            helpref[i].v += ref[i].a * dt;
            helpref[i].r += helpref[i].v * dt;
        }
        double E1 = update_E(helpref);
        Estep = abs(E1 - E);
    }
    return dt;
}



#endif
