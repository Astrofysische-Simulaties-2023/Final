#ifndef integrators
#define integrators

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Bodyclass.h"
#include "Vec.h"

using namespace std;


int update_pos_vel_acc_EU(vector<Body>& bodies, double dt) {

    update_acc(bodies);

    for (int i = 0; i < bodies.size(); i++) {
        // Update velocities based on accelerations
        bodies[i].v += bodies[i].a * dt;
        // Update positions based on new velocities
        bodies[i].r += bodies[i].v * dt;
        // Reset the acceleration for this body to zero
        bodies[i].a = Vec(0.0, 0.0, 0.0);
        bodies[i].t += dt;
    }
    return 1;
}


int update_pos_vel_acc_RK4(vector<Body>& bodies, double dt) {
    vector<Body> k1 = bodies; // Store the initial state
    vector<Body> k2 = bodies;
    vector<Body> k3 = bodies;
    vector<Body> k4 = bodies;

    // Step 1 of RK4: Calculate derivatives at the beginning
    update_acc(k1);

    // Step 2 of RK4: Calculate k2 using half-step
    for (int i = 0; i < bodies.size(); i++) {
        k2[i].v = k1[i].v + 0.5 * k1[i].a * dt;
        k2[i].r = k1[i].r + 0.5 * k1[i].v * dt;
    }
    update_acc(k2);

    // Step 3 of RK4: Calculate k3 using another half-step
    for (int i = 0; i < bodies.size(); i++) {
        k3[i].v = k1[i].v + 0.5 * k2[i].a * dt;
        k3[i].r = k1[i].r + 0.5 * k2[i].v * dt;
    }
    update_acc(k3);

    // Step 4 of RK4: Calculate k4 using a full step
    for (int i = 0; i < bodies.size(); i++) {
        k4[i].v = k1[i].v + k3[i].a * dt;
        k4[i].r = k1[i].r + k3[i].v * dt;
    }
    update_acc(k4);
    
    // Update positions and velocities using the weighted average of k1, k2, k3, and k4
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].v = k1[i].v + (k1[i].a + 2 * k2[i].a + 2 * k3[i].a + k4[i].a) * dt / 6;
        bodies[i].r = k1[i].r + (k1[i].v + 2 * k2[i].v + 2 * k3[i].v + k4[i].v) * dt / 6;
        bodies[i].a = Vec(0.0, 0.0, 0.0);
        bodies[i].t += dt;
    }
    return 4;
};


int update_pos_vel_acc_VV(vector<Body>& bodies, double dt) {  

    // Update positions using full step and velocity using first half step
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].v += bodies[i].a * dt/2;
        bodies[i].r += bodies[i].v * dt ;
        bodies[i].a = Vec(0.0, 0.0, 0.0);
    }

    update_acc(bodies);

    // Update velocities using second half-step
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].v += bodies[i].a * dt/2;
        bodies[i].t += dt;
    }
    return 1;
}


int update_pos_vel_acc_FR(vector<Body>& bodies, double dt) {
    const double c = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));

    // Update positions using first half-step
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].r += 0.5 * bodies[i].v * c * dt;
    }
    update_acc(bodies);

    // Update velocities using first full step
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].v += bodies[i].a * c * dt;
        bodies[i].r += 0.5 * bodies[i].v * (1.0 - c) * dt;
        bodies[i].a = Vec(0.0, 0.0, 0.0);
    }
    
    update_acc(bodies);

    // Update velocities using second full step
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].v += bodies[i].a * (1.0 - 2.0 * c) * dt;
        bodies[i].r += 0.5 * bodies[i].v * (1.0 - c) * dt;
        bodies[i].a = Vec(0.0, 0.0, 0.0);
    }

    update_acc(bodies);

    // Update velocities using third full step
    for (int i = 0; i < bodies.size(); i++) {
        bodies[i].v += bodies[i].a * c * dt;
        bodies[i].r += 0.5 * bodies[i].v * c * dt;
        bodies[i].t += dt;
        bodies[i].a = Vec(0.0, 0.0, 0.0);
    }
    return 3;
};

#endif
