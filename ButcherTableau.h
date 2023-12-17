#ifndef butchertableau
#define butchertableau


//Butcher tableau for a general integrator
//Author: Melvin de Bosscher
//Date: 29/11/2023

#define _USE_MATH_DEFINES
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "Bodyclass.h"
#include "Vec.h"
#include "Integrators.h"
#include <numeric>
using namespace std;

class ButcherTableau {
private:
    int s;
    vector<vector<double>> a;
    vector<double> b, c;

    // Utility function to print a vector
    void printVector(const vector<double>& vec) const {
        for (const auto& element : vec) {
            cout << element << " ";
        }
        cout << endl;
    }

    // Utility function to print a matrix
    void printMatrix(const vector<vector<double>>& mat) const {
        for (const auto& row : mat) {
            printVector(row);
        }
    }

public:
    // Constructor
    ButcherTableau(vector<vector<double>> coefficients, vector<double> weights, vector<double> nodes, const int stages = -1)
        : a(coefficients), b(weights), c(nodes), s(stages) {
        // if s isn't specified
        if (s == -1) {
            s = nodes.size();
        }
    }
    
    // Public member functions to access a, b, c and s
    const vector<vector<double>>& getA() const {
        return a;
    }

    const vector<double>& getB() const {
        return b;
    }

    const vector<double>& getC() const {
        return c;
    }

    const int& getS() const {
        return s;
    }
};

// UPDATE FUNCTIONS
// updates an explicit integrator using butcher tableaus
int update_integrator(ButcherTableau integrator, vector<Body>& bodies, double h) {
    vector<vector<double>> A = integrator.getA();  // The coefficients are part of a matrix A
    vector<double> B = integrator.getB(); // The weight is a vector B
    int S = integrator.getS(); // The amount of stages is an integer S
    vector<vector<Body>> k(S); // k is a vector of a vector containing the information of multiple bodies

    int no_driver_functions = 0;
    for (int m = 0; m < S; m++) {   // for S stages
        k[m] = bodies; // store initial state
        for (int i = 0; i < bodies.size(); i++) { // for n bodies
            k[m][i].v = k[0][i].v;
            k[m][i].r = k[0][i].r;
            for (int n = 0; n < m; n++) {
                k[m][i].v += h * A[m][n] * k[n][i].a; // update velocity of k_m
                k[m][i].r += h * A[m][n] * k[n][i].v; // update position of k_m
            }
        }
        update_acc(k[m]); // update acceleration of k_m
        no_driver_functions += 1;
    }

    // Update positions and velocities using the weighted average of k's
    for (int i = 0; i < bodies.size(); i++) {   // for n bodies
        bodies[i].v = k[0][i].v;
        bodies[i].r = k[0][i].r;
        for (int m = 0; m < S; m++) {   // for S stages
            bodies[i].v += h * B[m] * k[m][i].a; // update velocity of body
            bodies[i].r += h * B[m] * k[m][i].v; // update position of body
        }
        bodies[i].a = Vec(0.0, 0.0, 0.0); // update acceleration of body
        bodies[i].t += h;   // update time of body
    }
    return no_driver_functions;
};


// returns the k's in order to give an initial guess to the implicit integrator
vector<vector<Body>> initial_guess(ButcherTableau integrator, vector<Body>& bodies, double h) {

    vector<vector<double>> A = integrator.getA(); // Coefficients of the initial guess
    int s = integrator.getS(); // The amount of stages is an integer 's'

    // k is a vector of a vector containing the information of multiple bodies
    vector<vector<Body>> k(s); // k is the initial guess
    
    for (int m = 0; m < s; m++) {   // for s stages
        k[m] = bodies; // store initial state
        for (int i = 0; i < bodies.size(); i++) { // for n bodies
            k[m][i].v = k[0][i].v;
            k[m][i].r = k[0][i].r;
            for (int n = 0; n < m; n++) { // up to m-1
                k[m][i].v += h * A[m][n] * k[n][i].a; // update velocity of k_m
                k[m][i].r += h * A[m][n] * k[n][i].v; // update position of k_m
            }
        }
        update_acc(k[m]); // update acceleration of k_m
    }
    return k;
}

// updates an implicit integrator using the butcher tableau using an initial guess of choice
int update_implicit_integrator(ButcherTableau integrator, ButcherTableau guess, vector<Body>& bodies, double h) {
    vector<vector<double>> A = integrator.getA();  // The coefficients are part of a matrix 'a'
    vector<double> B = integrator.getB(); // The weight is a vector 'b'
    int s = integrator.getS(); // The amount of stages is an integer 's'
    vector<vector<Body>> l(s); // l is a vector of a vector containing the different bodies

    vector<vector<Body>> k = initial_guess(guess, bodies, h);
    int no_driver_functions = 0;

    // implicit method
    for (int m = 0; m < s; m++) {   // for s stages
        l[m] = bodies; // store initial state
        for (int i = 0; i < bodies.size(); i++) { // for n bodies
            l[m][i].v = k[0][i].v;
            l[m][i].r = k[0][i].r;
            for (int n = 0; n < s; n++) { // up to s-1
                l[m][i].v += h * A[m][n] * k[n][i].a; // update velocity of l_m
                l[m][i].r += h * A[m][n] * k[n][i].v; // update position of l_m
            }
        }
        update_acc(l[m]); // update acceleration of l_m
        no_driver_functions += 2;
    }

    // Update positions and velocities using the weighted average of k's
    for (int i = 0; i < bodies.size(); i++) {   // for n bodies
        bodies[i].v = k[0][i].v;
        bodies[i].r = k[0][i].r;
        for (int m = 0; m < s; m++) {   // for s stages
            bodies[i].v += h * B[m] * l[m][i].a; // update velocity of body
            bodies[i].r += h * B[m] * l[m][i].v; // update position of body
        }
        bodies[i].a = Vec(0.0, 0.0, 0.0); // update acceleration of body
        bodies[i].t += h;   // update time of body
    }
    return no_driver_functions;
}

// EXPLICIT METHODS
// Butcher Tableau for RK4 method
ButcherTableau RK4() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0, 0.0}, 
                                            {0.5, 0.0, 0.0, 0.0}, 
                                            {0.0, 0.5, 0.0, 0.0}, 
                                            {0.0, 0.0, 1.0, 0.0}};
    vector<double> weights = {1.0/6, 1.0/3, 1.0/3, 1.0/6};  // Weights
    vector<double> nodes = {0.0, 0.5, 0.5, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Explicit midpoint method
ButcherTableau Explicit_midpoint() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0},
                                            {0.5, 0.0}};
    vector<double> weights = {0.0, 1.0};  // Weights
    vector<double> nodes = {0.0, 0.5};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Heun's method
ButcherTableau Heun() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0},
                                            {1.0, 0.0}};
    vector<double> weights = {0.5, 0.5};  // Weights
    vector<double> nodes = {0.0, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Ralston's method
ButcherTableau Ralston() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0},
                                            {2.0/3, 0.0}};
    vector<double> weights = {1.0/4, 3.0/4};  // Weights
    vector<double> nodes = {0.0, 2.0/3};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for generic second-order method
ButcherTableau second_order(double alpha) {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0},
                                            {alpha, 0.0}};
    vector<double> weights = {1.0 - (1.0/(2.0 * alpha)), 1.0/(2.0 * alpha)};  // Weights
    vector<double> nodes = {0.0, alpha};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for RK3 method
ButcherTableau RK3() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0},
                                            {0.5, 0.0, 0.0},
                                            {-1.0, 2.0, 0.0}};
    vector<double> weights = {1.0/6, 2.0/3, 1.0/6};  // Weights
    vector<double> nodes = {0.0, 0.5, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Heun's thirs-order method
ButcherTableau Heun3() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0},
                                            {1.0/3, 0.0, 0.0},
                                            {0.0, 2.0/3, 0.0}};
    vector<double> weights = {1.0/4, 0.0, 3.0/4};  // Weights
    vector<double> nodes = {0.0, 1.0/3, 2.0/3};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Van der Houwen's/Wray third-order method
ButcherTableau VDH() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0},
                                            {8.0/15, 0.0, 0.0},
                                            {1.0/4, 5.0/12, 0.0}};
    vector<double> weights = {1.0/4, 0.0, 3.0/4};  // Weights
    vector<double> nodes = {0.0, 8.0/15, 2.0/3};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Ralston's third-order method
ButcherTableau Ralston3() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0},
                                            {0.5, 0.0, 0.0},
                                            {0.0, 3.0/4, 0.0}};
    vector<double> weights = {2.0/9, 1.0/3, 4.0/9};  // Weights
    vector<double> nodes = {0.0, 0.5, 3.0/4};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Third-order Strong Stability Preserving Runge-Kutta method
ButcherTableau SSPRK3() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0},
                                            {1.0, 0.0, 0.0},
                                            {1.0/4, 1.0/4, 0.0}};
    vector<double> weights = {1.0/6, 1.0/6, 2.0/3};  // Weights
    vector<double> nodes = {0.0, 1.0, 0.5};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for a generic third-order method
ButcherTableau third_order(double x) {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0},
                                            {x, 0.0, 0.0},
                                            {1.0 + (1.0 - x) / (x * (3.0 * x - 2.0)), -(1.0 - x) / (x * (3.0 * x - 2.0)), 0.0}};
    vector<double> weights = {0.5 - 1.0/(6.0*x), 1.0/((6.0*x)*(1.0-x)), (2.0-3.0*x)/((6.0)*(1.0-x))};  // Weights
    vector<double> nodes = {0.0, x, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for 3/8-method
ButcherTableau Three_Eighth() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0, 0.0}, 
                                            {1.0/3, 0.0, 0.0, 0.0}, 
                                            {-1.0/3, 1.0, 0.0, 0.0}, 
                                            {1.0, -1.0, 1.0, 0.0}};
    vector<double> weights = {1.0/8, 3.0/8, 3.0/8, 1.0/8};  // Weights
    vector<double> nodes = {0.0, 1.0/3, 2.0/3, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Ralston's fourth-order method
ButcherTableau Ralston4() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0, 0.0}, 
                                            {0.4, 0.0, 0.0, 0.0}, 
                                            {.29697761, .15875964, 0.0, 0.0}, 
                                            {.21810040, -3.05096516, 3.83286476, 0.0}};
    vector<double> weights = {.17476028, -.55148066, 1.20553560, .17118478};  // Weights
    vector<double> nodes = {0.0, 0.4, .45573725, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for RK5
ButcherTableau RK5(){
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                        { 1.0/4, 0.0, 0.0, 0.0, 0.0, 0.0},
                                        {1.0/8, 1.0/8, 0.0, 0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.5, 0.0, 0.0, 0.0},
                                        {3.0/16, -3.0/8, 3.0/8, 9.0/16, 0.0, 0.0},
                                        {-3.0/7, 8.0/7, 6.0/7, -12.0/7, 8.0/7, 0.0}};                   
    vector<double> weights {7.0/90, 0.0, 32.0/90, 12.0/90, 32.0/90, 7.0/90};
    vector<double> nodes = {0.0, 1.0/4, 1.0/4, 0.5, 3.0/4, 1.0};
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
}

// Butcher Tableau for RK6
ButcherTableau RK6(){
    vector<vector<double>> coefficients = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                        { 1.0/3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                                        {0.0, 2.0/3, 0.0, 0.0, 0.0, 0.0, 0.0},
                                        {1.0/12, 1.0/3, -1.0/12, 0.0, 0.0, 0.0, 0.0},
                                        {-1.0/16, 9.0/8, -3.0/16, -3.0/8, 0.0, 0.0, 0.0},
                                        {0.0, 9.0/8, -3.0/8, -3.0/4, 0.5, 0.0, 0.0},
                                        {9.0/44, -9.0/11, -63.0/44, 18.0/11, 0.0, -16.0/11, 0.0}};                   
    vector<double> weights {11.0/120, 0.0, 27.0/40, 27.0/40, -4.0/15, -4.0/15, 11.0/120};
    vector<double> nodes = {0.0, 1.0/3, 2.0/3, 1.0/3, 0.5, 0.5, 1.0};
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
}

// IMPLICIT METHODS
// Butcher Tableau for Symmetric and Symplectic tree-stage
ButcherTableau SRK3() {
        // Coefficients
    vector<vector<double>> coefficients = {{5.0/36, 2.0/9, 5.0/36-sqrt(15.0)/10}, 
                                            {5.0/36, 2.0/9, 5.0/36}, 
                                            {5.0/36+sqrt(15.0)/10, 2.0/9, 5.0/36}};
    vector<double> weights = {5.0/18, 4.0/9, 5.0/18};  // Weights
    vector<double> nodes = {0.5-sqrt(15.0)/10, 0.5, 0.5+sqrt(15.0)/10};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for Symmetric and Symplectic two-stage
ButcherTableau SRK2() {
        // Coefficients
    vector<vector<double>> coefficients = {{1.0/4, 1.0/4-sqrt(3.0)/6},  
                                            {1.0/4+sqrt(3.0)/6, 1.0/4}};
    vector<double> weights = {0.5, 0.5};  // Weights
    vector<double> nodes = {1.0/2-sqrt(3.0)/6, 1.0/2+sqrt(3.0)/6};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for the third-order Radau IA method
ButcherTableau RIA3() {
    // Coefficients
    vector<vector<double>> coefficients = {{1.0/4, -1.0/4},  
                                            {1.0/4,5.0/12}};
    vector<double> weights = {1.0/4, 3.0/4};  // Weights
    vector<double> nodes = {0.0, 2.0/3};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};

// Butcher Tableau for the Crank-Nicolson method
ButcherTableau Crank_Nicolson() {
    // Coefficients
    vector<vector<double>> coefficients = {{0.0, 0.0},  
                                            {0.5, 0.5}};
    vector<double> weights = {0.5, 0.5};  // Weights
    vector<double> nodes = {0.0, 1.0};    // Nodes
    ButcherTableau tableau(coefficients, weights, nodes);
    return tableau;
};
#endif
