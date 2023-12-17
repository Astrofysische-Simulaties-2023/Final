#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <chrono> // Include the chrono library
#include <random>
#include <string>

#include "Bodyclass.h"
#include "Integrators.h"
#include "Vec.h"
#include "ButcherTableau.h"
#include "Choose_integrator.h"

using namespace std;

int main() {

    //===== Declare the bodies in our system =====
    vector<Body> bodies; //1977-Sep-06 00:00:00.0000 TDB
    bodies.emplace_back(332959.936, Vec(2.374359334804984e-3, -4.397509852228428e-3, -6.403727192979392e-5), Vec(8.067422670429710e-6, 3.672787845010753e-8, -2.197171878619823e-7)); // sun 
    bodies.emplace_back(0.815195, Vec(1.274926136407149E-01, 7.048882120716604E-01,2.362234183952712E-03), Vec(-1.998001176601727E-02, 3.412330568007654E-03, 1.200315569581762E-03)); // venus
    bodies.emplace_back(0.01230536872, Vec(9.699068156158732e-1,-2.862430064941408e-1, -2.871929438588917e-4), Vec(4.043950282610457e-3, 1.655300393597739e-2, -2.280842502739660e-5)); //moon 
    bodies.emplace_back(1, Vec(9.693415425602092e-1, -2.888761757432423e-1, -7.064680860064746e-5), Vec(4.590127153526758e-3, 1.643984412273940e-2, 8.985875982417776e-7)); //earth 
    bodies.emplace_back(0.107449696, Vec(8.823367576696816e-1, 1.186540676671558, 3.201641663305391e-3), Vec(-1.071179986117151e-2, 9.504915586236624e-3, 4.628247779093462e-4)); // mars
    bodies.emplace_back(317.8377145, Vec(7.023069539386351e-1, 5.044885561956213, -3.650778429292263e-2), Vec(-7.563035000680561e-3, 1.388374806513800e-3, 1.636331181022656e-4)); // jupiter
    bodies.emplace_back(95.16442042, Vec(-7.191297133554489, 5.703064303254036, 1.860207681470850e-1), Vec(-3.766914204757593e-3, -4.385970133194619e-3, 2.266657950446760e-4)); // saturn
    bodies.emplace_back(14.5362, Vec(-13.88051911770789, -12.36710557923142, 1.341795439021756e-1), Vec(2.587073048286647e-3, -3.120394186515559e-3, -4.508153960231010e-5)); // uranus
    // spacecrafts 1977-Sep-06 00:00:00.0000 TDB
    bodies.emplace_back(1.208769e-22, Vec(9.703073590971770e-1, -2.864904318339970e-1, 1.066427909013482e-4), Vec(6.691403977335315e-3, 2.202627662224282e-2, 4.216890028880711e-4)); // Voyager 1
    bodies.emplace_back(1.208769e-22, Vec(9.909656753661980e-1, -2.017144553381927e-1, 3.033912016204687e-2), Vec(5.834315231766972e-3, 2.163241124839557e-2, 1.807574968257737e-3)); // Voyager 2


    //2023-Jan-01 
    //bodies.emplace_back(333000, Vec(-9.055145044460892e-3, 9.495765501773467e-5, 2.101594374427968e-4), Vec(9.450438941100025e-7,-8.999856517644221e-6, 5.175772303295750e-8)); // sun 
    //bodies.emplace_back(0.815195, Vec(5.521118597473940e-1, -4.622937696473001e-1, -3.851853938626681e-2), Vec(1.273466166130406e-2, 1.552337997072891e-2, -5.214742484350385e-4)); // venus
    //bodies.emplace_back(1, Vec(-1.793111147418941e-1,9.685804332587541e-1, 1.612991618128246e-4), Vec(-1.721950145453741e-2, -3.058649113379992e-3, -3.172877676994421e-7)); //earth 
    //bodies.emplace_back(0.0123, Vec(-1.771353884440425e-1,9.700106846238624e-1, 1.283977887164111e-4), Vec(-1.751100643076462e-2, -2.558486718724795e-3, 5.036175128214648e-5)); //moon 
    //bodies.emplace_back(0.10745, Vec(5.1823e-2, 1.56285, 3.147e-2), Vec(-1.3452e-2, 1.72456e-3, 3.664e-4)); // mars
    //bodies.emplace_back(317.8, Vec(4.829169403525787, 1.045331598448545, -1.123782328236326e-1), Vec(-1.682521975834704e-3, 7.729815369519151e-3, 5.582305732884232e-6)); // jupiter 
    //bodies.emplace_back(95.16576, Vec(8.138331649624233, -5.505029974575755, -2.283049437130539e-1), Vec(2.814258148133007e-3, 4.610041364679036e-3, -1.921624943094593e-4)); // saturn
    //bodies.emplace_back(0.02253, Vec(8.139815554731646, -5.497947804173316, -2.321048818296638e-1), Vec(-3.177500951150720e-4, 5.338783850066382e-3, -2.558360732475707e-4)); // Titan, moon saturn
    //bodies.emplace_back(14.5362, Vec(13.36068864772597, 14.42715963368727, -1.195070374691794e-1), Vec(-2.914435855150269e-3, 2.489162061236050e-3, 4.701464111661541e-5)); // Uranus
    //bodies.emplace_back(17.147646, Vec(29.75180913013332, -2.942716536143852, -6.250608682820173e-1), Vec(2.884027527271555e-4, 3.142669796707304e-3, -7.137268895852339e-5)); // Neptunus
    

    // 2038-Dec-08
    //bodies.emplace_back(333000, Vec(0.005637051637579521858, -0.0051392429845828228527, -0.00013812429071061654373), Vec(5.2186034737313050552e-06, 6.7145126859188918367e-06, -1.6448510404081506365e-07)); // sun
    //bodies.emplace_back(0.815195, Vec(0.21723721852153748957, -0.7011156750412541605, -0.02193017346876506915), Vec(0.019221076400216426106, 0.0058215919607687936083, -0.0010286159333626144954)); // venus
    //bodies.emplace_back(1, Vec(0.23629895966615929725, 0.9524973128926194077, -0.00022030166036129269721), Vec(-0.016997330415421379318, 0.0039648862632960732433, -6.2305182597280914753e-07)); //earth 
    //bodies.emplace_back(0.0123, Vec(0.23875584194411275663, 0.95360983776370322751, -0.00044809616005550480678), Vec(-0.017217405680334394585, 0.0044788517138313766602, 1.9074438536576441384e-05)); //moon 
    //bodies.emplace_back(0.10745, Vec(-0.82282417660557516381, -1.2943034602250287346, -0.0068734333014412530086), Vec(0.012305301326269196815, -0.0063572513494782514204, -0.00043485451391617446206)); // mars
    //bodies.emplace_back(317.8, Vec(-4.1146114650062335372, 3.3932698469983391121, 0.077924897417690180967), Vec(-0.0048902244022227875017, -0.0054706548578263951288, 0.0001321252489386841734)); // jupiter
    //bodies.emplace_back(95.16576, Vec(-9.3194751466275516094, 1.3088309165180604943, 0.34823790463205417378), Vec(-0.0010785846694082925224, -0.0055369665438612086089, 0.00013916413729544085606)); // saturn
    //bodies.emplace_back(14.5362, Vec(-8.0830065618934519023, 16.798255477084886422, 0.16706743080264901957), Vec(-0.003573241720737723847, -0.0018908607204464562487, 3.9249280327545366315e-05)); // Uranus
    //bodies.emplace_back(17.147646, Vec(25.897769578438747828, 14.756312791692359099, -0.90077435372806147207), Vec(-0.0015734695422524277216, 0.0027465947338189879851, -2.0307548080665057967e-05)); // Neptunus
    // new spacecrafts
    //bodies.emplace_back(1.25e-22, Vec(0.229, 0.96, -0.00018264), Vec(-2.0e-2, 1.295765e-2, 4.9040e-4)); // Voyager 3
    //bodies.emplace_back(1.25e-22, Vec(0.229, 0.96, -0.00018264), Vec(-2.0e-2, 1.295760e-2, 4.9040e-4)); // Voyager 4
    //bodies.emplace_back(1.25e-22, Vec(0.229, 0.96, -0.00018264), Vec(-2.0e-2, 1.295755e-2, 4.9040e-4)); // Voyager 5


    //Burrau's problem
    //bodies.emplace_back(5000, Vec(0.0,0.0,0.0), Vec(0.0,0.0,0.0));  
    //bodies.emplace_back(4000, Vec(3.0,0.0,0.0), Vec(0.0,0.0,0.0));  
    //bodies.emplace_back(3000, Vec(0.0,4.0,0.0), Vec(0.0,0.0,0.0)); 
    

    int integr_num = choose_integrator();                       //Choose the integrator at user input
    Integrator_struct integr_struct = integrator(integr_num);   //Struct with name and integrator function
    IntegratorFunction integrator_func = integr_struct.integrator_func;
    int numberbodies = bodies.size();

    // Initial Values
    double timestep = 0.01; // Max allowed time step
    double dt = timestep; // time step (for adaptive time steps: dt != timestep)
    double max_t = 3000;
    double numsteps = max_t/dt;
    int stepCounter = 0;  
    int outf = 50; // outfile every ... steps 

    //Energy and reference to initial configuration
    vector<Body> ref = bodies;  // bodies at t=0
    update_acc(ref);
    int no_driver_functions = 0; // number of driver function evaluations
    double E_0 = update_E(bodies);
    double E = E_0;
    double dE = 0.0;

    //adaptive update_dt
    double power  = 0.5; //power law (update_dt)
    double min_dt = 0.005;  // Minimum allowed time step
    double Emax = 1e-12; // Max allowed energy error for each dt step (update_dt2)


    //When using velocity verlet, calculate the acceleration a first time.
    if (integr_num == 3) {  
        update_acc(bodies);
        no_driver_functions += 1;}

    //Write to text file
    string filename = make_filename(numberbodies, integr_struct.name, dt); 
    cout << "Writing to file '" << filename << "'" << endl;
    ofstream outfile(filename);
    outfile << setprecision(15);
    outfile << "Time" << ' ' << "Xpos" << ' ' << "Ypos" << ' ' << "Zpos" << ' ' << "E" << ' ' << "dE" << ' ' << "dt" << '\n';

    
    auto start_time = chrono::high_resolution_clock::now(); // Start timing 


    while (max_t - bodies[0].t > dt/2) {
        
        // Output the state of the bodies every ... steps.
        if (stepCounter % outf == 0) {
            for (int i = 0; i < bodies.size(); i++) {
                outfile << bodies[i].t << ' ' << bodies[i].r.x() << ' ' << bodies[i].r.y() << ' ' << bodies[i].r.z() << ' ' << E << ' ' << dE << ' ' << dt << '\n';
            }
        }

        // Increment the step counter
        stepCounter++;


        //===== Uncomment this if you want to use adaptive =====
        //dt =  update_dt(bodies, ref, timestep, min_dt, power, no_driver_functions); // Use this when using change in acceleration (power law of the normalized acceleration of the bodies)
        //dt = update_dt2(bodies, timestep, min_dt, E, Emax, no_driver_functions); // Use this when using energy accuracy (when energy error for a specific time step is too big, divide the time step by 2)

        // Update positions velocities and acceleration of all bodies
        no_driver_functions += integrator_func(bodies, dt);
        E = update_E(bodies);
        dE = (E-E_0)/E_0;
    }

    auto end_time = chrono::high_resolution_clock::now(); // End timing
    chrono::duration<double> elapsed_time = end_time - start_time;

    outfile << bodies.size() << ' ' << timestep << ' ' << elapsed_time.count() << ' ' << elapsed_time.count()/numsteps << ' '<< no_driver_functions << ' ' << no_driver_functions/elapsed_time.count() << ' '<< no_driver_functions/max_t << ' ' <<'\n';
    outfile.close();
    cout << "Done writing." << endl << endl;

    cout << "It took " << elapsed_time.count() << " seconds for " << bodies.size() << " bodies." << endl;
    cout << "It took " << elapsed_time.count()/numsteps << " seconds per integration step." << endl;

    cout << "Total number of driver function evaluations: " << no_driver_functions << endl;
    cout << "Function evalutions per second: " << no_driver_functions/elapsed_time.count() << endl;
};
