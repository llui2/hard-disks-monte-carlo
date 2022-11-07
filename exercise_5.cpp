#include <iostream> // terminal input-output
#include <random> // random numbers
#include <vector> // vectors
#include <fstream> // files input-output
#include <string> // strings
#include <ctime> // calculate cpu times
#include <time.h> // print current hour
#include <algorithm> // min

int sign(float x){
    /* This function is equivalent to the sing(x) math function,
    returns +1 if x is positive and -1 if x is negative. */
    if (signbit(x)) {
        return 1;
    }
    else {
        return -1;
    }
}

double metric(double x, double x_0, double y, double y_0) {
    /* This function calculates distance between two points 
    in the plane (x,y) (x_0,y_0) */

    return sqrt(pow((x-x_0),2) + pow((y-y_0),2));
}

void generate_config(int N, double L_x, double L_y, double sigma,std::vector< std::vector<double> > &p_config) {
    /* This function generates a initial configuration of N solid discs of radius sigma
    in a LxL box */

    std::random_device seed; // generate a random seed (random integer)
    std::mt19937 generator(seed()); // initialise mt19937 random number generator with the seed
    std::uniform_real_distribution<float> uniform_one(-1, 1); // function to generate a uniform real random number between -1 and 1

    for (int i=0; i<N; i++) {
overlap:
        p_config[i][0] = (L_x/2-sigma/2)*uniform_one(generator);
        p_config[i][1] = (L_y/2-sigma/2)*uniform_one(generator);
        for (int j=0; j < i; j++) {
            if (metric(p_config[i][0],p_config[j][0],p_config[i][1],p_config[j][1]) < sigma) {
                goto overlap;
            }
        }
    }
}

bool metropolis(double N, double L_x, double L_y, double T, double sigma, double g, double delta,
    std::vector< std::vector<double> > &p_config_old, std::vector< std::vector<double> > &p_config_new) {
    /* This function does a Monte-Carlo step */
    
    std::random_device seed; // generate a random seed (random integer)
    std::mt19937 generator(seed()); // initialise mt19937 random number generator with the seed
    std::uniform_int_distribution<> rand_int(0, N-1); // generate a uniformly distributed random integer number between 0 and N-1
    std::uniform_real_distribution<> epsilon(0, 1); // generate a uniformly distributed random integer number between 0 and 1

    int i;
    bool valid = true;

    // choose random particle
    i = rand_int(generator);
    // update position
    p_config_new[i][0] = p_config_old[i][0] + delta*(epsilon(generator)-0.5);
    p_config_new[i][1] = p_config_old[i][1] + delta*(epsilon(generator)-0.5);
    // box check
    if (p_config_new[i][0] < -L_x/2 + sigma/2 ||  p_config_new[i][0] > L_x/2 - sigma/2) {
        valid = false;
    }
    if (p_config_new[i][1] < -L_y/2 + sigma/2 ||  p_config_new[i][1] > L_y/2 - sigma/2) {
        valid = false;
    }
    // check overlap
    for (int j=0; j < N; j++) {
        if (i!=j & metric(p_config_new[i][0],p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]) < sigma) {
            valid = false;
        } 
    }

    if (epsilon(generator) > std::min(1.0,exp(-(g*(p_config_new[i][1]-p_config_old[i][1]))/T))) {
        valid = false;
    }

    return valid;
}

std::string current_time(){
    time_t now = time(NULL);
    struct tm tstruct;
    char buf[40];
    tstruct = *localtime(&now);
    //format: HH:mm:ss
    strftime(buf, sizeof(buf), "%X", &tstruct);
    return buf;
}

int main() {

    // PARAMETERS
    int N = 800; // number of particles
    double sigma = 1; // particle diameter

    double g = 0.5; // gravity force
    double T = 1;

    double phi = 0.1; // area fraction
    const double pi = acos(-1);
    double delta = sigma/3; // random displacement lenght
    double L_x = 15;//sqrt(pi*pow(sigma,2)*N/(16*phi)); // box size
    double L_y = 10*L_x;

    std::cout << "L_x = " << L_x << '\n';

    std::vector< std::vector<double> > p_config_0(N, std::vector<double> (2, 0)); // particle configuration
    std::vector< std::vector<double> > p_config_old(N, std::vector<double> (2, 0)); // particle configuration with PBE
    std::vector< std::vector<double> > p_config_new(N, std::vector<double> (2, 0)); // particle configuration for the new step with PBE

    double time, est_time;;
    std::clock_t time0;

//  Generate initial configuration
    generate_config(N, L_x, L_y, sigma, p_config_0);

//  Monte-Carlo simulation
    p_config_old = p_config_0;
    p_config_new = p_config_0;
    int MCS = 50000;
    int Iini = 25000;
    bool valid;

    std::ofstream file;
    file.open("evolution.dat"); 

    // Monte-Carlo
    for (int IMC=0; IMC<MCS; IMC++){

        // save initialization time
        if (IMC == 0) time0 = std::clock();

        // save configuration
        if ((IMC)%((MCS-Iini)/200) == 0 & IMC > Iini) {
            for (int i=0; i<N; i++) {
                file << p_config_new[i][0] << " " << p_config_new[i][1] << '\n';
            }
            file << '\n';
            file << '\n';
        }

        // 1 Monte-Carlo steps
        for (int IPAS=0; IPAS<N; IPAS++) {
            // generate a Monte-Carlo step
            valid = metropolis(N,L_x,L_y,T,sigma,g,delta,p_config_old,p_config_new);
            // update accordingly
            if (valid==true) p_config_old = p_config_new;
            if (valid==false) p_config_new = p_config_old;
        }

        // estimate cpu time
        if (IMC == 0) {
                est_time = MCS * (std::clock()-time0)/CLOCKS_PER_SEC;
                time = est_time;
                std::cout << "Monte-Carlo CPU ESTIMATED TIME: " << (int)time/60 << " min " << (int)((time/60-(int)(time/60))*60) << " sec" << " at " << current_time() << "\n";
        }

        if (IMC%1 == 0 and IMC != 0) {
            time = est_time - (std::clock()-time0)/CLOCKS_PER_SEC;
            std::cout << " IMC = " << IMC << ", ESTIMATED TIME REMAINING: " << (int)time/60 << " min " << (int)((time/60-(int)(time/60))*60) << " sec   " << "\r";
            std::cout.flush();
        }
    }
    file.close();

    time = (std::clock()-time0)/CLOCKS_PER_SEC;
    std::cout << '\r';
    std::cout.flush();
    std::cout << "Monte-Carlo CPU TIME: " << (int)time/60 << " min " << (int)((time/60-(int)(time/60))*60) << " sec" << " at " << current_time() << "      \n";

    return 0;
}
