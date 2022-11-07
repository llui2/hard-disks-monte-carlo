#include <iostream> // terminal input-output
#include <random> // random numbers
#include <vector> // vectors
#include <fstream> // files input-output
#include <string> // strings
#include <ctime> // calculate cpu times
#include <math.h>  // modf

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

void generate_triangular_config(int N, double L, double sigma,std::vector< std::vector<double> > &p_config) {
    /* This function generates an initial configuration of N solid discs of radius sigma
    in a LxL box arranged in a triangular lattice */

    // trick to ajust properly for every N
    double M = sqrt(N);
    double Mint, Mrest;
    int l = 0;
resize:
    Mrest = modf(M, &Mint);
    if (Mrest != 0) {
        l += 1;
        M = sqrt(N+l);
        goto resize;
    }

    double d = sigma;
    int i = 0;
    int j = 0;
    double t = 0;

    for (int k=0; k<N; k++) {
        p_config[k][0] = -M/2*sigma + d*i;
        p_config[k][1] = -M/2*sigma + d*j + t*d/2;
        
        if (i == M-1) {
            j += 1;
            i = 0;
        } else i += 1;

        if (i%2 != 0) t = 1;
        if (i%2 == 0) t = 0;
    }
}

bool metropolis(double N, double L, double sigma, double delta,
    std::vector< std::vector<double> > &p_config_old, std::vector< std::vector<double> > &p_config_new) {
    /* This function does a Monte-Carlo step */
    
    std::random_device seed; // generate a random seed (random integer)
    std::mt19937 generator(seed()); // initialise mt19937 random number generator with the seed
    std::uniform_int_distribution<> rand_int(0, N-1); // generate a uniformly distributed random integer number between 0 and N-1
    std::uniform_int_distribution<> epsilon(0, 1); // generate a uniformly distributed random integer number between 0 and 1

    int i;
    bool valid = true;
    // choose random particle
    i = rand_int(generator);
    // update position
    p_config_new[i][0] = p_config_old[i][0] + delta*(epsilon(generator)-0.5);
    p_config_new[i][1] = p_config_old[i][1] + delta*(epsilon(generator)-0.5);
    // PBC check
    if (p_config_new[i][0] < -L/2 ||  p_config_new[i][0] > L/2) {
        p_config_new[i][0] += L*sign(p_config_new[i][0]);
    }
    if (p_config_new[i][1] < -L/2 ||  p_config_new[i][1] > L/2) {
        p_config_new[i][1] += L*sign(p_config_new[i][1]);
    }
    // check overlap with PBC
    for (int j=0; j < N; j++) {
        if (i!=j & metric(p_config_new[i][0],p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]) < sigma) {
            valid = false;
        } else if (i!=j & metric(p_config_new[i][0],p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]) > L-sigma) {
            if (metric(p_config_new[i][0],p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]+L) < sigma) valid = false;
            if (metric(p_config_new[i][0],p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]-L) < sigma) valid = false;
            if (metric(p_config_new[i][0],p_config_new[j][0]+L,p_config_new[i][1],p_config_new[j][1]) < sigma) valid = false;
            if (metric(p_config_new[i][0],p_config_new[j][0]-L,p_config_new[i][1],p_config_new[j][1]) < sigma) valid = false;

            if (metric(p_config_new[i][0]+L,p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]+L) < sigma) valid = false;
            if (metric(p_config_new[i][0]+L,p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]-L) < sigma) valid = false;
            if (metric(p_config_new[i][0]-L,p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]+L) < sigma) valid = false;
            if (metric(p_config_new[i][0]-L,p_config_new[j][0],p_config_new[i][1],p_config_new[j][1]-L) < sigma) valid = false;
        }
    }
    return valid;
}

int main() {

    std::random_device seed; // generate a random seed (random integer)
    std::mt19937 generator(seed()); // initialise mt19937 random number generator with the seed
    std::uniform_real_distribution<float> uniform_one(-1, 1); // function to generate a uniform real random number between -1 and 1

    // PARAMETERS
    int N = 1000; // number of particles
    double sigma = 1; // diameter
    double phi = 0.5; // area fraction (critical value 0.78)
    const double pi = acos(-1);
    //double delta[] = {sigma*0.001,sigma*0.003,sigma*0.01,sigma*0.03,sigma*0.1,sigma*0.3}; // random displacement lenght
    double delta = 0.3; // random displacement lenght
    double L = sqrt(pi*pow(sigma,2)*N/(4*phi)); // box size

    std::vector< std::vector<double> > p_config_0(N, std::vector<double> (2, 0)); // particle configuration
    std::vector< std::vector<double> > p_config_old(N, std::vector<double> (2, 0)); // particle configuration with PBE
    std::vector< std::vector<double> > p_config_new(N, std::vector<double> (2, 0)); // particle configuration for the new step with PBE

    double time;
    std::clock_t time0;

//  Generate initial configuration
    generate_triangular_config(N, L, sigma, p_config_0);

//  Monte-Carlo simulation
    std::uniform_int_distribution<> rand_int(0, N-1); // generate a uniformly distributed random integer number between 0 and N-1
    std::uniform_int_distribution<> epsilon(0, 1); // generate a uniformly distributed random integer number between 0 and 1
    p_config_old = p_config_0;
    p_config_new = p_config_0;
    int MCS = 3000;
    bool valid;

    double MSD; // mean-squared displacement
    double dx; // variable to store x_f - x_o
    double dy; // variable to store y_f - y_o

    std::ofstream file;
    file.open("evolution.dat"); 
    std::ofstream file1;
    file1.open("data_3.dat"); 

    // Monte-Carlo
    for (int IMC=0; IMC<MCS; IMC++){

        // save initialization time
        if (IMC == 0) time0 = std::clock();

        // save configuration
        if (IMC%(MCS/250) == 0) {
            for (int i=0; i<N; i++) {
                file << p_config_new[i][0] << " " << p_config_new[i][1] << '\n';
            }
            file << '\n';
            file << '\n';
        }

        // 1 Monte-Carlo steps
        for (int IPAS=0; IPAS<N; IPAS++) {
            // generate a Monte-Carlo step
            valid = metropolis(N,L,sigma,delta,p_config_old,p_config_new);
            // update accordingly
            if (valid==true) p_config_old = p_config_new;
            if (valid==false) p_config_new = p_config_old;
        }

        // save MSD
        if (IMC%1 == 0) {
            MSD = 0;
            for (int i=0; i<N; i++) {
                // compute MSD with PBC
                dx = 0;
                dy = 0;
                dx = abs(p_config_new[i][0] - p_config_0[i][0]);
                if (dx > L/2) dx -= L;
                dy = abs(p_config_new[i][1] - p_config_0[i][1]);
                if (dy > L/2) dy -= L;
                MSD += pow(dx,2) + pow(dy,2);
            }
            file1 << IMC << ' ' << MSD/N << '\n';
        }

        // estimate cpu time
        if (IMC == 0) {
                time = MCS * (std::clock()-time0)/CLOCKS_PER_SEC;
                std::cout << "Monte-Carlo CPU ESTIMATED TIME: " << (int)time/60 << " min " << (int)((time/60-(int)(time/60))*60) << " sec" << "\n";
        }
    }
    file.close();
    file1.close();    

    time = (std::clock()-time0)/CLOCKS_PER_SEC;
    std::cout << "Monte-Carlo CPU TIME: " << (int)time/60 << " min " << (int)((time/60-(int)(time/60))*60) << " sec" << "\n";
    
    std::string str = "gnuplot figure_3.gnu";
    std::cout << str << '\n';
    system(str.c_str());

    std::string radius = std::to_string(sigma/2);
    std::string lenght = std::to_string(L);
    std::string str1 = "gnuplot -e L=" + lenght + " -e r=" + radius + " evolution.gnu";
    std::cout << str1 << '\n';
    system(str1.c_str());

    return 0;
}


