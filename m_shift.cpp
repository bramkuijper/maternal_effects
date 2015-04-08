// the evolution of maternal effects in a sinusoidal environment
// 
// Bram Kuijper & Rebecca B. Hoyle
//
// this code is published according to the GNU Public license v3
// https://www.gnu.org/licenses/gpl.html 
//
// Kuijper, B & Hoyle, R. B. 2015
// When to rely on maternal effects and when to rely on phenotypic plasticity?
// Evolution, http://dx.doi.org/10.1111/evo.12635 
//
// You may also find the following paper interesting:
// Kuijper, B.; Johnstone, R. A. & Townley, S. (2014). The evolution of 
// multivariate maternal effects. PLoS Comp. Biol. 10: e1003550. 
// http://dx.doi.org/10.1371/journal.pcbi.1003550
//  


#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>

// random number generation
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// various functions, such as unique filename creation
#include "bramauxiliary.h"

//#define NDEBUG
//
// the compilation sign should only be turned on when one wants
// to assess the complete distribution of phenotypes
// 
//#define DISTRIBUTION

using namespace std;

// number of generations
const int NumGen = 50000;

// population size
const int Npop = 5000; 

// number of generations to skip when outputting data
const int skip = 10;

// track number of survivors
int NSurv = 0;

// indicator variable if we are printing stats for this generation
bool do_stats = 0;

double epsilon = 0; // the value of the environment
double epsilon_sens = 0; // the perceived value of the environment


double theta = 0;  // phenotypic optimum
double omega2 = 0; // width of the selection function
double omega_b_2 = 0; // width of the selection function for plasticity
double omega_m_2 = 0; // width of the selection function for maternal effects
double wmin = 0.0; // minimal survival probability
double sigma_e = 1.0; // variance of developmental noise
double sigma_ksi = 0.1; // variance of the autocorrelated process
double rho_t = 0.5; // temporal autocorrelation
double mu_g 	  = 0.05;            // mutation rate
double sdmu_g         = 0.05;			 // standard deviation mutation size
double mu_m 	  = 0.05;            // mutation rate
double sdmu_m         = 0.05;			 // standard deviation mutation size
double mu_b 	  = 0.05;            // mutation rate
double sdmu_b         = 0.05;			 // standard deviation mutation size
double ksi = 0;			 // standard deviation mutation size
double A = 0.0; // environmental intercept
double B = 2.0; // environmental amplitude of change
double delta_U = 10; // steepness of change
double tau = 0.0;       // developmental time lag
double init_g = 0; // initial values for g,m,b
double init_m = 0;
double init_b = 0;

const int n_alleles_b = 2; // number of alleles underlying genetic architecture
const int n_alleles_g = 2; // number of alleles underlying genetic architecture
const int n_alleles_m = 2; // number of alleles underlying genetic architecture

int offspring_control = 0; // control over genetic loci, offspring vs mother

// keep track of the current generation number
int generation = 0;

// random seed
unsigned seed = 0;

// gnu scientific library random number generator initialization
// http://www.gnu.org/software/gsl/ 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual struct
struct Individual
{
    double g[n_alleles_g]; 
    double b[n_alleles_b]; 
    double m[n_alleles_m];
    double phen; // an individual's phenotype, z
    double phen_m; //m
    double phen_b; // b
    double phen_g; // g
};

// allocate a population and a population of survivors
typedef Individual Population[Npop];
Population Pop;
Population Survivors;

// generate a unique filename for the output file
string filename("sim_evolving_m");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  

#ifdef DISTRIBUTION
// generate a filename for the phenotype distribution file
string filename_new2(create_filename("sim_evolving_m_dist"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION

// initialize simulations from command line arguments
void initArguments(int argc, char *argv[])
{
	omega2 = atof(argv[1]);
	mu_g = atof(argv[2]);
	mu_m = atof(argv[3]);
	mu_b = atof(argv[4]);
	sdmu_g = atof(argv[5]);
	sdmu_m = atof(argv[6]);
	sdmu_b = atof(argv[7]);
    B = atof(argv[8]);
	sigma_e = sqrt(atof(argv[9]));
	sigma_ksi = sqrt(atof(argv[10]));
	wmin = atof(argv[11]);
	rho_t = atof(argv[12]);
	omega_b_2 = atof(argv[13]);
	omega_m_2 = atof(argv[14]);
	tau = atof(argv[15]);
    init_g = atof(argv[16]);
    init_m = atof(argv[17]);
    init_b = atof(argv[18]);
    offspring_control = atof(argv[19]);
}


// mutation according to a continuum of alleles model
void MutateG(double &G)
{
	G += gsl_rng_uniform(r)<mu_g ? gsl_ran_gaussian(r, sdmu_g) : 0;
}

void MutateM(double &G)
{
	G += gsl_rng_uniform(r)<mu_m ? gsl_ran_gaussian(r,sdmu_m) : 0;
}

void MutateB(double &G)
{
	G += gsl_rng_uniform(r)<mu_b ? gsl_ran_gaussian(r,sdmu_b) : 0;
}

// write the parameters (typically at the end of the output file)
void WriteParameters()
{
    string offspring_control_text = "";

    switch(offspring_control)
    {
        case 0:
            offspring_control_text = "offspring_control_g_m";
            break;
        case 1:
            offspring_control_text = "offspring_control_g_maternal_control_m";
            break;
        case 2:
            offspring_control_text = "maternal_control_g_m";
            break;
        default:
            break;
    }

	DataFile << endl
		<< endl
		<< "type:;" << "evolve_m_sinusoidal" << ";" << endl
		<< "control:;" << offspring_control_text << ";" << endl
        << "mu_g:;" << mu_g << ";" << endl
        << "mu_m:;" << mu_m << ";" << endl
        << "mu_b:;" << mu_b << ";" << endl
        << "sdmu_g:;" << sdmu_g << ";" << endl
        << "sdmu_m:;" << sdmu_m << ";" << endl
        << "sdmu_b:;" << sdmu_b << ";" << endl
        << "omega2:;" << omega2 << ";" << endl
        << "omega_b_2:;" << omega_b_2 << ";" << endl
        << "omega_m_2:;" << omega_m_2 << ";" << endl
        << "wmin:;" << wmin << ";" << endl
        << "init_g:;" << init_g << ";" << endl
        << "init_m:;" << init_m << ";" << endl
        << "init_b:;" << init_b << ";" << endl
        << "A:;" << A << ";" << endl
        << "B:;" << B << ";" << endl
        << "delta_U:;" << delta_U << ";" << endl
        << "sigma_e:;" << sigma_e << ";" << endl
        << "sigma_ksi:;" << sigma_ksi << ";" << endl
        << "rho_t:;" << rho_t << ";" << endl
        << "tau:;" << tau << ";" << endl
		<< "seed:;" << seed << ";"<< endl;
}

// initialize the simulation
// by giving all the individuals 
// genotypic values
//
// and doing some other stuff (e.g., random seed)
void Init()
{
    // get the timestamp (with nanosecs)
    // to initialize the seed
	seed = get_nanoseconds();
    
    // set the seed to the random number generator
    // stupidly enough, for gsl this can only be done by setting
    // a shell environment parameter
    stringstream s;
    s << "GSL_RNG_SEED=" << setprecision(10) << seed;
    putenv(const_cast<char *>(s.str().c_str()));

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);


	// initialize the whole populatin
	for (int i = 0; i < Npop; ++i)
	{
        Pop[i].phen = 0;
        Pop[i].phen_m = 0;

        for (int j = 0; j < n_alleles_g; ++j)
        {
            Pop[i].g[j] = init_g/n_alleles_g;
        }

        for (int j = 0; j < n_alleles_b; ++j)
        {
            Pop[i].b[j] = init_b/n_alleles_b;
        }

        for (int j = 0; j < n_alleles_m; ++j)
        {
            Pop[i].m[j] = init_m/n_alleles_m;
        }
	}
}

// create an offspring
void Create_Kid(int mother, int father, Individual &kid)
{
    double sum_g = 0; // sum over all the breeding values of the offspring coding for the actual phenotype
    double sum_b = 0; // sum over all the breeding values of the offspring coding for the norm of reaction
    double sum_m = 0; // sum over all the breeding values of the offspring coding for the maternal effect

    // we assume all loci are unlinked 
    for (int i = 0; i < n_alleles_g;++i)
    {
        kid.g[i] = i % 2 == 0 ? Survivors[mother].g[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].g[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateG(kid.g[i]);
        sum_g += kid.g[i];
    }

    for (int i = 0; i < n_alleles_b; ++i)
    {
        kid.b[i] = i % 2 == 0 ? Survivors[mother].b[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].b[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateB(kid.b[i]);
        sum_b += kid.b[i];
    }

    for (int i = 0; i < n_alleles_m; ++i)
    {
        kid.m[i] = i % 2 == 0 ? Survivors[mother].m[i + gsl_rng_uniform_int(r, 2)] : Survivors[father].m[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateM(kid.m[i]);
        sum_m += kid.m[i];
    }

    kid.phen_m = sum_m;
    kid.phen_g = sum_g;
    kid.phen_b = sum_b;

    // phenotype determination according to Hoyle & Ezard 2012 Interface
    if (offspring_control == 0)
    {
        // complete offspring control over a,b,m
        kid.phen = kid.phen_g + gsl_ran_gaussian(r,sigma_e) + kid.phen_b * epsilon_sens + kid.phen_m * Survivors[mother].phen;
    } 
    else if (offspring_control == 1)
    {
        // offspring control over a,b but not m
        kid.phen = kid.phen_g + gsl_ran_gaussian(r,sigma_e) + kid.phen_b * epsilon_sens + Survivors[mother].phen_m * Survivors[mother].phen;
    }
    else
    {
        // complete maternal control
        kid.phen = Survivors[mother].phen_g + gsl_ran_gaussian(r,sigma_e) + Survivors[mother].phen_b * epsilon_sens + Survivors[mother].phen_m * Survivors[mother].phen;
    }

    assert(isnan(kid.phen) == 0);
}


// Survival of juveniles to reproductive adults
void Survive()
{
    double W;

    // shift after 50000 generations
    if (generation >= 50000)
    {
        delta_U = 10;
    }


    double theta = A + B * epsilon;

    NSurv = 0;

    for (int i = 0; i < Npop; ++i)
    {
        W = wmin + (1.0 -  wmin) * exp(-.5 * (
                                            pow((Pop[i].phen - theta),2.0)/omega2 
                                            + pow(Pop[i].phen_b,2.0)/omega_b_2
                                            + pow(Pop[i].phen_m,2.0)/omega_m_2
                                        )
                                    );


        assert(isnan(W) == 0);

        if (gsl_rng_uniform(r) < W)
        {
            Survivors[NSurv++] = Pop[i];
        }
    }
    
    if (NSurv == 0)
    {
        WriteParameters();
        exit(1);
    }

    // update the environment for the next generation
    // as an autocorrelated gaussian random variable
    ksi = rho_t*ksi + gsl_ran_gaussian(r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);

    epsilon = delta_U + ksi;

    // in the likely case there is a developmental timelag, tau,
    // update the environment for a number of 'sub' timesteps
    // to achieve a 'sensed' (rather than real) value of epsilon
    if (tau > 0)
    {
        int timesteps = rint(1.0 / tau);

        for (int time_i = 0; time_i < timesteps; ++time_i)
        {
            ksi = rho_t*ksi + gsl_ran_gaussian(r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);
        }
    }

    // update the value of the sensed environment
    epsilon_sens = delta_U + ksi;

    for (int i = 0; i < Npop; ++i)
    {
        Individual Kid;

        Create_Kid(gsl_rng_uniform_int(r,NSurv), gsl_rng_uniform_int(r,NSurv), Kid);

        Pop[i] = Kid;
    }
}


// write down summary statistics
void WriteData()
{
    double meanphen = 0;
    double meanphen_m = 0;
    double meang = 0;
    double ssg = 0;
    double meanm = 0;
    double ssm = 0;
    double meanb = 0;
    double ssb = 0;

    // get stats from the population
    for (int i =  0; i < Npop; ++i)
    {
        // stats for m
        meang += Pop[i].phen_g;
        ssg += Pop[i].phen_g * Pop[i].phen_g;

        // stats for m
        meanm += Pop[i].phen_m;
        ssm += Pop[i].phen_m * Pop[i].phen_m;

        meanb += Pop[i].phen_b;
        ssb += Pop[i].phen_b * Pop[i].phen_b;

        meanphen += Pop[i].phen;
        meanphen_m += Pop[i].phen_m;
    }

    DataFile << generation << ";" << epsilon << ";" << NSurv << ";" << ksi << ";";

    DataFile 
            << (meanphen/Npop) << ";"
            << (meanphen_m/Npop) << ";"
            << (meang/(Npop)) << ";"
            << (ssg/(Npop) - pow(meang/Npop,2.0)) << ";"
            << (meanm/(Npop)) << ";"
            << (ssm/Npop - pow(meanm/Npop,2.0)) << ";" 
            << (meanb/Npop) << ";"
            << (ssb/Npop - pow(meanb/Npop,2.0)) << ";"  << endl;
}

// write the headers of a datafile
void WriteDataHeaders()
{
    DataFile << "generation;epsilon;nsurv;ksi;meanz;meanphen_m;meang;varg;meanm;varm;meanb;varb;" << endl;
}


// the guts of the code
int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		Survive();

        if (do_stats)
		{
			WriteData();
		}
	}

	WriteParameters();
}
