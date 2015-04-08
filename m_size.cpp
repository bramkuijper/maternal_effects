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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "bramauxiliary.h"

//#define NDEBUG
//#define DISTRIBUTION

using namespace std;

const int NumGen = 50000;
const int Npop = 5000;
const int skip = 10;
int Noffspring = 0;
int N = 0;
int Nsurv = 0;
int total_offspring = 0;
bool do_stats = 0;

double epsilon = 0; // the value of the environment
double epsilon_sens = 0; // the perceived value of the environment

double size_min = 1.0;
double init_size = 0.1;

int offspring_control = 0;

double theta = 0; 
double omega_b_2 = 0; // width of the selection function
double omega_m_2 = 0; // width of the selection function
double wmin = 1.0; // width of the selection function
double envt_intercept = 0.5; // width of the selection function
double sigma_e = 1.0; // developmental noise
double sigma_ksi = 0.1; // variance of the autocorrelated process
double rho_t = 0.5; // temporal autocorrelation
double freq = 0; // temporal autocorrelation
double mu_g 	  = 0.05;            // mutation rate
double sdmu_g         = 0.4;			 // standard deviation mutation size
double mu_m 	  = 0.05;            // mutation rate
double sdmu_m         = 0.4;			 // standard deviation mutation size
double mu_b 	  = 0.05;            // mutation rate
double sdmu_b         = 0.4;			 // standard deviation mutation size
double ksi = 0;			 // standard deviation mutation size
double ampl = 2.0;
double tau = 0.0;
double mean_surv = 0.0;
double total_surv = 0.0;

const int nloci_b = 2; // number of alleles underlying genetic architecture
const int nloci_g = 2; // number of alleles underlying genetic architecture
const int nloci_m = 2; // number of alleles underlying genetic architecture

int generation = 0;
unsigned seed = 0;



// gsl initialization
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

// the individual struct
struct Individual
{
    // assuming haploid inheritance 
    // and 50 loci coding for g to get proper values of 
    // the genetic variance covariance matrix
    double g[nloci_g]; 
    double b[nloci_b]; 
    double m[nloci_m];
    double phen;
    double phen_m;
    double phen_b;
    double phen_g;
    double cumsurv;
    double surv;
    int generation;
};

typedef Individual Population[Npop];
typedef Individual Population2[Npop*500];
Population Pop;
Population2 Offspring;

string filename("sim_evolving_m");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

#ifdef DISTRIBUTION
string filename_new2(create_filename("sim_evolving_m_dist"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION

// get some of the parameter values
// as command line arguments
void initArguments(int argc, char *argv[])
{
	mu_g = atof(argv[1]);
	mu_m = atof(argv[2]);
	mu_b = atof(argv[3]);
	sdmu_g = atof(argv[4]);
	sdmu_m = atof(argv[5]);
	sdmu_b = atof(argv[6]);
	freq = atof(argv[7]);
	sigma_e = sqrt(atof(argv[8]));
	sigma_ksi = sqrt(atof(argv[9]));
	wmin = atof(argv[10]);
	rho_t = atof(argv[11]);
	omega_b_2 = atof(argv[12]);
	omega_m_2 = atof(argv[13]);
	tau = atof(argv[14]);
	envt_intercept = atof(argv[15]);
    init_size = atof(argv[16]);
    offspring_control = atoi(argv[17]);
}

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
		<< "type:;" << "evol_m_size" << ";" << endl
		<< "control:;" << offspring_control_text << ";" << endl
        << "size_min:;" << size_min << ";" << endl
        << "mu_g:;" << mu_g << ";" << endl
        << "mu_m:;" << mu_m << ";" << endl
        << "mu_b:;" << mu_b << ";" << endl
        << "sdmu_g:;" << sdmu_g << ";" << endl
        << "sdmu_m:;" << sdmu_m << ";" << endl
        << "sdmu_b:;" << sdmu_b << ";" << endl
        << "omega_b_2:;" << omega_b_2 << ";" << endl
        << "omega_m_2:;" << omega_m_2 << ";" << endl
        << "wmin:;" << wmin << ";" << endl
        << "ampl:;" << ampl << ";" << endl
        << "freq:;" << freq << ";" << endl
        << "sigma_e:;" << sigma_e << ";" << endl
        << "sigma_ksi:;" << sigma_ksi << ";" << endl
        << "rho_t:;" << rho_t << ";" << endl
        << "tau:;" << tau << ";" << endl
        << "envt_intercept:;" << envt_intercept << ";" << endl
		<< "seed:;" << seed << ";"<< endl;
}


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
    double ssphen = 0;

    // get stats from the population
    for (int i =  0; i < N; ++i)
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

        ssphen += Pop[i].phen * Pop[i].phen;
    }

    DataFile << generation << ";" << epsilon << ";" << total_offspring << ";" << Nsurv << ";" << mean_surv << ";" << ksi << ";";

    DataFile 
            << (meanphen/N) << ";"
            << (ssphen/N - pow(meanphen/N,2.0)) << ";"
            << (meanphen_m/N) << ";"
            << (meang/(N)) << ";"
            << (ssg/(N) - pow(meang/N,2.0)) << ";"
            << (meanm/(N)) << ";"
            << (ssm/N - pow(meanm/N,2.0)) << ";" 
            << (meanb/N) << ";"
            << (ssb/N - pow(meanb/N,2.0)) << ";"  << endl;
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
        Pop[i].phen = init_size;
        Pop[i].phen_m = init_size;

        for (int j = 0; j < nloci_g; ++j)
        {
            Pop[i].g[j] = .5*init_size;
        }

        for (int j = 0; j < nloci_b; ++j)
        {
            Pop[i].b[j] = 0;
        }

        for (int j = 0; j < nloci_m; ++j)
        {
            Pop[i].m[j] = 0;
        }

        Pop[i].generation = generation - 1;
	}

    epsilon = envt_intercept+envt_intercept*sin(freq * (generation)) + ksi;
    N = Npop;
}

// create an offspring
void Create_Kid(int mother, int father, Individual &kid)
{
    double sum_g = 0; // sum over all the breeding values of the offspring coding for the actual phenotype
    double sum_b = 0; // sum over all the breeding values of the offspring coding for the norm of reaction
    double sum_m = 0; // sum over all the breeding values of the offspring coding for the maternal effect

    assert(mother >= 0 && mother <= N);
    assert(father >= 0 && father <= N);

    // we assume all loci are unlinked (i.e., parent of origin (poi)
    // of one allele is independent of poi of any other alles
    for (int i = 0; i < nloci_g;++i)
    {
        kid.g[i] = i % 2 == 0 ? Pop[mother].g[i + gsl_rng_uniform_int(r, 2)] : Pop[father].g[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateG(kid.g[i]);
        sum_g += kid.g[i];
    }

    for (int i = 0; i < nloci_b; ++i)
    {
        kid.b[i] = i % 2 == 0 ? Pop[mother].b[i + gsl_rng_uniform_int(r, 2)] : Pop[father].b[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateB(kid.b[i]);
        sum_b += kid.b[i];
    }

    for (int i = 0; i < nloci_m; ++i)
    {
        kid.m[i] = i % 2 == 0 ? Pop[mother].m[i + gsl_rng_uniform_int(r, 2)] : Pop[father].m[i - 1 + gsl_rng_uniform_int(r, 2)];
        MutateM(kid.m[i]);
        sum_m += kid.m[i];
    }

    kid.phen_m = sum_m;
    kid.phen_g = sum_g;
    kid.phen_b = sum_b;
    kid.generation = generation;

    if (offspring_control == 0)
    {
        kid.phen = kid.phen_g + gsl_ran_gaussian(r,sigma_e) + kid.phen_b * epsilon_sens + kid.phen_m * Pop[mother].phen;
    } 
    else if (offspring_control == 1)
    {
        kid.phen = kid.phen_g + gsl_ran_gaussian(r,sigma_e) + kid.phen_b * epsilon_sens + Pop[mother].phen_m * Pop[mother].phen;
    }
    else
    {
        kid.phen = Pop[mother].phen_g + gsl_ran_gaussian(r,sigma_e) + Pop[mother].phen_b * epsilon_sens + Pop[mother].phen_m * Pop[mother].phen;
    }

    assert(isnan(kid.phen) == 0);
}


// Survival of juveniles to reproductive adults
void Reproduce()
{
    double c = epsilon;

    assert(c >= 0);

    Noffspring = 0;
    int Noffspring_total = 0;

    total_surv = 0;

    mean_surv = 0;

    double mean_size = 0;

    // allow all individuals to reproduce
    // as they are hermaphrodites
    for (int i = 0; i < N; ++i)
    {
        double clutch_size = 200;

        int father = -1;

        // choose a father (no selfing)
        do {
            father = gsl_rng_uniform_int(r, N);
        } while (father == i);

        int nkids_missed = 0;

        // create offspring
        do {
            Individual Kid;

            Create_Kid(i, father, Kid);

            Kid.surv = wmin + 
                (1.0 -  wmin) * exp(-.5 * (pow(Kid.phen_b,2.0)/omega_b_2
                                                + pow(Kid.phen_m,2.0)/omega_m_2))
                * (1.0 - exp(-c*(Kid.phen - size_min)));
            

            // remove inviable offspring        
            // if parents consistently produce those offspring
            // end their reproduction
            if (Kid.phen < size_min)
            {
                ++nkids_missed;
                if (nkids_missed > 10)
                {
                    break;
                }
                    
                continue;
            }

            if (Kid.phen <= clutch_size)
            {
                clutch_size -= Kid.phen;
        
                mean_surv += Kid.surv;
                ++Noffspring_total;


                if (gsl_rng_uniform(r) < Kid.surv)
                {
                    Offspring[Noffspring++] = Kid;
                }
            } 
            else if (gsl_rng_uniform(r) < clutch_size/Kid.phen)
            {
                mean_surv += Kid.surv;
                ++Noffspring_total;
                // one last offspring from the remainder of resources...
                //
                if (gsl_rng_uniform(r) < Kid.surv)
                {
                    Offspring[Noffspring++] = Kid;
                }
                break;
            }
            else
            {
                // out of resources
                break;
            }
        }
        while (clutch_size > 0);
    }

    mean_surv /= Noffspring_total;

    if (Noffspring <= 1)
    {
        cout << (total_surv == 0 ? "no survivors. " : "no offspring produced. ") << endl;
        WriteData();
        exit(1);
    }

    total_offspring = Noffspring;

}


void Survive()
{
    int random_offspring;

    N = Noffspring < Npop ? Noffspring : Npop;

    for (int i = 0; i < N; ++i)
    {
        random_offspring = gsl_rng_uniform_int(r, Noffspring);
        Pop[i] = Offspring[random_offspring];

        Offspring[random_offspring] = Offspring[Noffspring-1];
        --Noffspring;
    }

    epsilon = envt_intercept+envt_intercept*sin(freq * (generation+1)) + ksi;

    if (tau > 0)
    {
        int timesteps = rint(1.0 / tau);

        for (int time_i = 0; time_i < timesteps; ++time_i)
        {
            ksi = rho_t*ksi + gsl_ran_gaussian(r, sqrt(1.0-rho_t*rho_t)*sigma_ksi);
        }
    }

    epsilon_sens = sin(freq * (generation-tau+1)) + ksi;
}


void WriteDataHeaders()
{
    DataFile << "generation;epsilon;noffspring;nsurv;mean_surv;ksi;meanz;varz;meanphen_m;meang;varg;meanm;varm;meanb;varb;" << endl;
}

int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
        do_stats = generation % skip == 0;

		Reproduce();
		Survive();

        if (do_stats)
		{
			WriteData();
		}
	}

	WriteParameters();
}
