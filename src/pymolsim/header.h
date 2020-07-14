#ifndef HEADER_H
#define HEADER_H

//-------------------------------------- other header files------------------------------
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <vector>
#include <float.h>
#include "num2str.h"

//-------------------------------------------workspace------------------------------------
using namespace std;


//---------------------------------------some useful definitions--------------------------
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN(a,b,tmp)    (  (tmp=a) < b  ? tmp : b )
#define MAX(a,b,tmp)    (  (tmp=a) > b  ? tmp : b )
#define update_average(a,b) {a.now = b; a.sum += a.now; a.sumsq += a.now*a.now; a.n++;}


//velocity histograms
#define MBIN_v 1000
#define RESHISTO_v 100

//volume histogram
#define MBIN_Vol 1000
#define RESHISTO_Vol 2

//maximum number of particles
#define MAX_PART 20000

//velocity autocorrelation function
#define VACF_MAXT 1000		//maximum number of steps for the VACF
#define VACF_MAXT0 200		//maximum number of time origins
#define VACF_SAMPLE 4		//no. of time steps in between sampling VACF
#define VACF_FREQT0 50		//how often a new time origin is assigned

//-------------------------------------------classes---------------------------------------

class Potential
{
	public:
		double sigma,
		       epsilon;
		
		double c1,
		       c2,
		       c3,
		       c4,
		       c5;
		
		double rmin2,
		       rmax2;
		
		double virial;		//virial part of the pressure
		double hypervirial;	//hypervirial to calculate thermodynamic properties
		double press_kin;	//kinetic part of the pressure
		double energy;
		vector<double> stress {0, 0, 0, 0, 0, 0, 0, 0, 0};	//instantaneous stress tensor
		vector<double> fi{0., 0., 0.};
		vector<double> fj{0., 0., 0.};
		void forces(vector<double> );
		void potential_energy(vector<double>);
		void init();
};


class Particle
{
	public:
		vector<double> r{0., 0., 0.};
		vector<double> v{0., 0., 0.};
		vector<double> f{0., 0., 0.};
		double energy;
		int type;
		double mass;

		vector<double> c{0., 0.};	//to calculate system state
				//c[0] = c_r
				//c[1] = c_r * c_alpha
		double c_z;	// c_z
		double phi;	//to calculate collective variable

		double vt0[3][VACF_MAXT0];	//velocities at different time origins
};

class Thermostat
{
	public:
		double lgamma,		//Langevin friction
		       lc1,		//Langevin c1 coefficient
		       lc2;		//Langevin c2 coefficient

		double anu;		//Andersen collision frequency

		double nhlgamma,  	//Nose-Hoover-Langevin friction
		       nhlmu, 		//       "             mass
		       nhlc1,  		//       "             c1 coefficient
		       nhlc2,     	//       "             c2 coefficient
		       nhlxi;      	//       "             extended variable
		void init(double, double);
		void init_nhl(double, double, int);
		double dElangevin;
};

class Barostat
{
	public:
		bool isotropic;		//isotropic or anisotropic barostat

		double pv,		//piston velocity
		       pmass;		//piston mass

		double lgamma,		//Langevin friction for piston
		       lc1,		//Langevin c1 coefficient
		       lc2;		//Langevin c2 coefficient
		void init(double, double);
};


class Average
{
	public:
		string name;

		int64_t n,
			in;

		long double now,
		            sum,
		            sumsq;
};


class Statistic
{
	public:
		Average Ekin,		//kinetic energy
			T,
			p[3];		//momentum
		
		//constructor
		Statistic();
};


class System
{
	public:
		int64_t ncycle1,
			ncycle2;
		
		int swcount,
		    swplot,	//plot or not
		    plotrate,	//rate of plotting
		    outputrate;	//rate of writing output

		int nparticles;	//number of particles
		int readstruc;	//read in or create structure
		double E0;	//initial energy


		double rho;	//density
		vector<double> box_relative {0, 0, 0};	//relative box dimension

		double beta,		//inverse temperatur (input)
		       pressure,	//pressure (input)
		       dt,		//time step
		       time,		//total time
		       tot_energy;	//total energy for NVE

		double pe;
		double ke;
		vector<double> p{0., 0., 0.};


			//negative sum of increments of energy contributions
					//from Langevin thermostat

		double Tinst,		//instantaneous temperature
		       Pinst,		//instantaneous pressure
		       Pinst_1,		//instantaneous pressure without T fluctuations!
		       virial,		//virial part of the pressure
		       hypervirial,	//hypervirial to calculate thermodynamic properties
		       press_kin;	//kinetic part of the pressure
		vector<double> stress {0, 0, 0, 0, 0, 0, 0, 0, 0};	//instantaneous stress tensor

		long double c_v,		//heat capacity
		       c_p,		//heat capacity
		       gamma_v,		//thermal pressure coefficient
		       beta_S,		//adiabatic compressibility
		       beta_T,		//isothermal compressibility
		       alpha_p;		//thermal expansion coefficient

		
		vector<double> box {0, 0, 0};		//boxlength in x,y,z
		vector<double> box_2 {0, 0, 0};	// 1/2 of the boxlenght in x,y,z
		double volume;


		int integrator;		//define integration scheme
					// 0 = velocity verlet
					// 1 = langevin thermostat


		int64_t *rdf,		//array to sample rdf
			rdfcount;		//counter for rdf

		int rdfBIN,		//number of bins for the rdf
		    rdfRES;		//resolution of rdf;  0 < r < rdfBIN/rdfRES 
		
		double rmax_rdf;	//maximum distance that can be sampled in rdf

		FILE *fprdf;		//rdf file pointer

		int64_t histo_v[4][MBIN_v],	//velocity distribution
			vcount[4];
		int64_t histo_Vol[MBIN_Vol],	//volume distribution
			Volcount;
		double histo_Vol0;		//initial volume

		double histo_vacf[VACF_MAXT];		//histogram for VACF

		int64_t vacfcount[VACF_MAXT];		//counter for normalisation

		int vacf_indext0,			// index of the corresponding time origin
		    vacf_time,				// 'time' for the vacf
		    vacf_t0time[VACF_MAXT0],	//t0 for different time origins
		    vacf_countt0;			//count no. of time origins


		//holds individual properties
		vector<Particle> particles;
		Potential potential;
		Thermostat thermostat;
		Barostat barostat;
		vector<Average> pAv;
		
		vector<Particle> gparticles();
		void sparticles(vector<Particle> );

		Potential gpotential();
		void spotential(Potential);

		Thermostat gthermostat();
		void sthermostat(Thermostat);

		Barostat gbarostat();
		void sbarostat(Barostat);

		vector<Average> gaverage();
		void saverage(vector<Average>);

		//pyscal integration to read in files
		//will implement on the python side
		//force methods
		void forces();
		void potential_energy();
		void kinetic_energy();
		void total_energy();
		void total_momentum();
		void rescale_velocities();


		//averaging methods
		vector<double> image_distance(int, int);
		void remap();
		void average();
		void sample_average();
		void init_average();

		void init();

		//md methods
		void md_step();
		void md_verlet();
		void md_langevin();
		void md_andersen();
		void md_nosehooverlangevin_NVT();
		void md_andersen_NPH();
		void md_andersen_stochastic_NPT();
		void md_andersen_stochastic_nhlthermo_NPT();
		void langevin_thermo();
		void langevin_baro();
		void propagate_momenta_half();
		void propagate_position_half();
		void propagate_momenta_xi();
		void langevin_xi();
		void propagate_position_half_scale(double *);
		void rescale_position_momenta(double *);



};



//-------------------------------------------functions-------------------------------------------------------

//DEPRECATED
// in readinput.cpp
//void readinput(System &, Potential &, Thermostat &, Barostat &, int &,int argc,char *argv[]);
void readinput(System &, Potential &, Thermostat &, Barostat &, int &);
void read_initstructure(System &,Particle *);
void read_initlammpsdump(System &, Particle *);


//THIS CAN REMAIN THE SAME - oustide use is not required
//in ran3.cpp
double ran3(int=1);
double RandomVelocity(double,double);
double gauss();


//ALL THESE METHODS WOULD BE MOVED TO POTENTIAL/SYSTEM
//in forces.cpp
//specific
void kinetic_stress(Particle *, System &);

//in print.cpp
void print_lammps_dump(System &, Particle *, ofstream &);
void init_print_lammps_trajectory(System &, Particle *, ofstream &);
void print_lammps_trajectory(System &, Particle *, ofstream &);
void print_output(System &, Particle *, Average *, Potential &,int ,int );


//in main.cpp
int main_program();
void terminate_block(System);
void final_block(System);




//-------------------------------------------Global variables------------------------------------------------

extern Statistic block_stats;
extern Statistic final_stats;
extern char g_flname[FILENAME_MAX+1];  //file name
extern bool g_display;		//run with or without display
extern float g_Distance;	//set zoom in display

extern ofstream fpscreen;	//capture screen output

extern FILE *fp[2];		//output file pointer





#endif // HEADER_H

