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
};


class Average
{
	public:
		char *name;

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

		double dElangevin;	//negative sum of increments of energy contributions
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
		vector<Particle> pparticle;
		


};



//-------------------------------------------functions-------------------------------------------------------

// in readinput.cpp
//void readinput(System &, Potential &, Thermostat &, Barostat &, int &,int argc,char *argv[]);
void readinput(System &, Potential &, Thermostat &, Barostat &, int &);
void read_initstructure(System &,Particle *);
void read_initlammpsdump(System &, Particle *);

//in ran3.cpp
double ran3(int=1);
double RandomVelocity(double,double);
double gauss();

//in init.cpp
void init(Particle *,System,Potential);
void init_ene(Particle *,System &,Potential &);
void relax_steep(Particle *,System &,Potential &);
void init_structure(System &, Particle *, Potential &);
void init_lj(Potential &);
void init_langevin(Thermostat &, System &);
void init_langevin_piston(Barostat &, System);
void init_nhl(Thermostat &, System);


//in forces.cpp
void forces(Particle *,System &,Potential &);
double potential_energy(Particle *,System &,Potential &);
double kinetic_energy(Particle *,System &);
double total_energy(Particle *,System &,Potential &);
void total_momentum(Particle *,double *,System &);
void rescale_velocities(Particle *,System &);
void kinetic_stress(Particle *, System &);


//in analyse.cpp
void image_distance(Particle &,Particle &,double *,System &);
void remap(System &, Particle *);
void average(System &, Particle *);
void sample_average(System &, Particle *, Average *, Potential &);
void init_average(Average *);

void init_rdf(System &);
void sample_rdf(System &, Particle *);
void print_rdf(System &);
void init_histo(System &);
void sample_histo(System &, Particle *);
void print_histo(System &);
void init_vacf(System &, Particle *);
void sample_vacf(System &, Particle *);
void print_vacf(System &);



//in draw3D.cpp draw3D.h
extern void draw_stuff();
extern void init_graphics(int arc, char *argv[]);
extern void Idle();
extern void Display();
extern void Init_Graphics(int argc, char *argv[]);
void draw_box(double *,double); 
void update_display(System &,Particle *);
void convert_coordinates(System &,Particle &,double *);
void display_sph1(double *);
void print_test(char *);


//in md.cpp
void md_step(System &, Particle *,Potential &, Thermostat &, Barostat &);

void md_verlet(System &, Particle *,Potential &);
void md_langevin(System &, Particle *,Potential &, Thermostat &);
void md_andersen(System &, Particle *,Potential &, Thermostat &);
void md_nosehooverlangevin_NVT(System &, Particle *,Potential &, Thermostat &);
void md_andersen_NPH(System &, Particle *, Potential &, Barostat &);
void md_andersen_stochastic_NPT(System &, Particle *, Potential &, Barostat &);
void md_andersen_stochastic_nhlthermo_NPT(System &,Particle *,Potential &,Barostat &,Thermostat &);

void langevin_thermo(System &, Particle *, Thermostat &);
void langevin_baro(Barostat &);
void propagate_momenta_half(System &, Particle *);
void propagate_position_half(System &, Particle *);
void propagate_momenta_xi(System &, Particle *, Thermostat &);
void langevin_xi(Thermostat &);
void propagate_position_half_scale(System &, Particle *, double *);
void rescale_position_momenta(System &, Particle *, double *);


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

