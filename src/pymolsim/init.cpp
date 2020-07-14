#include "header.h"


//--------------------------------initialize energies and velocities for read in structure---------------------------


void System::init(){

	int i,j;
	double V,K;
	double c;

	//initialize mass of particles
	for(i=0;i<nparticles;i++){
		particles[i].mass = 1.0;
	}

	//initialize momenta
	for(i=0;i<3;i++){
		p[i] = 0.0;
	}

	//momenta from Gauss distribution for inverse temperature beta
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			//pparticle[i].v[j] = RandomVelocity(psystem.beta,pparticle[i].mass); 
			particles[i].v[j] = RandomVelocity(beta, particles[i].mass);  // to see the effect of equipartition initialize with T
			p[j] += particles[i].v[j]*particles[i].mass;
		}
	}

	kinetic_energy();

	//reset the momenta so that ptot_x and ptot_y = 0.0
	//excess momentum per particle
	for(i=0;i<3;i++){
		p[i] /= ((double) nparticles);
	}

	//substract excess momentum for each particle
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] -= p[j]/particles[i].mass;
		}
	}

	//get potential and kinetic energy
	potential_energy();
	kinetic_energy();

	//rescaling value to set total energy to input value
	c = (tot_energy - V)/K;

	//rescale velocities so that E_tot = K + V
	c = sqrt(c);
	c=1.0;
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] *= c;
		}
	}


}


vector<Particle> System::gparticles(){
	return particles;
}
void System::sparticles(vector<Particle> sparticles){
	particles = sparticles;
}
Potential System::gpotential(){
	return potential;
}
void System::spotential(Potential spot){
	potential = spot;
}

Thermostat System::gthermostat(){
	return thermostat;
}
void System::sthermostat(Thermostat stherm){
	thermostat = stherm;
}

Barostat System::gbarostat(){
	return barostat;
}
void System::sbarostat(Barostat sbar){
	barostat = sbar;
}

vector<Average> System::gaverage(){
	return pAv;
}
void System::saverage(vector<Average> svg){
	pAv = svg;
}
//--------------------------initialise Lennard-Jones potential---------------------------------------------------------

void Potential::init()
{

	c1 = 0.016132*epsilon;
	c2 = 3136.6*epsilon;
	c3 = -68.069*epsilon;
	c4 = -0.083312*epsilon;
	c5 = 0.74689*epsilon;

	rmin2 = SQR(2.3*sigma);
	rmax2 = SQR(2.5*sigma);

}


//--------------------------------------initialise Langevin coefficients-----------------------------------------------

void Thermostat::init(double dt, double beta)
{
	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//coefficents are
	// c1 = exp(-gamma*dt/2.0)
	// c2 = sqrt((1-c1**2)*m*kT) = sqrt(m) c2'
	// to propagate velocities instead of momenta => devide by m
	// !! for c2' sqrt(m)/m = 1/sqrt(m) has to be added individually for each atom during the MD step!!!

	lc1 = exp(-lgamma*dt/2.0);
	lc2 = sqrt((1.0-SQR(lc1))/beta);        //this is c2'

	//additional energy contribution due to langevin substracted from total energy!
	dElangevin = 0.0;

}

//--------------------------------------initialise Langevin coefficients for barostat piston---------------------------

void Barostat::init(double dt, double beta)
{
	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//coefficents are
	// c1 = exp(-gamma*dt/2.0)
	// c2 = sqrt((1-c1**2)*m*kT) = sqrt(m) c2'
	// to propagate velocities instead of momenta => devide by m
	// !! for c2' sqrt(m)/m = 1/sqrt(m) has to be added individually for each atom during the MD step!!!

	lc1 = exp(-lgamma*dt/2.0);
	lc2 = sqrt((1.0-SQR(lc1))/beta);        //this is c2'

	//only one mass for the piston !
	lc2 = lc2/sqrt(pmass);
}

//--------------------------------------initialise Nose-Hoover-Langevin------------------------------------------------

void Thermostat::init_nhl(double dt, double beta, int integrator)
{
	//Ref:  Ben Leimkuhler (A gentle thermostat)
	//coefficents are
	// c1 = exp(-gamma*dt/2.0)
	// c2 = sqrt((1-c1**2)*kT*2/mu)

	if(integrator == 3){
		nhlgamma = 0.0;
	}

	nhlc1 = exp(-nhlgamma*dt/2.0);
	nhlc2 = sqrt((1.0-SQR(nhlc1))*2.0/beta/nhlmu);        

	//initialise xi
	nhlxi = RandomVelocity(beta,nhlmu);   //maybe use mu/2 to initialize (?)
}

