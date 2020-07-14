#include "header.h"

//----------------------------------initialize the system------------------------------------

void init(Particle *pparticle,System system,Potential ljpot){
	
	int i,j,k;
	int na;
	double a0,dnn;
	double n,ncell[3];
	double ptot[3];
	double V,K;
	double c;

	//initilize configuration
	//set up an fcc type lattice, considering relative box dimension
	//no. of particles in x,y,z direction
	for(i=0;i<3;i++){
		ncell[i] = (system.nparticles/4.0*CUBE(system.box_relative[i])/(system.box_relative[0]*system.box_relative[1]*system.box_relative[2]));
		ncell[i] = pow(ncell[i],1.0/3.0);
	}
	n = ncell[0]*ncell[1]*ncell[2];
	fpscreen<<endl;
	fpscreen<<"no. of fcc cells \t = "<<n<<endl;
	fpscreen<<"no. of cells in x-direction \t = "<<ncell[0]<<endl;
	fpscreen<<"no. of cells in y-direction \t = "<<ncell[1]<<endl;
	fpscreen<<"no. of cells in z-direction \t = "<<ncell[2]<<endl;

	for(i=0;i<3;i++){
		a0 = system.box[i]/ncell[i];
		fpscreen<<"a0 \t\t\t = "<<a0<<endl;
	}
	dnn = a0/sqrt(2.0);
	fpscreen<<"NN distance from a0 \t = "<<dnn<<endl;
	dnn = sqrt(2.0)/system.rho;
	dnn = pow(dnn,1.0/3.0);
	fpscreen<<"NN distance from rho \t = "<<dnn<<endl;
	fpscreen<<endl;

	na=0;
	for(i=0;i<ncell[2];i++){
		for(j=0;j<ncell[1];j++){
			for(k=0;k<ncell[0];k++){
				pparticle[na].r[0] = k*a0; 
				pparticle[na].r[1] = j*a0; 
				pparticle[na].r[2] = i*a0;
				if((pparticle[na].r[0] < system.box[0]) && (pparticle[na].r[1] < system.box[1]) && (pparticle[na].r[2] < system.box[2])) 
					na++;
				if(na>=system.nparticles) break;

				pparticle[na].r[0] = k*a0+0.5*a0; 
				pparticle[na].r[1] = j*a0+0.5*a0; 
				pparticle[na].r[2] = i*a0;
				if((pparticle[na].r[0] < system.box[0]) && (pparticle[na].r[1] < system.box[1]) && (pparticle[na].r[2] < system.box[2]))
					na++;
				if(na>=system.nparticles) break;

				pparticle[na].r[0] = k*a0+0.5*a0; 
				pparticle[na].r[1] = j*a0; 
				pparticle[na].r[2] = i*a0+0.5*a0;
				if((pparticle[na].r[0] < system.box[0]) && (pparticle[na].r[1] < system.box[1]) && (pparticle[na].r[2] < system.box[2]))
					na++;
				if(na>=system.nparticles) break;
				
				pparticle[na].r[0] = k*a0; 
				pparticle[na].r[1] = j*a0+0.5*a0; 
				pparticle[na].r[2] = i*a0+0.5*a0;
				if((pparticle[na].r[0] < system.box[0]) && (pparticle[na].r[1] < system.box[1]) && (pparticle[na].r[2] < system.box[2]))
					na++;
				if(na>=system.nparticles) break;
				
			}
			if(na>=system.nparticles) break;
		}
		if(na>=system.nparticles) break;
	}

	fpscreen<<"No of atoms distributed on the lattice \t = "<<na<<endl;
	fpscreen<<"No of atoms originally put in          \t = "<<system.nparticles<<endl;
	fpscreen<<endl;
	//print out particle postitions
//	fpscreen<<endl;
//	for(i=0;i<system.nparticles;i++){
//		fpscreen<<i<<fixed<< setprecision(5)<<"\t"<<pparticle[i].r[0]<<"\t\t"<<pparticle[i].r[1]<<"\t\t"<<pparticle[i].r[2]<<endl;
//	}


	//initialize mass of particles
	for(i=0;i<system.nparticles;i++){
		pparticle[i].mass = 1.0;
	}
	fpscreen<<"Particle mass hardcoded to 1.0!\n\n";

	//initialize momenta
	for(i=0;i<3;i++){
		ptot[i] = 0.0;
	}

	//momenta from Gauss distribution for inverse temperature beta
	for(i=0;i<system.nparticles;i++){
		for(j=0;j<3;j++){
			//pparticle[i].v[j] = RandomVelocity(system.beta/2.0,pparticle[i].mass); 
			pparticle[i].v[j] = RandomVelocity(system.beta,pparticle[i].mass); // to see the effect of equipartition initialize with T
			ptot[j] += pparticle[i].v[j]*pparticle[i].mass;
		}
	}
	K = kinetic_energy(pparticle,system);
	fpscreen<<"Temperature before ptot correction \t = "<<(2.0*K/((3.0*(double)system.nparticles)-3.0))<<endl;

	//reset the momenta so that ptot_x and ptot_y = 0.0
	//excess momentum per particle
	for(i=0;i<3;i++){
		ptot[i] /= ((double) system.nparticles);
	}

	//substract excess momentum for each particle
	for(i=0;i<system.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] -= ptot[j]/pparticle[i].mass;
		}
	}

	//get potential and kinetic energy
	V = potential_energy(pparticle,system,ljpot);
	K = kinetic_energy(pparticle,system);
	fpscreen<<"Temperature after ptot correction \t = "<<(2.0*K/((3.0*(double)system.nparticles)-3.0))<<endl;

	//rescaling value to set total energy to input value
	c = (system.tot_energy - V)/K;

	//check if potential energy of the system is smaller than total energy
	if(c < 0.0){
		fpscreen<<"\n\n\t!ERROR!! Potential energy of the system is larger than total energy!!\n";
		fpscreen<<"E_pot = "<<V<<endl;
		fpscreen<<"E_tot = "<<system.tot_energy<<endl;
		exit(0);
	}

	//rescale velocities so that E_tot = K + V
	c = sqrt(c);
	fpscreen<<"\n\tNo rescaling, kinetic energy from setting beta\n\n";
	c=1.0;
	for(i=0;i<system.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] *= c;
		}
	}

	
	

}


//--------------------------------initialize energies and velocities for read in structure---------------------------


void init_ene(Particle *pparticle,System &psystem,Potential &pljpot){

	int i,j;
	double ptot[3];
	double V,K;
	double c;


	//initialize mass of particles
	for(i=0;i<psystem.nparticles;i++){
		pparticle[i].mass = 1.0;
	}
	fpscreen<<"Particle mass hardcoded to 1.0!\n\n";

	//initialize momenta
	for(i=0;i<3;i++){
		ptot[i] = 0.0;
	}

	//momenta from Gauss distribution for inverse temperature beta
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			//pparticle[i].v[j] = RandomVelocity(psystem.beta,pparticle[i].mass); 
			pparticle[i].v[j] = RandomVelocity(psystem.beta,pparticle[i].mass);  // to see the effect of equipartition initialize with T
			ptot[j] += pparticle[i].v[j]*pparticle[i].mass;
		}
	}
	K = kinetic_energy(pparticle,psystem);
	fpscreen<<"Temperature before ptot correction \t = "<<(2.0*K/((3.0*(double)psystem.nparticles)-3.0))<<endl;

	//reset the momenta so that ptot_x and ptot_y = 0.0
	//excess momentum per particle
	for(i=0;i<3;i++){
		ptot[i] /= ((double) psystem.nparticles);
	}

	//substract excess momentum for each particle
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] -= ptot[j]/pparticle[i].mass;
		}
	}

	//get potential and kinetic energy
	V = potential_energy(pparticle,psystem,pljpot);
	K = kinetic_energy(pparticle,psystem);
	fpscreen<<"Temperature after ptot correction \t = "<<(2.0*K/((3.0*(double)psystem.nparticles)-3.0))<<endl;

	//rescaling value to set total energy to input value
	c = (psystem.tot_energy - V)/K;

	//check if potential energy of the system is smaller than total energy
	if(c < 0.0){
		fpscreen<<"\n\n\t!ERROR!! Potential energy of the system is larger than total energy!!\n";
		fpscreen<<"E_pot = "<<V<<endl;
		fpscreen<<"E_tot = "<<psystem.tot_energy<<endl;

		fpscreen<<"Doing initial steepest descent relaxation\n";
		relax_steep(pparticle,psystem,pljpot);
		//exit(0);
	}

	//rescale velocities so that E_tot = K + V
	c = sqrt(c);
	fpscreen<<"\n\tNo rescaling, kinetic energy from setting beta\n\n";
	c=1.0;
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] *= c;
		}
	}


}

//---------------------------do some initial steepest descent steps----------------------------------------------------

void relax_steep(Particle *pparticle,System &psystem,Potential &pljpot){

	double scale;	//scale the steepest descent
	double dummy,tmp;
	double Epot, Epot_old;
	int i,j,k;
	Particle *old_part;
//	double old_f[3];
	double m[3];

	dummy = 1.0/psystem.rho;
	scale = pow(dummy,(1.0/3.0));
	scale = 0.02*scale;
	
	old_part = new Particle [psystem.nparticles];

	Epot_old = 0.0;
//	for(k=0;k<3;k++){
//		old_r[k] = 0.0;
//	}

	for(i=0;i<100;i++){
		if(i==0){
			Epot = potential_energy(pparticle,psystem,pljpot);
			Epot_old = Epot;
			fpscreen<<"\nRelaxing initial configuration\n";
			fpscreen<<"Initial Energy = "<<Epot<<endl;
			forces(pparticle,psystem,pljpot);
		}

		//calculate maximum downhill gradient
		tmp = 0.0;
		for(j=0;j<psystem.nparticles;j++){
			for(k=0;k<3;k++){
				old_part[j].r[k] = pparticle[j].r[k];
//				old_f[k] = pparticle[j].f[k];
				m[k] = fabs(pparticle[j].f[k]);

				if(m[k] > tmp){
					tmp = m[k];
				}
			}
		}

		fpscreen<<"i = "<<i<<"\tmax force = "<<tmp<<endl;
		tmp = scale/tmp;
		//fpscreen<<"i = "<<i<<"\ttmp = "<<tmp<<endl;

		//calculate improved positons
		for(j=0;j<psystem.nparticles;j++){
			for(k=0;k<3;k++){
				pparticle[j].r[k] += tmp*pparticle[j].f[k];
			}
		}
		
       		remap(psystem,pparticle);
		Epot = potential_energy(pparticle,psystem,pljpot);
		
		if(Epot < Epot_old){
			Epot_old = Epot;
			//if(i>30){
				scale *= 1.2;
			//}
			for(k=0;k<3;k++){
				if(scale > psystem.box_2[k]){
					scale = psystem.box_2[k];
				}
			}
		}
		else{
			for(j=0;j<psystem.nparticles;j++){
				for(k=0;k<3;k++){
					pparticle[j].r[k] = old_part[j].r[k];
				}
			}
			scale *= 0.1;
		}
		forces(pparticle,psystem,pljpot);
		fpscreen<<"i = "<<i<<"\tEpot = "<<Epot<<endl;
		//fpscreen.setf(ios::fixed);
		//fpscreen.precision(6);
		//for(j=0;j<psystem.nparticles;j++){
		//	fpscreen<<"\t"<<pparticle[j].r[0]<<"\t"<<pparticle[j].r[1]<<"\t"<<pparticle[j].r[2];
		//	fpscreen<<"\t"<<pparticle[j].f[0]<<"\t"<<pparticle[j].f[1]<<"\t"<<pparticle[j].f[2]<<endl;
		//}
		//fpscreen<<endl;

	}
	
	Epot = potential_energy(pparticle,psystem,pljpot);
	fpscreen<<"Final Energy = "<<Epot<<endl;
			


}



//--------------------------initialise Lennard-Jones potential---------------------------------------------------------

void init_lj(Potential &pljpot)
{

	pljpot.c1 = 0.016132*pljpot.epsilon;
	pljpot.c2 = 3136.6*pljpot.epsilon;
	pljpot.c3 = -68.069*pljpot.epsilon;
	pljpot.c4 = -0.083312*pljpot.epsilon;
	pljpot.c5 = 0.74689*pljpot.epsilon;

	pljpot.rmin2 = SQR(2.3*pljpot.sigma);
	pljpot.rmax2 = SQR(2.5*pljpot.sigma);

	fpscreen<<"LJ potential coefficients\n";
	fpscreen<<"C1 \t\t = "<<pljpot.c1<<endl;
	fpscreen<<"C2 \t\t = "<<pljpot.c2<<endl;
	fpscreen<<"C3 \t\t = "<<pljpot.c3<<endl;
	fpscreen<<"C4 \t\t = "<<pljpot.c4<<endl;
	fpscreen<<"C5 \t\t = "<<pljpot.c5<<endl;
	fpscreen<<"LJ potential lower and upper cutoff:\n";
	fpscreen<<"rmin \t\t = "<<sqrt(pljpot.rmin2)<<endl;
	fpscreen<<"rmax \t\t = "<<sqrt(pljpot.rmax2)<<endl;
	fpscreen<<endl;
}


//--------------------------------------initialise Langevin coefficients-----------------------------------------------

void init_langevin(Thermostat &pthermostat, System &psystem)
{
	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//coefficents are
	// c1 = exp(-gamma*dt/2.0)
	// c2 = sqrt((1-c1**2)*m*kT) = sqrt(m) c2'
	// to propagate velocities instead of momenta => devide by m
	// !! for c2' sqrt(m)/m = 1/sqrt(m) has to be added individually for each atom during the MD step!!!

	pthermostat.lc1 = exp(-pthermostat.lgamma*psystem.dt/2.0);
	pthermostat.lc2 = sqrt((1.0-SQR(pthermostat.lc1))/psystem.beta);        //this is c2'

	fpscreen<<"Langevin thermostat coefficients\n";
	fpscreen<<"gamma \t\t = "<<pthermostat.lgamma<<endl;
	fpscreen<<"c1    \t\t = "<<pthermostat.lc1<<endl;
	fpscreen<<"c2    \t\t = "<<pthermostat.lc2<<endl;
	fpscreen<<endl;

	//additional energy contribution due to langevin substracted from total energy!
	psystem.dElangevin = 0.0;

}

//--------------------------------------initialise Langevin coefficients for barostat piston---------------------------

void init_langevin_piston(Barostat &pbaro, System system)
{
	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//coefficents are
	// c1 = exp(-gamma*dt/2.0)
	// c2 = sqrt((1-c1**2)*m*kT) = sqrt(m) c2'
	// to propagate velocities instead of momenta => devide by m
	// !! for c2' sqrt(m)/m = 1/sqrt(m) has to be added individually for each atom during the MD step!!!

	pbaro.lc1 = exp(-pbaro.lgamma*system.dt/2.0);
	pbaro.lc2 = sqrt((1.0-SQR(pbaro.lc1))/system.beta);        //this is c2'

	//only one mass for the piston !
	pbaro.lc2 = pbaro.lc2/sqrt(pbaro.pmass);

	fpscreen<<"Langevin piston coefficients\n";
	fpscreen<<"gamma piston\t\t = "<<pbaro.lgamma<<endl;
	fpscreen<<"c1 piston    \t\t = "<<pbaro.lc1<<endl;
	fpscreen<<"c2 piston   \t\t = "<<pbaro.lc2<<endl;
	fpscreen<<endl;
}

//--------------------------------------initialise Nose-Hoover-Langevin------------------------------------------------

void init_nhl(Thermostat &pthermostat, System system)
{
	//Ref:  Ben Leimkuhler (A gentle thermostat)
	//coefficents are
	// c1 = exp(-gamma*dt/2.0)
	// c2 = sqrt((1-c1**2)*kT*2/mu)

	if(system.integrator == 3){
		pthermostat.nhlgamma = 0.0;
	}

	pthermostat.nhlc1 = exp(-pthermostat.nhlgamma*system.dt/2.0);
	pthermostat.nhlc2 = sqrt((1.0-SQR(pthermostat.nhlc1))*2.0/system.beta/pthermostat.nhlmu);        

	//initialise xi
	pthermostat.nhlxi = RandomVelocity(system.beta,pthermostat.nhlmu);   //maybe use mu/2 to initialize (?)

	fpscreen<<"Nose-Hoover-Langevin thermostat coefficients\n";
	fpscreen<<"gamma \t\t = "<<pthermostat.nhlgamma<<endl;
	fpscreen<<"c1    \t\t = "<<pthermostat.nhlc1<<endl;
	fpscreen<<"c2    \t\t = "<<pthermostat.nhlc2<<endl;
	fpscreen<<"xi    \t\t = "<<pthermostat.nhlmu<<endl;
	fpscreen<<endl;
}

