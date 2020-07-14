/**************************************************************************************************

  			               Program "MD for LJ"
			         do NVE, velocity verlet MD for LJ 
			LJ: J.Q. Broughton, G.H. Gilmer, JCP 79, 5095 (1983) 
			with CORRECTION!! cf. C.A. Becker et al., PRB 79, 054109 (2009)
				      by Jutta Rogal
				      VERSION 130523

**************************************************************************************************/			

#include "header.h"

//global variables
Statistic block_stats;
Statistic final_stats;
char g_flname[FILENAME_MAX+1]="input.dat";  //default file name
bool g_display=false;
ofstream fpscreen;
FILE *fp[2];

//---------------------------------------------------C++ main program--------------------------------------------------

int main(int argc, char *argv[]){

	//open global screen output file
	fpscreen.open("screen.log");

	//assign input file name if given
    	if (argc > 1) strcpy(g_flname,argv[1]); 

	if((argc >2) && (0==strncmp(argv[2],"--display",9))){
		fpscreen<<"Starting main program and activate graphics output...\n\n";
		g_display = true;
		//try some graphics stuff
		Init_Graphics(argc,argv);
	}
	else{
		fpscreen<<"Starting the main program without graphics output...\n\n";
		g_display = false;
		main_program();
	}

	fpscreen.close();


	return 0;
}

//----------------------------------------------------this is the actual main program----------------------------------

int main_program()
{
	int i,j;
	int k;
	int seed;			//seed for the random number generator
	int64_t count2;			//loop counter for output
	int64_t count_av,sample_av;	//how often to sample averages
	int64_t countplot;		//loop counter for plotting
	int count_vacf;
	System system;			//contain all system variables
	Particle *particle;		//contain particle positions, velocities, forces
	Potential ljpot;		//contains parameters for LJ potential
	Thermostat thermostat;		//contains parameters for different thermostats
	Barostat barostat;		//contains parameters for the differen barostats
	double p[3];			//total momenta
	double max1,max2,dummy;		//to determin max value of something


	Average Av[15];		//get some averages
				// Av[0] = T
	 			// Av[1] = dUtot
				// Av[2] = pressure, P
				// Av[3] = volume
				// Av[4] = pressure without T fluctuation, P'
				// Av[5] = hypervirial
				// Av[6] = kinetic energy
				// Av[7] = potential energy
				// Av[8] = total energy
				// Av[9] = P * Ekin
				// Av[10] = P * Epot
				//  ( Av[11] = Etot + pV )  NOT NEEDED
				//  ( Av[12] = (Etot + pV) * V ) NOT NEEDED
				// Av[13] = Etot + p_fixed * V
				// Av[14] = (Etot + p_fixed * V) * V
	ofstream fplammps;

	fp[0] = fopen("averages.dat","w");
	fp[1] = fopen("ensemble.dat","w");


	fpscreen<<"\n------------------------------------------------------------------------------\n";
        fpscreen<<"Doing MD simulation\n";
        fpscreen<<"Lennard-Jones potential\n";
	fpscreen<<"Version-130523\n";

	//read input file
	readinput(system,ljpot,thermostat,barostat,seed);

	//Print out integrator
	if(system.integrator == 0){
		fpscreen<<"Doing NVE using Velocity Verlet\n";
	}
	else if(system.integrator == 1){
		fpscreen<<"Doing NVT using Langevin thermostat\n";
	}
	else if(system.integrator == 2){
		fpscreen<<"Doing NVT using Andersen thermostat\n";
	}
	else if(system.integrator == 3){
		fpscreen<<"Doing NVT using Nose-Hoover thermostat\n";
	}
	else if(system.integrator == 4){
		fpscreen<<"Doing NVT using Nose-Hoover-Langevin thermostat\n";
	}
	else if(system.integrator == 5){
		fpscreen<<"Doing NPH using Andersen barostat, isotropic\n";
	}
	else if(system.integrator == 6){
		fpscreen<<"Doing NPT using Andersen barostat with Langevin friction (Langevin piston), isotropic\n";
	}
	else if(system.integrator == 7){
		fpscreen<<"Doing NPT using Andersen barostat with Langevin friction (Langevin piston) + NHL thermostat, isotropic\n";
	}








	ran3(seed); //initialize random number generator
        fpscreen<<"\n\n\t First random-number generated: "<<ran3()<<"\n";



	//initialize LJ potential
	init_lj(ljpot);

	//initialise coefficients for different thermostat
	if(system.integrator == 1){		//langevin
		init_langevin(thermostat,system);
	}
	if(system.integrator == 3 || system.integrator == 4){
		init_nhl(thermostat,system);
	}
	//initialise coefficients for Langevin piston barostat
	if(system.integrator == 6){
		init_langevin_piston(barostat,system);
	}
	//langevin piston and nhl thermostat
	if(system.integrator == 7){
		init_nhl(thermostat,system);
		init_langevin_piston(barostat,system);
	}





	//initialize the system
	//read in initial structure from file resp. check box size
	if(system.readstruc<1){
		Particle *dummypart;

		dummypart = new Particle [MAX_PART];
		if(system.readstruc==0){
			read_initlammpsdump(system,dummypart);
		}
		else{
			fpscreen<<"\n\n\tSomething weird, index for readstruc should be 0, but it's "<<system.readstruc<<endl;
			fpscreen<<"\tExiting program...!"<<endl;
		}	

		particle = new Particle [system.nparticles];
		for(i=0;i<system.nparticles;i++){
			for(j=0;j<3;j++){
				particle[i].r[j] = dummypart[i].r[j];
				particle[i].v[j] = dummypart[i].v[j];
			}
			particle[i].mass = dummypart[i].mass;
		}
		delete [] dummypart;
	}
	else{
		//box size, use relative dimensions of the box
		for(i=0;i<3;i++){
			system.box[i] = system.nparticles/system.rho*CUBE(system.box_relative[i])/(system.box_relative[0]*system.box_relative[1]*system.box_relative[2]);
			system.box[i] = pow(system.box[i],1.0/3.0);
			system.box_2[i] = system.box[i]/2.0;
		}
	}
	system.volume = system.box[0] * system.box[1] * system.box[2];

	fpscreen<<"\nSystem dimensions:\n";
	fpscreen<<"length x \t = "<<system.box[0]<<endl;
	fpscreen<<"length y \t = "<<system.box[1]<<endl;
	fpscreen<<"length z \t = "<<system.box[2]<<endl;
	fpscreen<<"volume \t\t = "<<system.volume<<endl;
	fpscreen<<"calculated density \t = "<<(system.nparticles/system.volume)<<endl;



	if(system.readstruc==1){			//setup completely new structure
		//allocate array of particles
		particle = new Particle [system.nparticles];
		init(particle,system,ljpot);		//initialize also particle positions
	}


	//Minimum image condition
	if((system.box[0] < 2.0 * sqrt(ljpot.rmax2)) || (system.box[1] < 2.0 * sqrt(ljpot.rmax2)) || (system.box[2] < 2.0 * sqrt(ljpot.rmax2))){
		fpscreen<<"\n\n\t!WARNING!! Minimum image condition is violated!\n";
		fpscreen<<"Box size x = "<<system.box[0]<<endl;
		fpscreen<<"Box size y = "<<system.box[1]<<endl;
		fpscreen<<"Box size z = "<<system.box[2]<<endl;
		fpscreen<<"2*r_cut = "<<sqrt(ljpot.rmax2)*2.0<<endl;
	}




	//set zoom for display
	if(g_display==true){
		max1 = MAX(system.box[0],system.box[1],dummy);
		max2 = MAX(max1,system.box[2],dummy);
		g_Distance = -max2*4.0;
	}



	//check total momentum
	total_momentum(particle,p,system);
	fpscreen.setf(ios::fixed);
	fpscreen<<"\nTotal momentum px = "<<setprecision(10)<<p[0]<<"\t py = "<<p[1]<<"\t pz = "<<p[2]<<endl;
	fpscreen.unsetf(ios::fixed);
	fpscreen.precision(6);
	fpscreen<<"Total momentum px = "<<p[0]<<"\t py = "<<p[1]<<"\t pz = "<<p[2]<<endl<<endl;


	//initialize counter
        countplot=count2=count_av=count_vacf=0;
        if(system.ncycle2 <= 1000){
		sample_av = 1;
        }
        else{
	      sample_av =  system.ncycle2/1000;
        }
        if(system.outputrate < sample_av){
	       fpscreen<<"\n\n\tSystem output rate must be larger than sample rate!\n";
	       fpscreen<<"\tsystem output rate = "<<system.outputrate<<endl;
	       fpscreen<<"\tsample rate for average = "<<sample_av<<endl;
	       fpscreen<<"\tExiting program...\n\n";
	       exit(1);
        }

       //initialize averages
	init_average(Av);
              

       //initialize forces
       forces(particle,system,ljpot);

       //remap system to be sure everything is in the box
       remap(system,particle);

       //initialize radial distribution function
       init_rdf(system);
       //initialize velocity and volume histograms
       init_histo(system);
       //initialize velocity autocorrelation function
       init_vacf(system,particle);
  
       fprintf(fp[0], "# time \t\t\t Etot/N \t\t <Etot>/N \t\t Ekin/N \t\t <Ekin>/N \t  Epot/N \t\t  <Epot>/N \t  dE/N \t\t  <dE>/N\n");
       fprintf(fp[1], "# time \t\t T_inst  \t T_av  \t\t P_inst  \t P_av  \t\t V_inst  \t V_av  \t\t rho_inst  \t rho_av \n");

       system.E0 = total_energy(particle,system,ljpot);
       fpscreen<<"\nInitial total energy \t = "<<system.E0<<endl<<endl<<endl;
       double K;
       K = kinetic_energy(particle,system);
       fpscreen<<"Temperature before MD starts \t = "<<(2.0*K/((3.0*(double)system.nparticles)-3.0))<<endl;


       fplammps.open("traj.dat");
       init_print_lammps_trajectory(system,particle,fplammps);

       system.time = 0.0; 
      

	//---START MAIN LOOP---

       for(i=0;i<system.ncycle1;i++){
		for(j=0;j<system.ncycle2;j++){
			//do MD step
			md_step(system,particle,ljpot,thermostat,barostat);
			//remap system to the box
			remap(system,particle);

			//average kinetic energy, momenta and state of the system, instantaneous T
			//also update block averages
			average(system,particle);

			//update all other average with some frequency
			if(count_av == sample_av){
				sample_average(system,particle,Av,ljpot);				
				count_av = 0; 
			}
			count_av ++;

			//update display
			if((g_display == true) && (countplot == system.plotrate)){
				update_display(system,particle);
				countplot = 0;
			}
			countplot++;

			if(count2 == system.outputrate){
				print_output(system,particle,Av,ljpot,i,j);
				//print out trajectory
				print_lammps_trajectory(system,particle,fplammps);

				count2 = 0;
			}
			count2++;
			system.time += system.dt;
			count_vacf++;
			if(count_vacf == VACF_SAMPLE){
				sample_vacf(system,particle);
				count_vacf=0;
			}

		}
		//get some statistics
		fpscreen<<"\n\n-----------------------------------------------------------------------------------";
		fpscreen<<"\n\nBLOCK NR. "<<i<<"\n\n";
		terminate_block(system);
		fpscreen<<"\n\n-----------------------------------------------------------------------------------\n";
		//print out running average of rdf
		print_rdf(system);

	}

	//---END MAIN LOOP---

	fpscreen<<"\n---------------------------FINAL STUFF-----------------------------------------\n";
	fpscreen<<endl;
	
	//some final statistics
	fpscreen<<"\n\n-----------------------------------------------------------------------------------";
	fpscreen<<"\n\nFINAL STATISTICS \n\n";
	final_block(system);
	fpscreen<<"\n\n-----------------------------------------------------------------------------------\n";
	
	//print out histograms
	print_histo(system);

	for(k=0;k<2;k++){
		fclose(fp[k]);
	}
	print_rdf(system);
	fclose(system.fprdf);
	delete [] system.rdf;

	//print vacf
	print_vacf(system);

	//close file for printing trajectory
	fplammps.close();




	return 0;
}







//-------------------------------------------terminate block-------------------------------------------------
void terminate_block(System system){

	int i;

	//print out some block statisticA
	block_stats.Ekin.sum /= block_stats.Ekin.n;
	block_stats.T.sum /= block_stats.T.n;
	for(i=0;i<3;i++){
		block_stats.p[i].sum /= block_stats.p[i].n;
	}
	
	fpscreen<<endl;
	fpscreen<<"<Ekin> \t = "<<block_stats.Ekin.sum<<endl;
	fpscreen<<"<T> \t = "<<block_stats.T.sum<<endl;
	fpscreen<<"<px> \t = "<<block_stats.p[0].sum<<endl;
	fpscreen<<"<py> \t = "<<block_stats.p[1].sum<<endl;
	fpscreen<<"<pz> \t = "<<block_stats.p[2].sum<<endl;
	fpscreen<<endl;

	fpscreen<<endl;

	//update final statistic
	update_average(final_stats.Ekin,block_stats.Ekin.sum);
	update_average(final_stats.T,block_stats.T.sum);
	for(i=0;i<3;i++){
		update_average(final_stats.p[i],block_stats.p[i].sum);
	}

	//reset block averages
	block_stats.Ekin.n = block_stats.T.n = 0;
	block_stats.Ekin.in = block_stats.T.in = 0;
	block_stats.Ekin.now = block_stats.T.now = 0.0;
	block_stats.Ekin.sum = block_stats.T.sum = 0.0;
	block_stats.Ekin.sumsq = block_stats.T.sumsq = 0.0;

	for(i=0;i<3;i++){
		block_stats.p[i].n = 0;
		block_stats.p[i].in = 0;
		block_stats.p[i].now = 0.0;
		block_stats.p[i].sum = 0.0;
		block_stats.p[i].sumsq = 0.0;
	}

}



	
//---------------------------------------final statistics----------------------------------------------------

void final_block(System system){

	int i;

	//print out some block statistics
	final_stats.Ekin.sum /= final_stats.Ekin.n;
	final_stats.Ekin.sumsq = sqrt(final_stats.Ekin.sumsq/final_stats.Ekin.n 
					- final_stats.Ekin.sum*final_stats.Ekin.sum);

	final_stats.T.sum /= final_stats.T.n;
	final_stats.T.sumsq = sqrt(final_stats.T.sumsq/final_stats.T.n 
					- final_stats.T.sum*final_stats.T.sum);

	for(i=0;i<3;i++){
		final_stats.p[i].sum /= final_stats.p[i].n;
		final_stats.p[i].sumsq = sqrt(final_stats.p[i].sumsq/final_stats.p[i].n 
						- final_stats.p[i].sum*final_stats.p[i].sum);
	}


	fpscreen<<endl;
	fpscreen<<"<Ekin> \t = "<<final_stats.Ekin.sum
		<<" +/- "<<final_stats.Ekin.sumsq
		<<" ("<<(final_stats.Ekin.sumsq/final_stats.Ekin.sum*100.0)<<"%)"<<endl;

	fpscreen<<"<T> \t = "<<final_stats.T.sum
		<<" +/- "<<final_stats.T.sumsq
		<<" ("<<(final_stats.T.sumsq/final_stats.T.sum*100.0)<<"%)"<<endl;

	fpscreen<<"<px> \t = "<<final_stats.p[0].sum
		<<" +/- "<<final_stats.p[0].sumsq
		<<" ("<<(final_stats.p[0].sumsq/final_stats.p[0].sum*100.0)<<"%)"<<endl;

	fpscreen<<"<py> \t = "<<final_stats.p[1].sum
		<<" +/- "<<final_stats.p[1].sumsq
		<<" ("<<(final_stats.p[1].sumsq/final_stats.p[1].sum*100.0)<<"%)"<<endl;
	fpscreen<<"<pz> \t = "<<final_stats.p[2].sum
		<<" +/- "<<final_stats.p[2].sumsq
		<<" ("<<(final_stats.p[2].sumsq/final_stats.p[2].sum*100.0)<<"%)"<<endl;
	fpscreen<<endl;


	fpscreen<<endl;

}

//-------------------------------------------constructor for class 'Statistic--------------------------------
Statistic::Statistic()
{
	int i;
	
	Ekin.n = T.n = 0;
	Ekin.in = T.in = 0;
	Ekin.now = T.now = 0.0;
	Ekin.sum = T.sum = 0.0;
	Ekin.sumsq = T.sumsq = 0.0;
	for(i=0;i<3;i++){
		p[i].n = 0;
		p[i].in = 0;
		p[i].now = 0.0;
		p[i].sum = 0.0;
		p[i].sumsq = 0.0;
	}

	//assign some names
	Ekin.name = "kinetic energy";
	T.name = "temperature";
	p[0].name = "momentum x";
	p[1].name = "momentum y";
	p[2].name = "momentum z";
}


