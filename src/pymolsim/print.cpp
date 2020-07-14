#include "header.h"


//--------------------------print slice in lammps format------------------------------------------
void print_lammps_dump(System &psystem, Particle *pparticle, ofstream &pfplammps){


	int i,j;
	pfplammps<<"ITEM: TIMESTEP"<<endl;
	pfplammps<<"0"<<endl;
	pfplammps<<"ITEM: NUMBER OF ATOMS"<<endl;
	pfplammps<<psystem.nparticles<<endl;
	pfplammps<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
	pfplammps<<"0 "<<psystem.box[0]<<endl;
	pfplammps<<"0 "<<psystem.box[1]<<endl;
	pfplammps<<"0 "<<psystem.box[2]<<endl;
	pfplammps<<"ITEM: ATOMS id type mass x y z vx vy vz"<<endl;
//	pfplammps.setf(ios::scientific);
//	pfplammps.precision(6);

	for(i=0;i<psystem.nparticles;i++){
		pfplammps<<i+1<<" 1 1.000 ";
		for(j=0;j<3;j++){
			pfplammps<<pparticle[i].r[j]<<" ";
		}
		for(j=0;j<3;j++){
			pfplammps<<pparticle[i].v[j]<<" ";
		}
		pfplammps<<endl;
	}

	pfplammps.flush();

}

//-------------------------init printing trajectory in lammps dump format--------------------------
void init_print_lammps_trajectory(System &psystem, Particle *pparticle, ofstream &pfplammps){

	int i,j;


	pfplammps<<"ITEM: TIMESTEP"<<endl;
	pfplammps<<"0"<<endl;
	pfplammps<<"ITEM: NUMBER OF ATOMS"<<endl;
	pfplammps<<psystem.nparticles<<endl;
	pfplammps<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
	pfplammps<<"0 "<<psystem.box[0]<<endl;
	pfplammps<<"0 "<<psystem.box[1]<<endl;
	pfplammps<<"0 "<<psystem.box[2]<<endl;
	pfplammps<<"ITEM: ATOMS id type mass x y z vx vy vz"<<endl;

	for(i=0;i<psystem.nparticles;i++){
		pfplammps<<i+1<<" 1 1.000 ";
		for(j=0;j<3;j++){
			pfplammps<<pparticle[i].r[j]<<" ";
		}
		for(j=0;j<3;j++){
			pfplammps<<pparticle[i].v[j]<<" ";
		}
		pfplammps<<endl;
	}

	pfplammps.flush();

}

//----------------------------- printing trajectory in lammps dump format--------------------------
void print_lammps_trajectory(System &psystem, Particle *pparticle, ofstream &pfplammps){

	int i,j;
	
	pfplammps<<"ITEM: TIMESTEP"<<endl;
	pfplammps<<psystem.time<<endl;
	pfplammps<<"ITEM: NUMBER OF ATOMS"<<endl;
	pfplammps<<psystem.nparticles<<endl;
	pfplammps<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
	pfplammps<<"0 "<<psystem.box[0]<<endl;
	pfplammps<<"0 "<<psystem.box[1]<<endl;
	pfplammps<<"0 "<<psystem.box[2]<<endl;
	pfplammps<<"ITEM: ATOMS id type mass x y z vx vy vz"<<endl;

	for(i=0;i<psystem.nparticles;i++){
		pfplammps<<i+1<<" 1 1.000 ";
		for(j=0;j<3;j++){
			pfplammps<<pparticle[i].r[j]<<" ";
		}
		for(j=0;j<3;j++){
			pfplammps<<pparticle[i].v[j]<<" ";
		}
		pfplammps<<endl;
	}

	pfplammps.flush();


}


//--------------------------------print some output-------------------------------------------------
void print_output(System &psystem, Particle *pparticle, Average *pAv, Potential &pljpot,int i,int j){

	//determine total energy and total momenta
	fpscreen.setf(ios::fixed);
	fpscreen<<"n = "<<(i*psystem.ncycle2+j)<<"\t Etot = "<<pAv[8].now
			<<"\t T_inst = "<<psystem.Tinst<<"\t T_av = "<<(pAv[0].sum/pAv[0].n)
			<<"\t dE_av = "<<(pAv[1].sum/pAv[1].n)
			<<"\t P_inst = "<<psystem.Pinst<<"\t P_av = "<<(pAv[2].sum/pAv[2].n)
			<<"\t V_inst = "<<psystem.volume<<"\t V_av = "<<(pAv[3].sum/pAv[3].n)
			<<endl;
	fflush(stdout);
	
	fprintf(fp[0], "%14.6f  %20.6Lf  %20.6Lf  %20.6Lf  %20.6Lf  %20.6Lf  %20.6Lf %14.4Le  %14.4Le\n",
			psystem.time,
			pAv[8].now/psystem.nparticles,(pAv[8].sum/pAv[8].n)/psystem.nparticles,
			pAv[6].now/psystem.nparticles,(pAv[6].sum/pAv[6].n)/psystem.nparticles,
			pAv[7].now/psystem.nparticles,(pAv[7].sum/pAv[7].n)/psystem.nparticles,
			pAv[1].now,(pAv[1].sum/pAv[1].n));
	fflush(fp[0]);

	fprintf(fp[1], "%14.6f  %14.6f  %14.6Lf  %14.6f  %14.6Lf  %14.6f  %14.6Lf  %14.6f  %14.6Lf\n",
			psystem.time,
			psystem.Tinst,(pAv[0].sum/pAv[0].n),
			psystem.Pinst,(pAv[2].sum/pAv[2].n),
			psystem.volume,(pAv[3].sum/pAv[3].n),
			(psystem.nparticles/psystem.volume),(psystem.nparticles/(pAv[3].sum/pAv[3].n)));
	fflush(fp[1]);

}



























