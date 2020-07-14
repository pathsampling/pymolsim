#include "header.h"

//-------------------------------calculate distance between two particles----------------------

vector<double> System::image_distance(int part_1, int part_2){

	int i;
	vector<double> dr{0., 0., 0.};
	//get x distance, x2-x1 = x21
	for(i=0;i<3;i++){
		dr[i] = particles[part_2].r[i] - particles[part_1].r[i];
		if(dr[i] > box_2[i]){
			dr[i] -= box[i];
		}
		else if (dr[i] < -box_2[i]){
			dr[i] += box[i];
		}
	}
	return dr;
	
}


//------------------------------------------remap the system into the box-------------------------------

void System::remap(){

	int i,j;
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			if(particles[i].r[j] > box[j])
				particles[i].r[j] -= box[j];
			if(particles[i].r[j] < 0.0)
				particles[i].r[j] += box[j];
		}
	}
}


//------------------------------------------do some averaging-------------------------------------------

void System::average(){
	
	int i;
	double T;


	//average kinetic energy
	kinetic_energy();
	update_average(block_stats.Ekin, ke);
	//average temperature
	T = 2.0*ke/((3.0*(double)nparticles)-3.0);
	update_average(block_stats.T, T);
	Tinst = T;
	
	//average total momenta
	total_momentum();
	for(i=0;i<3;i++){
		update_average(block_stats.p[i],p[i]);
	}

}

//--------------------------------------some more averaging---------------------------------------------

void System::sample_average(){
				
	
	double Ekin,Epot,E,dE;
	
	//running average for the temperature
	update_average(pAv[0], Tinst);
	
	//calculate energies and energy drift
	total_energy();
	dE = fabs((tot_energy - E0)/(double)nparticles);
	
	update_average(pAv[1],dE);
	update_average(pAv[6],ke);
	update_average(pAv[7],pe);
	update_average(pAv[8],tot_energy);

	//calculate kinetic part of the pressure (NKT/V) and do average
	//make sure Tinst has the proper value here ! calculated in 'average'
	Pinst = virial + (nparticles*Tinst/volume);
	update_average(pAv[2], Pinst);
	update_average(pAv[3], volume);
	update_average(pAv[9],(Pinst * ke));
	update_average(pAv[10],(Pinst * pe));
	update_average(pAv[13],(tot_energy + pressure*volume));
	update_average(pAv[14],((tot_energy + pressure*volume) * volume));
	//pressure without T fluctuations
	Pinst_1 = virial + (nparticles/volume/beta);
	update_average(pAv[4],Pinst_1);
	//update hypervirial
	update_average(pAv[5],hypervirial);

	//TODO : Update theese functions
	//update radial distribution function
	//sample_rdf(psystem,pparticle);
	//update histograms
	//sample_histo(psystem,pparticle);

}

//-----------------------------------initialize average counter-----------------------------------------
void System::init_average(){
	int i;

	for(i=0;i<15;i++){
	       pAv[i].n = 0;
	       pAv[i].in = 0;
	       pAv[i].now = 0.0;
	       pAv[i].sum = 0.0;
	       pAv[i].sumsq = 0.0;
       }

       pAv[0].name = "temperature";
       pAv[1].name = "energy drift";
       pAv[2].name = "pressure";
       pAv[3].name = "volume";
       pAv[4].name = "pressure with system T";
       pAv[5].name = "hypervirial";
       pAv[6].name = "kinetic energy";
       pAv[7].name = "potential energy";
       pAv[8].name = "total energy";
       pAv[9].name = "pressure * kinetic energy";
       pAv[10].name = "pressure * potential energy";
       pAv[11].name = "total energy + pV";                  //NOT NEEDED
       pAv[12].name = "(total energy + pV) * V";           // NOT NEEDED
       pAv[13].name = "total energy + p_fix * V";
       pAv[14].name = "(total energy + p_fix * V) * V";

}


