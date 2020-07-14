#include "header.h"

//-------------------------------calculate distance between two particles----------------------

void image_distance(Particle &part_1,Particle &part_2,double *dr,System &psystem){

	int i;
	//get x distance, x2-x1 = x21
	for(i=0;i<3;i++){
		dr[i] = part_2.r[i] - part_1.r[i];
		if(dr[i] > psystem.box_2[i]){
			dr[i] -= psystem.box[i];
		}
		else if (dr[i] < -psystem.box_2[i]){
			dr[i] += psystem.box[i];
		}
	}
	

}


//------------------------------------------remap the system into the box-------------------------------

void remap(System &psystem, Particle *pparticle){

	int i,j;

	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			if(pparticle[i].r[j] > psystem.box[j])
				pparticle[i].r[j] -= psystem.box[j];
			if(pparticle[i].r[j] < 0.0)
				pparticle[i].r[j] += psystem.box[j];
		}
	}
}


//------------------------------------------do some averaging-------------------------------------------

void average(System &psystem, Particle *pparticle){
	
	int i;
	double K;
	double p[3];
	double T;


	//average kinetic energy
	K = kinetic_energy(pparticle,psystem);
	update_average(block_stats.Ekin,K);
	//average temperature
	T = 2.0*K/((3.0*(double)psystem.nparticles)-3.0);
	update_average(block_stats.T,T);
	psystem.Tinst = T;
	
	//average total momenta
	total_momentum(pparticle,p,psystem);
	for(i=0;i<3;i++){
		update_average(block_stats.p[i],p[i]);
	}

}

//--------------------------------------some more averaging---------------------------------------------

void sample_average(System &psystem, Particle *pparticle, Average *pAv, Potential &pljpot){
				
	
	double Ekin,Epot,E,dE;
	
	//running average for the temperature
	update_average(pAv[0],psystem.Tinst);
	//calculate energies and energy drift
	Ekin = kinetic_energy(pparticle,psystem);
	Epot = potential_energy(pparticle,psystem,pljpot);
	E = Ekin + Epot;
	dE = fabs((E - psystem.E0)/psystem.nparticles);
	update_average(pAv[1],dE);
	update_average(pAv[6],Ekin);
	update_average(pAv[7],Epot);
	update_average(pAv[8],E);
	//calculate kinetic part of the pressure (NKT/V) and do average
	//make sure Tinst has the proper value here ! calculated in 'average'
	psystem.Pinst = psystem.virial + (psystem.nparticles*psystem.Tinst/psystem.volume);
	update_average(pAv[2],psystem.Pinst);
	update_average(pAv[3],psystem.volume);
	update_average(pAv[9],(psystem.Pinst * Ekin));
	update_average(pAv[10],(psystem.Pinst * Epot));
	update_average(pAv[13],(E + psystem.pressure*psystem.volume));
	update_average(pAv[14],((E + psystem.pressure*psystem.volume) * psystem.volume));
	//pressure without T fluctuations
	psystem.Pinst_1 = psystem.virial + (psystem.nparticles/psystem.volume/psystem.beta);
	update_average(pAv[4],psystem.Pinst_1);
	//update hypervirial
	update_average(pAv[5],psystem.hypervirial);

	//update radial distribution function
	sample_rdf(psystem,pparticle);
	//update histograms
	sample_histo(psystem,pparticle);

}

//-----------------------------------initialize average counter-----------------------------------------
void init_average(Average *pAv){
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


//------------------------------------initialize radial distribution function---------------------------

void init_rdf(System &psystem){

	double min,dummy;
	int i;
	//rdf can be sampled up to distances of 1/2 box_size
	//and I think it makes only sense if you take the shortest direction
	min = MIN(psystem.box_2[0],psystem.box_2[1],dummy);
	psystem.rmax_rdf = MIN(min,psystem.box_2[2],dummy);

	fpscreen<<"rdf max radius = "<<psystem.rmax_rdf<<endl;
	fpscreen<<endl;

	//determine number of bins necessary 
	psystem.rdfBIN = (int) ceil(psystem.rmax_rdf * psystem.rdfRES) + 1;
	//allocate array for rdf
	psystem.rdf = new int64_t [psystem.rdfBIN];
	//initialize array
	for(i=0;i<psystem.rdfBIN;i++){
		psystem.rdf[i] = 0;
	}
	psystem.rdfcount = 0;

	psystem.fprdf = fopen("rdf.dat","wa");
	fprintf(psystem.fprdf,"# r \t \t g(r)\n");

}

//-------------------------------------sample radial distribution function-----------------------------

void sample_rdf(System &psystem, Particle *pparticle){

	int i,j,m;
	double r;
	double dr[3];

	//update counter
	psystem.rdfcount++;

	//loop over all pairs
	for(i=0;i<psystem.nparticles-1;i++){
		for(j=i+1;j<psystem.nparticles;j++){
			//get distance between particles considering PBC
			// xji in this case
			image_distance(pparticle[i],pparticle[j],dr,psystem);
			r = sqrt(SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]));
			//update g(r)
			if(r < psystem.rmax_rdf){
				m = (int) rint(r*psystem.rdfRES);
				if(m>=psystem.rdfBIN){
                                        fpscreen<<"\nrdf sampling exceeds number of bins...exiting program\n\n";
                                        exit(1);
                                }
				psystem.rdf[m] += 2;
			}
		}
	}

}

//--------------------------------print radial distribution function------------------------------------

void print_rdf(System &psystem){

	int i;
	long double norm,gr;
	//loop over all bins
	for(i=0;i<psystem.rdfBIN-1;i++){
		//calculate number of particles in the bin volume of an ideal gas at this density rho
		norm = (4.0/3.0*M_PI) * (psystem.nparticles/psystem.volume) * (CUBE(i+1) - CUBE(i)) / CUBE(psystem.rdfRES);
		gr = psystem.rdf[i]/(psystem.rdfcount*psystem.nparticles*norm);
		fprintf(psystem.fprdf,"%6.4e  %6.4Le\n",(double)i/psystem.rdfRES,gr);
	}
	fprintf(psystem.fprdf,"\n");
	fflush(psystem.fprdf);
}




//finish rdf
//free memory, close fp
//currently done at the end of the main program

//------------------------------------initialize histograms---------------------------------------------

void init_histo(System &psystem){

	int i,j;
	//initialize array for velocity histogram
	for(i=0;i<4;i++){
		for(j=0;j<MBIN_v;j++){
			psystem.histo_v[i][j] = 0;
		}
		psystem.vcount[i] = 0;
	}

	//initialize array for volume histogram
	for(i=0;i<MBIN_Vol;i++){
		psystem.histo_Vol[i] = 0;
	}
	psystem.Volcount = 0;
	psystem.histo_Vol0 = psystem.volume;


}

//----------------------------------sample histograms----------------------------------------------------

void sample_histo(System &psystem, Particle *pparticle){

	int i,j;
	int n;
	double v;


	//sample velocities
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			n = (int) rint(pparticle[i].v[j] * RESHISTO_v + MBIN_v/2.0);
			if((n>=0) && (n<MBIN_v)){
				psystem.histo_v[j][n]++;
				psystem.vcount[j]++;
			}
		}
		v = sqrt(SQR(pparticle[i].v[0]) + SQR(pparticle[i].v[1]) + SQR(pparticle[i].v[2]));
		n = (int) rint(v*RESHISTO_v);
		if((n>=0) && (n<MBIN_v)){
			psystem.histo_v[3][n]++;
			psystem.vcount[3]++;
		}
	}

	//sample volume
	n = (int) rint((psystem.volume - psystem.histo_Vol0) * RESHISTO_Vol + MBIN_Vol/2.0);
	if((n>=0) && (n<MBIN_Vol)){
		psystem.histo_Vol[n]++;
		psystem.Volcount++;
	}


}

//---------------------------------------print histograms--------------------------------------------------

void print_histo(System &psystem){

	int i,j;
	FILE *fphisto[2];
	fphisto[0] = fopen("velocity_dist.dat","w");
	fphisto[1] = fopen("volume_dist.dat","w");

	fprintf(fphisto[0],"#vx/vy/vz \t P(vx) \t P(vy) \t P(vz) \t v \t P(v)\n");
	for(i=0;i<MBIN_v;i++){
		fprintf(fphisto[0],"%14.6f",(double)(i-MBIN_v/2.0)/RESHISTO_v);
		for(j=0;j<3;j++){
			fprintf(fphisto[0],"%14.6Lf",(long double)psystem.histo_v[j][i]/psystem.vcount[j]*RESHISTO_v);
		}
		fprintf(fphisto[0],"%14.6f",(double)i/RESHISTO_v);
		fprintf(fphisto[0],"%14.6Lf",(long double)psystem.histo_v[3][i]/psystem.vcount[3]*RESHISTO_v);
		fprintf(fphisto[0],"\n");
	}

	fprintf(fphisto[1],"# Vol \t P(Vol)\n");
	for(i=0;i<MBIN_Vol;i++){
		fprintf(fphisto[1],"%14.6f %14.6Lf\n",(double)(i-MBIN_Vol/2.0)/RESHISTO_Vol + psystem.histo_Vol0,
				(long double)psystem.histo_Vol[i]/psystem.Volcount*RESHISTO_Vol);
	}

	for(i=0;i<2;i++){
		fflush(fphisto[i]);
		fclose(fphisto[i]);
	}

}



//--------------------------------------init VACF--------------------------------------------------------------

void init_vacf(System &psystem, Particle *pparticle){

	int i,j,k;

	for(i=0;i<VACF_MAXT;i++){
		psystem.histo_vacf[i] = 0.0;
		psystem.vacfcount[i] = 0;
	}

	for(i=0;i<psystem.nparticles;i++){
		for(k=0;k<VACF_MAXT0;k++){
			for(j=0;j<3;j++){
				pparticle[i].vt0[j][k] = 0.0;
			}
		}
	}

	for(i=0;i<VACF_MAXT0;i++){
		psystem.vacf_t0time[i] = 0;
	}

	psystem.vacf_time = 0;
	psystem.vacf_countt0 = 0;

}

//-------------------------------------sample VACF-------------------------------------------------------------

void sample_vacf(System &psystem, Particle *pparticle){

	int i,j,k;
	int CorrelTime;
	int tmp;
	double vacf;

	//increase sampling time counte
	psystem.vacf_time++;

	//check if a new time origin should be set
	if((psystem.vacf_time % VACF_FREQT0) == 0){

		psystem.vacf_indext0 = psystem.vacf_countt0 % VACF_MAXT0;
		psystem.vacf_countt0++;
		psystem.vacf_t0time[psystem.vacf_indext0] = psystem.vacf_time;

		for(i=0;i<psystem.nparticles;i++){
			for(j=0;j<3;j++){
				pparticle[i].vt0[j][psystem.vacf_indext0] = pparticle[i].v[j];
			}
		}
	}

	//loop over all time origins and record vacf
	for(k=0;k<MIN(psystem.vacf_countt0,VACF_MAXT0,tmp);k++){
		CorrelTime = psystem.vacf_time - psystem.vacf_t0time[k];

		if(CorrelTime < VACF_MAXT){
			psystem.vacfcount[CorrelTime]++;
			vacf = 0.0;
			for (i=0;i<psystem.nparticles;i++){
				for(j=0;j<3;j++){
					vacf += pparticle[i].v[j]*pparticle[i].vt0[j][k];
				}
			}
			psystem.histo_vacf[CorrelTime] += vacf;
		}
	}
		
}

//---------------------------------------print vacf--------------------------------------------------------

void print_vacf(System &psystem){

	int i;
	FILE *fpvacf;
	double cumsum,deltat;
	fpvacf = fopen("vacf.dat","w");

	deltat = psystem.dt*VACF_SAMPLE;
	cumsum = 0.0;
	fprintf(fpvacf,"#time \t vacf \t sum_vacf\n");
	for(i=0;i<VACF_MAXT;i++){
		if(psystem.vacfcount[i]>0){
			psystem.histo_vacf[i] /= (double)(psystem.nparticles*psystem.vacfcount[i]);
		}
		else{
			psystem.histo_vacf[i] = 0.0;
		}
		cumsum += deltat*psystem.histo_vacf[i]/3.0;
		fprintf(fpvacf,"%14.6f %14.6f %14.6f\n",(i+1)*deltat,psystem.histo_vacf[i],cumsum);
	}

	fflush(fpvacf);
	fclose(fpvacf);

}

