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

