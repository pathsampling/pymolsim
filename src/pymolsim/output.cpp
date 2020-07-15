#include "header.h"


//--------------------------------initialize energies and velocities for read in structure---------------------------


void Sim::dump(int step){

	fstream fout;
	if (step == 0){
		fout.open(trajfile, fstream::out);
	}
	else{
		fout.open(trajfile, fstream::app);	
	}

	int i,j;
	fout<<"ITEM: TIMESTEP"<<endl;
	fout<<step<<endl;
	fout<<"ITEM: NUMBER OF ATOMS"<<endl;
	fout<<nparticles<<endl;
	fout<<"ITEM: BOX BOUNDS pp pp pp"<<endl;
	fout<<"0 "<<box[0]<<endl;
	fout<<"0 "<<box[1]<<endl;
	fout<<"0 "<<box[2]<<endl;
	fout<<"ITEM: ATOMS id type mass x y z vx vy vz"<<endl;
//	pfplammps.setf(ios::scientific);
//	pfplammps.precision(6);

	for(i=0;i<nparticles;i++){
		fout<<i+1<<" "<<particles[i].type<<" "<<particles[i].mass<<" ";
		for(j=0;j<3;j++){
			fout<<particles[i].r[j]<<" ";
		}
		for(j=0;j<3;j++){
			fout<<particles[i].v[j]<<" ";
		}
		fout<<endl;
	}

	fout.flush();
	fout.close();	

}




