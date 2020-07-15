#include "header.h"

//--------------------------------------------calculate forces----------------------------------------------------
void LJ::forces(vector<double> dr){
	//assign the 
	int i,j,k,m,n;
	double r2,r6i;
	double f,g;
	double sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;

	virial = 0.0;
	hypervirial = 0.0;
	for(m=0;m<3;m++){
		for(n=0;n<3;n++){
			stress[3*m+n] = 0.0;
		}
	}
	for(j=0;j<3;j++){
		fi[j] = 0.0;
		fj[j] = 0.0;
	}

	r2 = SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]);	
	if(r2 <= rmin2){
		r6i = 1.0/(r2*r2*r2);
		f = 24.0*epsilon*sigma6*r6i/r2* (2.0*sigma6*r6i - 1.0);
		g = 48.0*epsilon*sigma6*r6i * (13.0*sigma6*r6i - 3.5);

		virial += f*r2;
		hypervirial += (g-f*r2);
		for(m=0;m<3;m++){
			for(n=0;n<3;n++){
				stress[3*m+n] += f*dr[m]*dr[n];
			}
		}
		// force: \sum -dV/dr * 1/r * (xi - xj) = f * xij = -f * xji  
		for(k=0;k<3;k++){
			fi[k] -= f*dr[k];	// -f * xji
			fj[k] += f*dr[k];	/// f * xji
		}
	}
	else if(r2 < rmax2){
		r6i = 1.0/(r2*r2*r2);
		f = 2.0/r2 * (6.0*c2 * SQR(sigma6*r6i) + 3.0*c3*sigma6*r6i) - 2.0*c4/(sigma*sigma);
		// calculate g = d^2V/dr^2 * r^2
		g = 6.0*sigma6*r6i*(c2*26.0*sigma6*r6i + 7.0*c3) + 2.0*c4*r2/(sigma*sigma);
		//virial part of the pressure
		virial += f*r2;
		//hypervirial
		hypervirial += (g-f*r2);

		//stress tensor \sigma_mn
		for(m=0;m<3;m++){
			for(n=0;n<3;n++){
				stress[3*m+n] += f*dr[m]*dr[n];
			}
		}	

		// force: \sum -dV/dr * 1/r * (xi - xj) = f * xij = -f * xji  
		for(k=0;k<3;k++){
			fi[k] -= f*dr[k];
			fj[k] += f*dr[k];
		}

	}
	else{
		f = 0.0;
		g = 0.0;
		//virial part of the pressure
		virial += f*r2;
		//hypervirial
		hypervirial += (g-f*r2);

		//stress tensor \sigma_mn
		for(m=0;m<3;m++){
			for(n=0;n<3;n++){
				stress[3*m+n] += f*dr[m]*dr[n];
			}
		}	
		// force: \sum -dV/dr * 1/r * (xi - xj) = f * xij = -f * xji  
		for(k=0;k<3;k++){
			fi[k] -= 0.0;
			fj[k] += 0.0;
		}		
	}
	//TODO : scaling for virial is not yet done
}

void LJ::potential_energy(vector<double> dr)
{
	int i,j;
	double r2,r6i;
	double sigma6;

	energy = 0.0;
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;

	//loop over all particle pairs
	r2 = SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]);
			
	if(r2 <= rmin2){		// r <= 2.3 \sigma
		r6i = 1.0/(r2*r2*r2);
		energy += (4.0*epsilon*sigma6*r6i)*(sigma6*r6i - 1.0) + c1;
	}			
	else if(r2 < rmax2){	// 2.3 \sigma < r < 2.5 \sigma
		r6i = 1.0/(r2*r2*r2);
		energy += c2*SQR(sigma6*r6i) + c3*sigma6*r6i + c4*r2/(sigma*sigma) + c5;
	}
	else{				// r >= 2.5 \sigma
		energy += 0.0;
	}
}


//Now we need a equivalent wrapper implementation in system
void Sim::forces()
{
	//cout<<"I ran0";
	//initialize forces
	//TODO : Set nparticles
	for(int i=0;i<nparticles;i++){
		for(int j=0;j<3;j++){
			particles[i].f[j] = 0.0;
		}
	}
	
	//loop over all particle pairs
	for(int i=0;i<nparticles-1;i++){
		for(int j=i+1;j<nparticles;j++){
			//get distance between particles considering PBC
			// xji = xj - xi in this case
			// dr = rj - ri
			//TODO : Modify 
			auto dr = image_distance(i, j);
			//call forces here
			potential.forces(dr);
			//now we need to update the forces
			//cout<<potential.fi[0]<<endl;
			for(int c=0; c<3; c++){
				particles[i].f[c] += potential.fi[c];
				particles[j].f[c] += potential.fj[c];
			}
		}
	}
	//cout<<"I ran";
	for(int m=0;m<3;m++){
		for(int n=0;n<3;n++){
			stress[3*m+n] = potential.stress[3*m+n]/volume;
		}
	}

	//scale pressure
	virial = potential.virial/(3.0*volume);
	//scale hypervirial
	hypervirial = potential.hypervirial/9.0;

}

//---------------------------------------------calculate potential energy-----------------------------------------

void Sim::potential_energy()
{
	pe = 0.0;
	//loop over all particle pairs
	for(int i=0;i<nparticles-1;i++){
		for(int j=i+1;j<nparticles;j++){
			//get distance between particles considering PBC
			auto dr = image_distance(i, j);
			potential.potential_energy(dr);
			pe += potential.energy;
		}
	}
}

//-------------------------------------------------calculate kinetic energy---------------------------------------

void Sim::kinetic_energy()
{
	
	ke = 0.0;
	
	for(int i=0;i<nparticles;i++){
		ke += particles[i].mass*(SQR(particles[i].v[0]) + SQR(particles[i].v[1]) + SQR(particles[i].v[2]));
	}

	ke /= 2.0;

}



//-----------------------------------------------calculate total energy-------------------------------------------

void Sim::total_energy()
{
	double K,V;
	kinetic_energy();
	potential_energy();
	tot_energy = pe + ke;
}

//-----------------------------------------------calculate total momentum-----------------------------------------

void Sim::total_momentum()
{

	for(int i=0;i<3;i++){
		p[i] = 0.0;
	}

	for(int i=0;i<nparticles;i++){
		for(int j=0;j<3;j++){
			p[j] += particles[i].v[j]*particles[i].mass;
		}
	}
}

//-------------------------------------------------rescale velocities to desired temperature--------------------------
void Sim::rescale_velocities()
{
	int i,j;
	double T,T_aim,c;

	kinetic_energy();
	T = 2.0*ke/((3.0*(double)nparticles)-3.0);
	T_aim = 1.0/beta;
	c = sqrt(T_aim/T);
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] *= c;
		}
	}
}

