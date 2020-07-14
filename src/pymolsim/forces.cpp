#include "header.h"

//--------------------------------------------calculate forces----------------------------------------------------

void forces(Particle *pparticle, System &psystem, Potential &pljpot)
{
	int i,j,k,m,n;
	double dr[3];
	double r2,r6i;
	double f,g;
	double epsilon,sigma,sigma6;

	epsilon = pljpot.epsilon;
	sigma = pljpot.sigma;
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;

	//initialize forces
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].f[j] = 0.0;
		}
	}
	psystem.virial = 0.0;
	psystem.hypervirial = 0.0;
	for(m=0;m<3;m++){
		for(n=0;n<3;n++){
			psystem.stress[3*m+n] = 0.0;
		}
	}

	//loop over all particle pairs
	for(i=0;i<psystem.nparticles-1;i++){
		for(j=i+1;j<psystem.nparticles;j++){
			//get distance between particles considering PBC
			// xji = xj - xi in this case
			// dr = rj - ri
			image_distance(pparticle[i],pparticle[j],dr,psystem);
			r2 = SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]);
			if(r2 <= pljpot.rmin2){
				//calculate f = -dV/dr * 1/r
				r6i = 1.0/(r2*r2*r2);
				f = 24.0*epsilon*sigma6*r6i/r2* (2.0*sigma6*r6i - 1.0);
				// calculate g = d^2V/dr^2 * r^2
				g = 48.0*epsilon*sigma6*r6i * (13.0*sigma6*r6i - 3.5);
				//virial part of the pressure
				psystem.virial += f*r2;
				//hypervirial
				psystem.hypervirial += (g-f*r2);
				//stress tensor \sigma_mn
				for(m=0;m<3;m++){
					for(n=0;n<3;n++){
						psystem.stress[3*m+n] += f*dr[m]*dr[n];
					}
				}	

				// force: \sum -dV/dr * 1/r * (xi - xj) = f * xij = -f * xji  
				for(k=0;k<3;k++){
					pparticle[i].f[k] -= f*dr[k];	// -f * xji

					pparticle[j].f[k] += f*dr[k];	/// f * xji
				}
			}
			else if(r2 < pljpot.rmax2){
				r6i = 1.0/(r2*r2*r2);
				f = 2.0/r2 * (6.0*pljpot.c2 * SQR(sigma6*r6i) + 3.0*pljpot.c3*sigma6*r6i) - 2.0*pljpot.c4/(sigma*sigma);
				// calculate g = d^2V/dr^2 * r^2
				g = 6.0*sigma6*r6i* (pljpot.c2*26.0*sigma6*r6i + 7.0*pljpot.c3) + 2.0*pljpot.c4*r2/(sigma*sigma);
				//virial part of the pressure
				psystem.virial += f*r2;
				//hypervirial
				psystem.hypervirial += (g-f*r2);

				//stress tensor \sigma_mn
				for(m=0;m<3;m++){
					for(n=0;n<3;n++){
						psystem.stress[3*m+n] += f*dr[m]*dr[n];
					}
				}	

				// force: \sum -dV/dr * 1/r * (xi - xj) = f * xij = -f * xji  
				for(k=0;k<3;k++){
					pparticle[i].f[k] -= f*dr[k];
				
					pparticle[j].f[k] += f*dr[k];
				}
			}
			else{
				f = 0.0;
				g = 0.0;
				//virial part of the pressure
				psystem.virial += f*r2;
				//hypervirial
				psystem.hypervirial += (g-f*r2);

				//stress tensor \sigma_mn
				for(m=0;m<3;m++){
					for(n=0;n<3;n++){
						psystem.stress[3*m+n] += f*dr[m]*dr[n];
					}
				}	
				// force: \sum -dV/dr * 1/r * (xi - xj) = f * xij = -f * xji  
				for(k=0;k<3;k++){
					pparticle[i].f[k] -= 0.0;

					pparticle[j].f[k] += 0.0;
				}
			}
		}
	}
	//scale pressure
	psystem.virial /= 3.0*psystem.volume;
	//scale hypervirial
	psystem.hypervirial /= 9.0;
	//scale stress tensor
	for(m=0;m<3;m++){
		for(n=0;n<3;n++){
			psystem.stress[3*m+n] /= psystem.volume;
		}
	}

	

}

//---------------------------------------------calculate potential energy-----------------------------------------

double potential_energy(Particle *pparticle,System &psystem,Potential &pljpot)
{
        double energy;
	int i,j;
	double dr[3];
	double r2,r6i;
	double epsilon,sigma,sigma6;

	energy = 0.0;
	epsilon = pljpot.epsilon;
	sigma = pljpot.sigma;
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;

	//loop over all particle pairs
	for(i=0;i<psystem.nparticles-1;i++){
		for(j=i+1;j<psystem.nparticles;j++){
			//get distance between particles considering PBC
			image_distance(pparticle[i],pparticle[j],dr,psystem);
			r2 = SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]);
			if(r2 <= pljpot.rmin2){		// r <= 2.3 \sigma
				r6i = 1.0/(r2*r2*r2);
				energy += (4.0*epsilon*sigma6*r6i)*(sigma6*r6i - 1.0) + pljpot.c1;
			}
			else if(r2 < pljpot.rmax2){	// 2.3 \sigma < r < 2.5 \sigma
				r6i = 1.0/(r2*r2*r2);
				energy += pljpot.c2*SQR(sigma6*r6i) + pljpot.c3*sigma6*r6i + pljpot.c4*r2/(sigma*sigma) + pljpot.c5;
			}
			else{				// r >= 2.5 \sigma
				energy += 0.0;
			}
		}
	}

	return (energy);
}

//-------------------------------------------------calculate kinetic energy---------------------------------------

double kinetic_energy(Particle *pparticle, System &psystem)
{
	double energy;
	int i;

	energy = 0.0;
	
	for(i=0;i<psystem.nparticles;i++){
		energy += pparticle[i].mass*(SQR(pparticle[i].v[0]) + SQR(pparticle[i].v[1]) + SQR(pparticle[i].v[2]));
	}

	energy /= 2.0;

	return (energy);
}



//-----------------------------------------------calculate total energy-------------------------------------------

double total_energy(Particle *pparticle, System &psystem, Potential &pljpot)
{
	double K,V;
	K = kinetic_energy(pparticle,psystem);
	V = potential_energy(pparticle,psystem,pljpot);

	return (K+V);
}

//-----------------------------------------------calculate total momentum-----------------------------------------

void total_momentum(Particle *pparticle,double *p,System &psystem)
{

	int i,j;

	for(i=0;i<3;i++){
		p[i] = 0.0;
	}

	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			p[j] += pparticle[i].v[j]*pparticle[i].mass;
		}
	}
}

//-------------------------------------------------rescale velocities to desired temperature--------------------------
void rescale_velocities(Particle *pparticle,System &psystem)
{
	double K;
	int i,j;
	double T,T_aim,c;
	K = kinetic_energy(pparticle,psystem);
	T = 2.0*K/((3.0*(double)psystem.nparticles)-3.0);
	T_aim = 1.0/psystem.beta;
	c = sqrt(T_aim/T);
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] *= c;
		}
	}
}

