#include "header.h"


//--------------------------------------do an MD step with the chosen integrator----------------------------------

void md_step(System &psystem, Particle *pparticle, Potential &pljpot, Thermostat &pthermostat, Barostat &pbaro)
{
	if(psystem.integrator == 0){
		//velocity verlet, NVE
		md_verlet(psystem,pparticle,pljpot);
	}
	else if (psystem.integrator == 1){
		//langevin thermostat, NVT
		md_langevin(psystem,pparticle,pljpot,pthermostat);
	}
	else if (psystem.integrator == 2){
		//andersen thermostat, NVT
		md_andersen(psystem,pparticle,pljpot,pthermostat);
	}
	else if (psystem.integrator == 3){
		//nose-hoover thermostat, NVT (with c1 = 1.0 and c2 = 0.0)
		md_nosehooverlangevin_NVT(psystem,pparticle,pljpot,pthermostat);
	}
	else if (psystem.integrator == 4){
		//nose-hoover-langevin thermostat, NVT
		md_nosehooverlangevin_NVT(psystem,pparticle,pljpot,pthermostat);
	}
	else if (psystem.integrator == 5){
		//andersen barostat, NPH
		md_andersen_NPH(psystem,pparticle,pljpot,pbaro);
	}
	else if (psystem.integrator == 6){
		//andersen barostat with langevin piston, NPT
		md_andersen_stochastic_NPT(psystem,pparticle,pljpot,pbaro);
	}
	else if (psystem.integrator == 7){
		//andersen barostat with langevin piston + nose-hoover langevin thermostat, NPT
		md_andersen_stochastic_nhlthermo_NPT(psystem,pparticle,pljpot,pbaro,pthermostat);
	}
	else{
		fpscreen<<"\n\t\tNO PROPER INTEGRATOR DEFINED\n";
		fpscreen<<"\t\t...exiting program...\n\n";
		exit(0);
	}
}

//------------------------------------------do a velocity verlet move---------------------------------------------

void md_verlet(System &psystem, Particle *pparticle,Potential &pljpot)
{
	//do half step algorithm (a bit redundant, but just for proof that it works

	propagate_momenta_half(psystem,pparticle);
	propagate_position_half(psystem,pparticle);
	propagate_position_half(psystem,pparticle);
	
	forces(pparticle,psystem,pljpot);
	
	propagate_momenta_half(psystem,pparticle);

}

//---------------------------------------------do langevin integration step--------------------------------------------

void md_langevin(System &psystem, Particle *pparticle,Potential &pljpot, Thermostat &pthermo)
{

	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//first half of stochastic advance of momenta
	langevin_thermo(psystem,pparticle,pthermo);

	//update positions and velocities according to 'normal' velocity verlet
	md_verlet(psystem,pparticle,pljpot);

	//second half of stochastic propagation of momenta
	langevin_thermo(psystem,pparticle,pthermo);


}


//---------------------------------------------do langevin integration step--------------------------------------------

void md_andersen(System &psystem, Particle *pparticle,Potential &pljpot, Thermostat &pthermo)
{

	int i,j;
	//update positions and velocities according to 'normal' velocity verlet
	md_verlet(psystem,pparticle,pljpot);

	//couple to andersen thermostat
	for(i=0;i<psystem.nparticles;i++){
		if(ran3() < (pthermo.anu*psystem.dt)){
			for(j=0;j<3;j++){
				pparticle[i].v[j] = sqrt(1.0/(pparticle[i].mass*psystem.beta))*gauss();
			}
		}
	}



}

//------------------------------------------do a Nose-Hoover-Langevin move----------------------------------------------------------------

void md_nosehooverlangevin_NVT(System &psystem, Particle *pparticle,Potential &pljpot, Thermostat &pthermo)
{
	double K,Nf;
	//do half step algorithm
	Nf = (3.0*(double)psystem.nparticles)-3.0;

	//propagate xi first half time step
	pthermo.nhlxi = pthermo.nhlc1*pthermo.nhlxi + pthermo.nhlc2*gauss();
	//propagate momenta with xi
	propagate_momenta_xi(psystem,pparticle,pthermo);
	//1/2 Verlet step
	propagate_momenta_half(psystem,pparticle);
	propagate_position_half(psystem,pparticle);

	//propagate xi full time step with kinetic energy (middle step)
	K = kinetic_energy(pparticle,psystem);
	pthermo.nhlxi = pthermo.nhlxi + psystem.dt*(2.0*K - (Nf/psystem.beta))/pthermo.nhlmu;
	
	//second half step
	propagate_position_half(psystem,pparticle);
	//get forces for new positions to propagate momenta	
	forces(pparticle,psystem,pljpot);
	
	propagate_momenta_half(psystem,pparticle);
	propagate_momenta_xi(psystem,pparticle,pthermo);
	pthermo.nhlxi = pthermo.nhlc1*pthermo.nhlxi + pthermo.nhlc2*gauss();

}





//----------------------------------------do Andersen NPH integration step----------------------------------------------

void md_andersen_NPH(System &psystem, Particle *pparticle, Potential &pljpot, Barostat &pbaro)
{
	//constant pressure extended system approach by Andersen
	//Ref:  Andersen, JCP 72, 2384 (1980)
	//      Kolb, Duenweg, JCP 111, 4453 (1999) -> implementation, symplectic integrator


	int i,j;
	double dt;
	double K;
	double Vold,Vnew;
	double box[3];

	Vold = psystem.volume;
	dt = psystem.dt;
	//update momenta with f(t)
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] += 0.5*pparticle[i].f[j]/pparticle[i].mass*dt;
		}
	}

	//isotropic barostat
	if(pbaro.isotropic == true){
		//fpscreen<<"Isotropic barostat\n";
		
		//evaluate pressure with r(t), L(t) and the new velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		K = kinetic_energy(pparticle,psystem);
		psystem.Pinst = psystem.virial + 2.0*(double)psystem.nparticles*K/(((3.0*(double)psystem.nparticles)-3.0)*psystem.volume);

		//update half time step momentum of box variables/volume
		//pbaro.pp = pbaro.pp + (psystem.Pinst - psystem.pressure)*0.5*dt;
		//pbaro.pv = pbaro.pp/pbaro.pmass;
		pbaro.pv = pbaro.pv + (psystem.Pinst - psystem.pressure)*0.5*dt/pbaro.pmass;

		//update volume, half time step and get box side length
		Vnew = Vold + 0.5*dt*pbaro.pv;
		for(i=0;i<3;i++){
			box[i] = psystem.box_relative[i] * pow((Vnew/(psystem.box_relative[0]*psystem.box_relative[1]*psystem.box_relative[2])),(1.0/3.0));
		}
		//update positions
		for(i=0;i<psystem.nparticles;i++){
			for(j=0;j<3;j++){
				pparticle[i].r[j] += SQR(psystem.box[j]/box[j])*pparticle[i].v[j]*dt;
			}
		}

		//update second half step for volume and update box side length
		Vnew = Vnew + 0.5*dt*pbaro.pv;
		for(i=0;i<3;i++){
			box[i] = psystem.box_relative[i] * pow((Vnew/(psystem.box_relative[0]*psystem.box_relative[1]*psystem.box_relative[2])),(1.0/3.0));
		}

		//rescale positions and velocities and update system box size and volume
		for(i=0;i<psystem.nparticles;i++){
			for(j=0;j<3;j++){
				pparticle[i].r[j] = box[j]/psystem.box[j] * pparticle[i].r[j];  
				pparticle[i].v[j] = psystem.box[j]/box[j] * pparticle[i].v[j];  
			}
		}
		for(i=0;i<3;i++){
			psystem.box[i] = box[i];
			psystem.box_2[i] = psystem.box[i]/2.0;
		}
		psystem.volume = psystem.box[0]*psystem.box[1]*psystem.box[2];

		//evaluate forces and virial contribution to the pressure with new positions
		forces(pparticle,psystem,pljpot);

		//evaluate pressure with r(t+dt), L(t+dt) and the half updated/rescaled velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		K = kinetic_energy(pparticle,psystem);
		psystem.Pinst = psystem.virial + 2.0*(double)psystem.nparticles*K/(((3.0*(double)psystem.nparticles)-3.0)*psystem.volume);

		//update second half time step momentum of box variables/volume
		//pbaro.pp = pbaro.pp + (psystem.Pinst - psystem.pressure)*0.5*dt;
		//pbaro.pv = pbaro.pp/pbaro.pmass;
		pbaro.pv = pbaro.pv + (psystem.Pinst - psystem.pressure)*0.5*dt/pbaro.pmass;



	}
	else{
		// !! only update in z-direction !!

		//fpscreen<<"Anisotropic barostat\n";
		//evaluate pressure with r(t), L(t) and the new velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		K = kinetic_energy(pparticle,psystem);
		psystem.Pinst = psystem.virial + 2.0*(double)psystem.nparticles*K/(((3.0*(double)psystem.nparticles)-3.0)*psystem.volume);

		//update half time step momentum of box variables/volume
		pbaro.pv = pbaro.pv + (psystem.Pinst - psystem.pressure)*0.5*dt/pbaro.pmass;

		//update volume, half time step and get box side length
		Vnew = Vold + 0.5*dt*pbaro.pv;
		//update only box in z-direction
		box[2] = Vnew/(psystem.box[0]*psystem.box[1]);
		//update positions
		for(i=0;i<psystem.nparticles;i++){
			for(j=0;j<2;j++){
				pparticle[i].r[j] += pparticle[i].v[j]*dt;
			}
			pparticle[i].r[2] += SQR(psystem.box[2]/box[2])*pparticle[i].v[2]*dt;
		}

		//update second half step for volume and update box side length
		Vnew = Vnew + 0.5*dt*pbaro.pv;
		box[2] = Vnew/(psystem.box[0]*psystem.box[1]);

		//rescale positions and velocities and update system box size and volume
		for(i=0;i<psystem.nparticles;i++){
			pparticle[i].r[2] = box[2]/psystem.box[2] * pparticle[i].r[2];  
			pparticle[i].v[2] = psystem.box[2]/box[2] * pparticle[i].v[2];  
		}
		psystem.box[2] = box[2];
		psystem.box_2[2] = psystem.box[2]/2.0;
		psystem.volume = psystem.box[0]*psystem.box[1]*psystem.box[2];

		//evaluate forces and virial contribution to the pressure with new positions
		forces(pparticle,psystem,pljpot);

		//evaluate pressure with r(t+dt), L(t+dt) and the half updated/rescaled velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		K = kinetic_energy(pparticle,psystem);
		psystem.Pinst = psystem.virial + 2.0*(double)psystem.nparticles*K/(((3.0*(double)psystem.nparticles)-3.0)*psystem.volume);

		//update second half time step momentum of box variables/volume
		pbaro.pv = pbaro.pv + (psystem.Pinst - psystem.pressure)*0.5*dt/pbaro.pmass;


	}

	//final update of momenta with new forces
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] += 0.5*pparticle[i].f[j]/pparticle[i].mass*dt;
		}
	}

}


//----------------------------------------do Andersen-stochastic NPT integration step----------------------------------------------

void md_andersen_stochastic_NPT(System &psystem, Particle *pparticle, Potential &pljpot, Barostat &pbaro)
{

	// due to coupling of the piston to a heat bath, this is formally an NPT ensemble
	// since coupling only through one variable (the piston), coupling might be slow
	// Langevin piston method:  Feller, Zhang, Pastor, Brooks, JCP 103, 4613 (1995)


	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//first half of stochastic advance of piston momentum
	langevin_baro(pbaro);

	//update positions and velocities according to Andersen NPH
	md_andersen_NPH(psystem,pparticle,pljpot,pbaro);

	//second half of stochastic propagation of piston momentum
	langevin_baro(pbaro);
}

//----------------------------------------stochastic anderson barostat + NHL thermostat---------------------------------------------------
void md_andersen_stochastic_nhlthermo_NPT(System &psystem,Particle *pparticle,Potential &pljpot,Barostat &pbaro,Thermostat &pthermo)
{
	int i;
	double K,Nf;
	double dt;
	double Vnew;
	double box[3],scale[3];

	Vnew = psystem.volume;
	dt = psystem.dt;

	
	Nf = (3.0*(double)psystem.nparticles)-3.0;
	
	//propagte barostat velocity 1/2 time step
	langevin_baro(pbaro);
	//propagate xi 1/2 time step
	langevin_xi(pthermo);
	//propagate momenta with xi
	propagate_momenta_xi(psystem,pparticle,pthermo);
	//1/2 Verlet step momenta
	propagate_momenta_half(psystem,pparticle);

	//update kinetic energy and evaluate Pinst
	//evaluate pressure with r(t), L(t) and the new velocities
	//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
	K = kinetic_energy(pparticle,psystem);
	psystem.Pinst = psystem.virial + 2.0*(double)psystem.nparticles*K/(Nf*psystem.volume);

	//update half time step momentum of box variables/volume
	pbaro.pv = pbaro.pv + (psystem.Pinst - psystem.pressure)*0.5*dt/pbaro.pmass;
	//update volume, half time step and get box side length
	Vnew = Vnew + 0.5*dt*pbaro.pv;
	if(pbaro.isotropic == true){
		for(i=0;i<3;i++){
			box[i] = psystem.box_relative[i] * pow((Vnew/(psystem.box_relative[0]*psystem.box_relative[1]*psystem.box_relative[2])),(1.0/3.0));
		}
	}
	else{
		box[0] = psystem.box[0];
		box[1] = psystem.box[1];
		box[2] = Vnew/(psystem.box[0]*psystem.box[1]);
	}
	for(i=0;i<3;i++){
		scale[i] = SQR(psystem.box[i]/box[i]);
	}
	//propagate positions 1/2 time step
	propagate_position_half_scale(psystem,pparticle,scale);
	//propagate xi full time step with kinetic energy (middle step)
	pthermo.nhlxi = pthermo.nhlxi + psystem.dt*(2.0*K - (Nf/psystem.beta))/pthermo.nhlmu;

	//second half of symplectic integrator
	//propagate positions 1/2 time step
	propagate_position_half_scale(psystem,pparticle,scale);
	//update volume, half time step and get box side length
	Vnew = Vnew + 0.5*dt*pbaro.pv;
	if(pbaro.isotropic == true){
		for(i=0;i<3;i++){
			box[i] = psystem.box_relative[i] * pow((Vnew/(psystem.box_relative[0]*psystem.box_relative[1]*psystem.box_relative[2])),(1.0/3.0));
		}
	}
	else{
		box[0] = psystem.box[0];
		box[1] = psystem.box[1];
		box[2] = Vnew/(psystem.box[0]*psystem.box[1]);
	}
	for(i=0;i<3;i++){
		scale[i] = psystem.box[i]/box[i];
	}
	//rescale positions and velocities
	rescale_position_momenta(psystem,pparticle,scale);
	//update system box and volume
	for(i=0;i<3;i++){
		psystem.box[i] = box[i];
		psystem.box_2[i] = psystem.box[i]/2.0;
	}
	psystem.volume = psystem.box[0]*psystem.box[1]*psystem.box[2];

	//evaluate forces and virial contribution to the pressure with new positions
	forces(pparticle,psystem,pljpot);
	//evaluate pressure with r(t+dt), L(t+dt) and the half updated/rescaled velocities
	//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
	K = kinetic_energy(pparticle,psystem);
	psystem.Pinst = psystem.virial + 2.0*(double)psystem.nparticles*K/(Nf*psystem.volume);
	//update second half time step momentum of box variables/volume
	pbaro.pv = pbaro.pv + (psystem.Pinst - psystem.pressure)*0.5*dt/pbaro.pmass;
	//1/2 Verlet step momenta
	propagate_momenta_half(psystem,pparticle);
	//propagate momenta with xi
	propagate_momenta_xi(psystem,pparticle,pthermo);
	//propagate xi 1/2 time step
	langevin_xi(pthermo);
	//propagte barostat velocity 1/2 time step
	langevin_baro(pbaro);

}




//------------------------------------------langevin thermostat of velocities------------------------------------------

void langevin_thermo(System &psystem, Particle *pparticle, Thermostat &pthermo)
{
	int i,j;

	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//update velocities and account for langevin contribution to the kinetic energy (to be substracted later)
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			psystem.dElangevin += 0.5*pparticle[i].mass*SQR(pparticle[i].v[j]);
			pparticle[i].v[j] = pthermo.lc1*pparticle[i].v[j] + pthermo.lc2/sqrt(pparticle[i].mass)*gauss();
			psystem.dElangevin -= 0.5*pparticle[i].mass*SQR(pparticle[i].v[j]);
		}
	}
}


//------------------------------------------langevin thermostat of piston velocity-------------------------------------

void langevin_baro(Barostat &pbaro)
{

	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//update piston velocity 
	//is there a hidden contribution to the velocities somewhere??

	//only isotropic case for now
	pbaro.pv = pbaro.lc1*pbaro.pv + pbaro.lc2*gauss();
}



//---------------------------------------propagate momenta 1/2 time step------------------------------------------------------------------

void propagate_momenta_half(System &psystem, Particle *pparticle)
{
	double dt;
	int i,j;

	dt = psystem.dt;

	//loop over all particles
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] += 0.5*pparticle[i].f[j]/pparticle[i].mass*dt;
		}
	}

}

//----------------------------------------------propagate positions 1/2 time step---------------------------------------------------------

void propagate_position_half(System &psystem, Particle *pparticle)
{
	double dt;
	int i,j;

	dt = psystem.dt;

	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].r[j] += 0.5*dt*pparticle[i].v[j];

		}
	}

}

//---------------------------------------propagate momenta 1/2 time step with xi----------------------------------------------------------

void propagate_momenta_xi(System &psystem, Particle *pparticle, Thermostat &pthermo)
{
	double c;
	int i,j;

	c = exp(-pthermo.nhlxi*psystem.dt/2.0);

	//loop over all particles
	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].v[j] = pparticle[i].v[j]*c ;
		}
	}

}

//------------------------------------------langevin thermostat of nose-hover extended variable------------------------

void langevin_xi(Thermostat &pthermo)
{
	pthermo.nhlxi = pthermo.nhlc1*pthermo.nhlxi + pthermo.nhlc2*gauss();
}


//----------------------------------------------propagate positions 1/2 time step and scale-----------------------------------------------

void propagate_position_half_scale(System &psystem, Particle *pparticle, double *pscale)
{
	double dt;
	int i,j;

	dt = psystem.dt;

	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].r[j] += 0.5*dt*pparticle[i].v[j]*pscale[j];

		}
	}

}

//----------------------------------------------rescale positions and momenta for barostat------------------------------------------------

void rescale_position_momenta(System &psystem, Particle *pparticle, double *pscale)
{
	int i,j;

	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			pparticle[i].r[j] /= pscale[j];
			pparticle[i].v[j] *= pscale[j];
		}
	}

}




