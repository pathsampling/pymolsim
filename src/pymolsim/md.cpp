#include "header.h"


//--------------------------------------do an MD step with the chosen integrator----------------------------------

void Sim::md_step()
{
	if(integrator == 0){
		//velocity verlet, NVE
		//cout<<"----------------"<<endl;
		//cout<<particles[0].r[0]<<endl;
		md_verlet();
		//cout<<particles[0].r[0]<<endl;
	}
	else if (integrator == 1){
		//langevin thermostat, NVT
		md_langevin();
	}
	else if (integrator == 2){
		//andersen thermostat, NVT
		md_andersen();
	}
	else if (integrator == 3){
		//nose-hoover thermostat, NVT (with c1 = 1.0 and c2 = 0.0)
		md_nosehooverlangevin_NVT();
	}
	else if (integrator == 4){
		//nose-hoover-langevin thermostat, NVT
		md_nosehooverlangevin_NVT();
	}
	else if (integrator == 5){
		//andersen barostat, NPH
		md_andersen_NPH();
	}
	else if (integrator == 6){
		//andersen barostat with langevin piston, NPT
		md_andersen_stochastic_NPT();
	}
	else if (integrator == 7){
		//andersen barostat with langevin piston + nose-hoover langevin thermostat, NPT
		md_andersen_stochastic_nhlthermo_NPT();
	}
	else{
		exit(0);
	}
}

//------------------------------------------do a velocity verlet move---------------------------------------------

void Sim::md_verlet()
{
	//do half step algorithm (a bit redundant, but just for proof that it works

	propagate_momenta_half();
	propagate_position_half();
	propagate_position_half();
	
	forces();
	
	propagate_momenta_half();

}

//---------------------------------------------do langevin integration step--------------------------------------------

void Sim::md_langevin()
{

	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//first half of stochastic advance of momenta
	langevin_thermo();

	//update positions and velocities according to 'normal' velocity verlet
	md_verlet();

	//second half of stochastic propagation of momenta
	langevin_thermo();


}


//---------------------------------------------do langevin integration step--------------------------------------------

void Sim::md_andersen()
{

	int i,j;
	//update positions and velocities according to 'normal' velocity verlet
	md_verlet();

	//couple to andersen thermostat
	for(i=0;i<nparticles;i++){
		if(ran3() < (thermostat.anu*dt)){
			for(j=0;j<3;j++){
				particles[i].v[j] = sqrt(1.0/(particles[i].mass*beta))*gauss();
			}
		}
	}



}

//------------------------------------------do a Nose-Hoover-Langevin move----------------------------------------------------------------

void Sim::md_nosehooverlangevin_NVT()
{
	double K,Nf;
	//do half step algorithm
	Nf = (3.0*(double)nparticles)-3.0;

	//propagate xi first half time step
	thermostat.nhlxi = thermostat.nhlc1*thermostat.nhlxi + thermostat.nhlc2*gauss();
	//propagate momenta with xi
	propagate_momenta_xi();
	//1/2 Verlet step
	propagate_momenta_half();
	propagate_position_half();

	//propagate xi full time step with kinetic energy (middle step)
	kinetic_energy();
	thermostat.nhlxi = thermostat.nhlxi + dt*(2.0*ke - (Nf/beta))/thermostat.nhlmu;
	
	//second half step
	propagate_position_half();
	//get forces for new positions to propagate momenta	
	forces();
	
	propagate_momenta_half();
	propagate_momenta_xi();
	thermostat.nhlxi = thermostat.nhlc1*thermostat.nhlxi + thermostat.nhlc2*gauss();

}





//----------------------------------------do Andersen NPH integration step----------------------------------------------

void Sim::md_andersen_NPH()
{
	//constant pressure extended system approach by Andersen
	//Ref:  Andersen, JCP 72, 2384 (1980)
	//      Kolb, Duenweg, JCP 111, 4453 (1999) -> implementation, symplectic integrator


	int i,j;

	double K;
	double Vold,Vnew;
	double lbox[3];

	Vold = volume;

	//update momenta with f(t)
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] += 0.5*particles[i].f[j]/particles[i].mass*dt;
		}
	}

	//isotropic barostat
	if(barostat.isotropic == true){
		//fpscreen<<"Isotropic barostat\n";
		
		//evaluate pressure with r(t), L(t) and the new velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		kinetic_energy();
		Pinst = virial + 2.0*(double)nparticles*ke/(((3.0*(double)nparticles)-3.0)*volume);

		//update half time step momentum of box variables/volume
		//pbaro.pp = pbaro.pp + (psystem.Pinst - psystem.pressure)*0.5*dt;
		//pbaro.pv = pbaro.pp/pbaro.pmass;
		barostat.pv = barostat.pv + (Pinst - pressure)*0.5*dt/barostat.pmass;

		//update volume, half time step and get box side length
		Vnew = Vold + 0.5*dt*barostat.pv;
		for(i=0;i<3;i++){
			lbox[i] = box_relative[i] * pow((Vnew/(box_relative[0]*box_relative[1]*box_relative[2])),(1.0/3.0));
		}
		//update positions
		for(i=0;i<nparticles;i++){
			for(j=0;j<3;j++){
				particles[i].r[j] += SQR(box[j]/lbox[j])*particles[i].v[j]*dt;
			}
		}

		//update second half step for volume and update box side length
		Vnew = Vnew + 0.5*dt*barostat.pv;
		for(i=0;i<3;i++){
			lbox[i] = box_relative[i] * pow((Vnew/(box_relative[0]*box_relative[1]*box_relative[2])),(1.0/3.0));
		}

		//rescale positions and velocities and update system box size and volume
		for(i=0;i<nparticles;i++){
			for(j=0;j<3;j++){
				particles[i].r[j] = lbox[j]/box[j] * particles[i].r[j];  
				particles[i].v[j] = box[j]/lbox[j] * particles[i].v[j];  
			}
		}
		for(i=0;i<3;i++){
			box[i] = lbox[i];
			box_2[i] = box[i]/2.0;
		}
		volume = box[0]*box[1]*box[2];

		//evaluate forces and virial contribution to the pressure with new positions
		forces();

		//evaluate pressure with r(t+dt), L(t+dt) and the half updated/rescaled velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		kinetic_energy();
		Pinst = virial + 2.0*(double)nparticles*ke/(((3.0*(double)nparticles)-3.0)*volume);

		//update second half time step momentum of box variables/volume
		//pbaro.pp = pbaro.pp + (psystem.Pinst - psystem.pressure)*0.5*dt;
		//pbaro.pv = pbaro.pp/pbaro.pmass;
		barostat.pv = barostat.pv + (Pinst - pressure)*0.5*dt/barostat.pmass;



	}
	else{
		// !! only update in z-direction !!

		//fpscreen<<"Anisotropic barostat\n";
		//evaluate pressure with r(t), L(t) and the new velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		kinetic_energy();
		Pinst = virial + 2.0*(double)nparticles*ke/(((3.0*(double)nparticles)-3.0)*volume);

		//update half time step momentum of box variables/volume
		barostat.pv = barostat.pv + (Pinst - pressure)*0.5*dt/barostat.pmass;

		//update volume, half time step and get box side length
		Vnew = Vold + 0.5*dt*barostat.pv;
		//update only box in z-direction
		lbox[2] = Vnew/(box[0]*box[1]);
		//update positions
		for(i=0;i<nparticles;i++){
			for(j=0;j<2;j++){
				particles[i].r[j] += particles[i].v[j]*dt;
			}
			particles[i].r[2] += SQR(box[2]/lbox[2])*particles[i].v[2]*dt;
		}

		//update second half step for volume and update box side length
		Vnew = Vnew + 0.5*dt*barostat.pv;
		lbox[2] = Vnew/(box[0]*box[1]);

		//rescale positions and velocities and update system box size and volume
		for(i=0;i<nparticles;i++){
			particles[i].r[2] = lbox[2]/box[2] * particles[i].r[2];  
			particles[i].v[2] = box[2]/lbox[2] * particles[i].v[2];  
		}
		box[2] = lbox[2];
		box_2[2] = box[2]/2.0;
		volume = box[0]*box[1]*box[2];

		//evaluate forces and virial contribution to the pressure with new positions
		forces();

		//evaluate pressure with r(t+dt), L(t+dt) and the half updated/rescaled velocities
		//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
		kinetic_energy();
		Pinst = virial + 2.0*(double)nparticles*ke/(((3.0*(double)nparticles)-3.0)*volume);

		//update second half time step momentum of box variables/volume
		barostat.pv = barostat.pv + (Pinst - pressure)*0.5*dt/barostat.pmass;


	}

	//final update of momenta with new forces
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] += 0.5*particles[i].f[j]/particles[i].mass*dt;
		}
	}

}


//----------------------------------------do Andersen-stochastic NPT integration step----------------------------------------------

void Sim::md_andersen_stochastic_NPT()
{

	// due to coupling of the piston to a heat bath, this is formally an NPT ensemble
	// since coupling only through one variable (the piston), coupling might be slow
	// Langevin piston method:  Feller, Zhang, Pastor, Brooks, JCP 103, 4613 (1995)


	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//first half of stochastic advance of piston momentum
	langevin_baro();

	//update positions and velocities according to Andersen NPH
	md_andersen_NPH();

	//second half of stochastic propagation of piston momentum
	langevin_baro();
}

//----------------------------------------stochastic anderson barostat + NHL thermostat---------------------------------------------------
void Sim::md_andersen_stochastic_nhlthermo_NPT()
{
	int i;
	double K,Nf;

	double Vnew;
	double lbox[3],scale[3];

	Vnew = volume;

	
	Nf = (3.0*(double)nparticles)-3.0;
	
	//propagte barostat velocity 1/2 time step
	langevin_baro();
	//propagate xi 1/2 time step
	langevin_xi();
	//propagate momenta with xi
	propagate_momenta_xi();
	//1/2 Verlet step momenta
	propagate_momenta_half();

	//update kinetic energy and evaluate Pinst
	//evaluate pressure with r(t), L(t) and the new velocities
	//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
	kinetic_energy();
	Pinst = virial + 2.0*(double)nparticles*ke/(Nf*volume);

	//update half time step momentum of box variables/volume
	barostat.pv = barostat.pv + (Pinst - pressure)*0.5*dt/barostat.pmass;
	//update volume, half time step and get box side length
	Vnew = Vnew + 0.5*dt*barostat.pv;
	if(barostat.isotropic == true){
		for(i=0;i<3;i++){
			lbox[i] = box_relative[i] * pow((Vnew/(box_relative[0]*box_relative[1]*box_relative[2])),(1.0/3.0));
		}
	}
	else{
		lbox[0] = box[0];
		lbox[1] = box[1];
		lbox[2] = Vnew/(box[0]*box[1]);
	}
	for(i=0;i<3;i++){
		scale[i] = SQR(box[i]/lbox[i]);
	}
	//propagate positions 1/2 time step
	propagate_position_half_scale(scale);
	
	//propagate xi full time step with kinetic energy (middle step)
	thermostat.nhlxi = thermostat.nhlxi + dt*(2.0*ke - (Nf/beta))/thermostat.nhlmu;

	//second half of symplectic integrator
	//propagate positions 1/2 time step
	propagate_position_half_scale(scale);
	//update volume, half time step and get box side length
	Vnew = Vnew + 0.5*dt*barostat.pv;
	if(barostat.isotropic == true){
		for(i=0;i<3;i++){
			lbox[i] = box_relative[i] * pow((Vnew/(box_relative[0]*box_relative[1]*box_relative[2])),(1.0/3.0));
		}
	}
	else{
		lbox[0] = box[0];
		lbox[1] = box[1];
		lbox[2] = Vnew/(box[0]*box[1]);
	}
	for(i=0;i<3;i++){
		scale[i] = box[i]/lbox[i];
	}
	//rescale positions and velocities
	rescale_position_momenta(scale);
	//update system box and volume
	for(i=0;i<3;i++){
		box[i] = lbox[i];
		box_2[i] = box[i]/2.0;
	}
	volume = box[0]*box[1]*box[2];

	//evaluate forces and virial contribution to the pressure with new positions
	forces();
	//evaluate pressure with r(t+dt), L(t+dt) and the half updated/rescaled velocities
	//virial part of the pressure is done in the force routine and should be okay here, only new kinetic part
	kinetic_energy();
	Pinst = virial + 2.0*(double)nparticles*ke/(Nf*volume);
	//update second half time step momentum of box variables/volume
	barostat.pv = barostat.pv + (Pinst - pressure)*0.5*dt/barostat.pmass;
	//1/2 Verlet step momenta
	propagate_momenta_half();
	//propagate momenta with xi
	propagate_momenta_xi();
	//propagate xi 1/2 time step
	langevin_xi();
	//propagte barostat velocity 1/2 time step
	langevin_baro();

}




//------------------------------------------langevin thermostat of velocities------------------------------------------

void Sim::langevin_thermo()
{
	int i,j;

	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//update velocities and account for langevin contribution to the kinetic energy (to be substracted later)
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			thermostat.dElangevin += 0.5*particles[i].mass*SQR(particles[i].v[j]);
			particles[i].v[j] = thermostat.lc1*particles[i].v[j] + thermostat.lc2/sqrt(particles[i].mass)*gauss();
			thermostat.dElangevin -= 0.5*particles[i].mass*SQR(particles[i].v[j]);
		}
	}
}


//------------------------------------------langevin thermostat of piston velocity-------------------------------------

void Sim::langevin_baro()
{

	//Ref:  Bussi, Parrinello, PRE 75, 056707 (2007)	
	//update piston velocity 
	//is there a hidden contribution to the velocities somewhere??

	//only isotropic case for now
	barostat.pv = barostat.lc1*barostat.pv + barostat.lc2*gauss();
}



//---------------------------------------propagate momenta 1/2 time step------------------------------------------------------------------

void Sim::propagate_momenta_half()
{

	int i,j;

	//loop over all particles
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] += 0.5*particles[i].f[j]/particles[i].mass*dt;
		}
	}

}

//----------------------------------------------propagate positions 1/2 time step---------------------------------------------------------

void Sim::propagate_position_half()
{

	int i,j;

	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			//cout<<particles[i].r[j]<<" "<<particles[i].v[j]<<" "<<dt<<endl;
			particles[i].r[j] += 0.5*dt*particles[i].v[j];
			//cout<<particles[i].r[j]<<" "<<particles[i].v[j]<<endl;

		}
	}

}

//---------------------------------------propagate momenta 1/2 time step with xi----------------------------------------------------------

void Sim::propagate_momenta_xi()
{
	double c;
	int i,j;

	c = exp(-thermostat.nhlxi*dt/2.0);

	//loop over all particles
	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].v[j] = particles[i].v[j]*c ;
		}
	}

}

//------------------------------------------langevin thermostat of nose-hover extended variable------------------------

void Sim::langevin_xi()
{
	thermostat.nhlxi = thermostat.nhlc1*thermostat.nhlxi + thermostat.nhlc2*gauss();
}


//----------------------------------------------propagate positions 1/2 time step and scale-----------------------------------------------

void Sim::propagate_position_half_scale(double *pscale)
{

	int i,j;

	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].r[j] += 0.5*dt*particles[i].v[j]*pscale[j];

		}
	}

}

//----------------------------------------------rescale positions and momenta for barostat------------------------------------------------

void Sim::rescale_position_momenta(double *pscale)
{
	int i,j;

	for(i=0;i<nparticles;i++){
		for(j=0;j<3;j++){
			particles[i].r[j] /= pscale[j];
			particles[i].v[j] *= pscale[j];
		}
	}

}




