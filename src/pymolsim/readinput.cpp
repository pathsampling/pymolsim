#include "header.h"

//--------------------------------read input file-----------------------------------------------------------------


//void readinput(System &psystem, Potential &pljpot, Thermostat &pthermostat, Barostat &pbaro, int &pseed,int argc,char *argv[])
void readinput(System &psystem, Potential &pljpot, Thermostat &pthermostat, Barostat &pbaro, int &pseed)
{



	int i;
	//start reading the input file
    	//char flname[FILENAME_MAX+1]="input.dat";  //default file name
	string str;
	bool readstrucbool;
	readstrucbool = false;

	//parse the filename from first argument
    	//if (g_argc > 1) strcpy(flname,g_argv[1]); 
	
	//open input file
	ifstream fpinput (g_flname);
	if(!fpinput)
	{
		fpscreen<<"\t !ERROR! input file "<<g_flname<<" cannot be opened!\n";
		exit(1);
	}
	//Get the first 3 lines
	getline(fpinput,str);
	fpscreen<<"\n"<<str<<"\n";
	getline(fpinput,str);
	fpscreen<<str<<"\n";
	getline(fpinput,str);
	fpscreen<<str<<"\n";

	//Get rest of the input over a parser
	string string1;
	char instring[201],compstring[201];
	unsigned long linenumber = 0;
	getline(fpinput,string1);
	unsigned long run,runi,runstart;
	runstart=16;

	//do this for each line
	do{linenumber++;
		strcpy(instring,string1.c_str());
		char *stringpointer;
		unsigned short readcolumn=0;
		char *col[256];
		strcpy(compstring,instring);
		stringpointer = instring;
		//First 16 characters in col[0]
		//End the string at character no. 16
		instring[15]=0;
		//copy complete string to col[0]
		col[readcolumn]=stringpointer;
		//find the first not empty character to start parsing the rest
		for(runi=16;runi<201;runi++)
		{
			if(instring[runi]!=' '){
				stringpointer = &(instring[runi]);
				runstart=runi;
				break;
			}
			if(runi==200){
				fpscreen<<"\n\t!ERROR stringpointer could not be set\n\n";
				stringpointer=&(instring[16]);
				runstart=16;
			}
		}
		//stringpointer = &(instring[16]);  //this does not work if the 17th character is a empty!!
		readcolumn++;
		for(run=runstart;run<201;run++)
		{
			//if you are at the end of the line
			if(instring[run]==0)
			{
				col[readcolumn]=stringpointer;
				stringpointer = &(instring[run+1]);
				readcolumn++;
				break;
			}
			//get values separated by space ' '
			if(instring[run]==' ' && instring[run+1]!=' ')
			{
				instring[run]=0;
				col[readcolumn]=stringpointer;
				stringpointer=&(instring[run+1]);
				readcolumn++;
			}
		}

		// PROCESS THIS LINE (this needs to be done individually!!!)

  		// linenumber ... which line is processes
  		// readcolumn ... how many columns does this line have
  		// col[0] bis col[readcolumn-1] ... content of a respective column
  		// if column i contains a number: value is the given by atof(col[i-1]) ... as of type 'double'
		//                                                      atoi(col[i-1]) ... as of type 'int'

		if(0==strncmp(compstring,"ncycle1         ",16))
		{
			psystem.ncycle1=atoi(col[1]);
			fpscreen<<"ncycle1 \t\t = "<<psystem.ncycle1<<"\n";
			continue;
		}

		if(0==strncmp(compstring,"ncycle2         ",16))
		{
			psystem.ncycle2=atoi(col[1]);
			fpscreen<<"ncycle2 \t\t = "<<psystem.ncycle2<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"integrator      ",16))
		{
			psystem.integrator=atoi(col[1]);
			if(psystem.integrator == 0){
				fpscreen<<"integrator \t\t = velocity Verlet, NVE (option 0)\n";
			}
			else if(psystem.integrator == 1){
				fpscreen<<"integrator \t\t = Langevin thermostat, NVT (option 1)\n";
			}
			else if(psystem.integrator == 2){
				fpscreen<<"integrator \t\t = Andersen thermostat, NVT (option 2)\n";
			}
			else if(psystem.integrator == 3){
				fpscreen<<"integrator \t\t = Nose-Hoover thermostat, NVT (option 3)\n";
			}
			else if(psystem.integrator == 4){
				fpscreen<<"integrator \t\t = Nose-Hoover-Langevin thermostat, NVT (option 4)\n";
			}
			else if(psystem.integrator == 5){
				fpscreen<<"integrator \t\t = Andersen barostat, isotropic, NPH (option 5)\n";
				pbaro.isotropic = true;
			}
			else if(psystem.integrator == 6){
				fpscreen<<"integrator \t\t = Andersen barostat with Langevin friction piston, isotropic, NPT (option 6)\n";
				pbaro.isotropic = true;

			}
			else if(psystem.integrator == 7){
				fpscreen<<"integrator \t\t = Andersen barostat with Langevin friction piston + NHL thermostat, isotropic, NPT (option 7)\n";
				pbaro.isotropic = true;

			}
			else{
				fpscreen<<"\nNO PROPER INTEGRATOR CHOSEN\n";
				fpscreen<<"integrator variable = "<<psystem.integrator<<endl;
				fpscreen<<"\t\t...exiting program...\n\n";
				exit(0);
			}

			continue;
		}


		if(0==strncmp(compstring,"temperature     ",16))
		{
			psystem.beta=1.0/atof(col[1]);
			if(psystem.beta<=0.0){
				fpscreen<<"\n\n\t !ERROR! Reading input file: beta must be > 0.0!!\n";
				exit(1);
			}
			fpscreen<<"temperature \t\t = "<<1.0/psystem.beta<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"pressure        ",16))
		{
			psystem.pressure=atof(col[1]);
			fpscreen<<"pressure   \t\t = "<<psystem.pressure<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"energy          ",16))
		{
			psystem.tot_energy=atof(col[1]);
			//if(psystem.tot_energy<=0.0){
			//	fpscreen<<"\n\n\t !ERROR! Reading input file: total energy must be > 0.0!!\n";
			//	exit(1);
			//}
			fpscreen<<"total energy   \t\t = "<<psystem.tot_energy<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"switch plot     ",16))
		{
			psystem.swplot=atoi(col[1]);
			if(psystem.swplot<0 || psystem.swplot>1){
				fpscreen<<"\n\n\t !ERROR! Reading the input: plot switch must be 0 or 1!!\n";
				exit(1);
			}
			fpscreen<<"switch plot \t = "<<psystem.swplot<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"plot rate       ",16))
		{
			psystem.plotrate=atoi(col[1]);
			if(psystem.plotrate<0){
				fpscreen<<"\n\n\t !ERROR! Reading the input: plot rate must be larger than 0 !!\n";
				exit(1);
			}
			fpscreen<<"plot rate \t = "<<psystem.plotrate<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"outputrate      ",16))
		{
			psystem.outputrate=atoi(col[1]);
			if(psystem.outputrate<0){
				fpscreen<<"\n\n\t !ERROR! Reading the input: output rate must be larger than 0 !!\n";
				exit(1);
			}
			fpscreen<<"output rate \t\t = "<<psystem.outputrate<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"time_step       ",16))
		{
			psystem.dt=atof(col[1]);
			fpscreen<<"time step \t\t = "<<psystem.dt<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"particles       ",16))
		{
			if(readstrucbool){
				if(psystem.readstruc<1){
					fpscreen<<"no. of particles \t = structure read from input, see below\n";
				}
				else if(psystem.readstruc==1){
					psystem.nparticles=atoi(col[1]);
					fpscreen<<"no. of particles \t = "<<psystem.nparticles<<"\n";
				}
				else{
					fpscreen<<"\n\n\t!ERROR! Reading input file: something weired while reading nparticles\n";
					fpscreen<<"\tread_struc shoud be either 0 or 1, but it's = "<<psystem.readstruc<<endl;
					fpscreen<<"\tExiting program...\n";
					exit(1);
				}
			}
			else{
				fpscreen<<"\n\n\t!ERROR! Reading the input: must determine first, if structure is created or read in!\n";
				fpscreen<<"\tPut read_struc variable before particles, rho, and box!\n";
				fpscreen<<"\tExiting progam...\n";
				exit(1);
			}
			continue;
		}
		if(0==strncmp(compstring,"rho             ",16))
		{
			if(readstrucbool){
				if(psystem.readstruc<1){
					fpscreen<<"density rho \t\t = structure read from input, see below\n";
				}
				else if(psystem.readstruc==1){
					psystem.rho=atof(col[1]);
					fpscreen<<"density rho \t\t = "<<psystem.rho<<"\n";
				}
				else{
					fpscreen<<"\n\n\t!ERROR! Reading input file: something weired while reading rho\n";
					fpscreen<<"\tread_struc shoud be either 0 or 1, but it's = "<<psystem.readstruc<<endl;
					fpscreen<<"\tExiting program...\n";
					exit(1);
				}

			}
			else{
				fpscreen<<"\n\n\t!ERROR! Reading the input: must determine first, if structure is created or read in!\n";
				fpscreen<<"\tPut read_struc variable before particles, rho, and box!\n";
				fpscreen<<"\tExiting progam...\n";
				exit(1);
			}
			continue;

		}
		if(0==strncmp(compstring,"box             ",16))
		{
			if(readstrucbool){
				if(psystem.readstruc<1){
					fpscreen<<"relative box \t\t = structure read from input, see below\n";
				}
				else if(psystem.readstruc==1){
					for(i=0;i<3;i++){
						psystem.box_relative[i]=atof(col[i+1]);
					}
					fpscreen<<"relative box \t\t = "<<psystem.box_relative[0]<<"  "<<psystem.box_relative[1]<<"  "<<psystem.box_relative[2]<<"\n";
				}
				else{
					fpscreen<<"\n\n\t!ERROR! Reading input file: something weired while reading box\n";
					fpscreen<<"\tread_struc shoud be either 0 or 1, but it's = "<<psystem.readstruc<<endl;
					fpscreen<<"\tExiting program...\n";
					exit(1);
				}
			}
			else{
				fpscreen<<"\n\n\t!ERROR! Reading the input: must determine first, if structure is created or read in!\n";
				fpscreen<<"\tPut read_struc variable before particles, rho, and box!\n";
				fpscreen<<"\tExiting progam...\n";
				exit(1);
			}
			continue;
		}
		if(0==strncmp(compstring,"read_struc      ",16))
		{
			readstrucbool = true;
			psystem.readstruc=atoi(col[1]);
			fpscreen<<"read_struc \t\t = "<<psystem.readstruc;
			if(psystem.readstruc == 0){
				fpscreen<<"\tReading structure from lammps dump file format\n";
			}
			else if(psystem.readstruc == 1){
				fpscreen<<"\tCreating structure from nparticle, rho, and box\n";
			}			
			else{
				fpscreen<<"\n\n\t!ERROR! Reading the input: read_struc must be either 0, 1 or 2!\n";
				fpscreen<<"\tExiting program...\n";
				exit(1);
			}	
			continue;
		}
		if(0==strncmp(compstring,"epsilon         ",16))
		{
			pljpot.epsilon=atof(col[1]);
			fpscreen<<"LJ epsilon \t = "<<pljpot.epsilon<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"sigma           ",16))
		{
			pljpot.sigma=atof(col[1]);
			fpscreen<<"LJ sigma \t = "<<pljpot.sigma<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"Seed            ",16))
                {
                        pseed=atoi(col[1]);
                        fpscreen<<"Seed read from input \t = "<<pseed<<endl;
                        if(pseed>=0){
                                srand(time(NULL));
                                pseed=-(rand()%1000 +1);
                        }
                        fpscreen<<"Seed used \t\t = "<<pseed<<"\n";
                        continue;
                }
		if(0==strncmp(compstring,"langevin gamma  ",16))
		{
			pthermostat.lgamma=atof(col[1]);
			if(pthermostat.lgamma<=0.0){
				fpscreen<<"\n\n\t !ERROR! Reading input file: langevin gamma must be > 0.0!!\n";
				exit(1);
			}
			fpscreen<<"langevin gamma   \t\t = "<<pthermostat.lgamma<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"andersen nu     ",16))
		{
			pthermostat.anu=atof(col[1]);
			if(pthermostat.anu<=0.0){
				fpscreen<<"\n\n\t !ERROR! Reading input file: andersen nu must be > 0.0!!\n";
				exit(1);
			}
			fpscreen<<"andersen nu      \t\t = "<<pthermostat.anu<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"piston mass     ",16))
		{
			pbaro.pmass=atof(col[1]);
			fpscreen<<"piston mass      \t\t = "<<pbaro.pmass<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"piston gamma    ",16))
		{
			pbaro.lgamma=atof(col[1]);
			fpscreen<<"piston gamma      \t\t = "<<pbaro.lgamma<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"nhl gamma       ",16))
		{
			pthermostat.nhlgamma=atof(col[1]);
			if(pthermostat.nhlgamma<0.0){
				fpscreen<<"\n\n\t !ERROR! Reading input file: nose-hoover-langevin gamma must be > 0.0!!\n";
				exit(1);
			}
			fpscreen<<"NHL gamma   \t\t = "<<pthermostat.nhlgamma<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"nhl mu          ",16))
		{
			pthermostat.nhlmu=atof(col[1]);
			if(pthermostat.nhlmu<=0.0){
				fpscreen<<"\n\n\t !ERROR! Reading input file: nose-hoover-langevin mu must be > 0.0!!\n";
				exit(1);
			}
			fpscreen<<"NHL mu     \t\t = "<<pthermostat.nhlmu<<"\n";
			continue;
		}
		if(0==strncmp(compstring,"rdf resolution  ",16))
		{
			psystem.rdfRES=(int)rint(1.0/atof(col[1]));
			fpscreen<<"resolution of rdf      \t\t = "<<(double)(1.0/psystem.rdfRES)<<"\n";
			continue;
		}








		//print out text string
		if(0==strncmp(compstring,"#",1))
		{
			fpscreen<<"\n"<<compstring<<"\n";
		}




	}while (NULL!=getline(fpinput,string1));  //do until end of the file

	fpscreen.setf(ios::fmtflags(0),ios::floatfield);
	fpscreen<<"\nEND READING THE INPUT FILE '"<<g_flname<<"'\n";
	fpscreen<<"\n------------------------------------------------------------------------------------\n";



	return;
}

//-----------------------------------------read in initial structure-----------------------------------------------------------------------

void read_initstructure(System &psystem, Particle *pparticle){
	
	int i,j;
	double dummy;
	char dummy_char[256];
	ifstream fpinput("structure.in");
	bool directcoord;
	
	if(!fpinput)
	{
		fpscreen<<"\n\t !ERROR! input file structure.in cannot be opened!\n";
		exit(1);
	}


	fpinput.getline(dummy_char,256);
	fpinput.getline(dummy_char,256);
	//only for orthogonal boxes right now!
	fpinput >> psystem.box[0];
	fpinput >> dummy;
	fpinput >> dummy;
	
	fpinput >> dummy;
	fpinput >> psystem.box[1];
	fpinput >> dummy;

	fpinput >> dummy;
	fpinput >> dummy;
	fpinput >> psystem.box[2];

	fpinput >> psystem.nparticles;
	
	for(i=0;i<3;i++){
		psystem.box_2[i] = psystem.box[i]/2.0;
	}
	psystem.box_relative[0] = 1.0;
	psystem.box_relative[1] = psystem.box[1]/psystem.box[0];
	psystem.box_relative[2] = psystem.box[2]/psystem.box[0];
	//calculate density
	psystem.volume = psystem.box[0]*psystem.box[1]*psystem.box[2];
	psystem.rho = psystem.nparticles/psystem.volume;
	
	fpinput.getline(dummy_char,256);
	fpinput.getline(dummy_char,256);
	directcoord = false;
	if((0==strncmp(dummy_char,"C",1)) || (0==strncmp(dummy_char,"c",1))){
		fpscreen<<"\nC(c)artesien coordinates\n";
		//exit(1);
	}
	else if((0==strncmp(dummy_char,"D",1)) || (0==strncmp(dummy_char,"d",1))){
		fpscreen<<"\nD(d)irect coordinates\n";
		directcoord = true;
		//exit(1);
	}
	else{
		fpscreen<<"\n\n\t!ERROR! Reading input structure\n";
		fpscreen<<"\tInput format not recognized;  must be C(c)artesian or D(d)irect\n";
		fpscreen<<"\tExiting program...\n";
		exit(1);
	}
	//fpscreen<<"dummy char = "<<dummy_char<<endl;
	//fpscreen<<endl;

	//allocate array of particles
	//pparticle = new Particle [psystem.nparticles];

	if(psystem.nparticles > MAX_PART){
		fpscreen<<"\n\n\t !ERROR! Reading input structure: no of particles larger than MAX_PART\n";
		fpscreen<<"\t MAX_PART = "<<MAX_PART<<"\t nparticles = "<<psystem.nparticles<<endl;
		fpscreen<<"\tExiting program...\n";
		exit(1);
	}

	for(i=0;i<psystem.nparticles;i++){
		for(j=0;j<3;j++){
			fpinput >> pparticle[i].r[j];
		}
	}

	if(directcoord){
		fpscreen<<"converting postions\n";
		for(i=0;i<psystem.nparticles;i++){
			for(j=0;j<3;j++){
				pparticle[i].r[j] = pparticle[i].r[j]*psystem.box[j];
			}
		}
		//exit(1);
	}
	else{
		fpscreen<<"not converting positions\n";
		//exit(1);
	}
}


//-----------------------------------------read in initial structure in lammps dump format-------------------------------------------------

void read_initlammpsdump(System &psystem, Particle *pparticle){
	
	int i,j;
	double dummy;
	char dummy_char[256];
	ifstream fpinput("conf.dump");
	
	if(!fpinput)
	{
		fpscreen<<"\n\t !ERROR! input file conf.dump cannot be opened!\n";
		fpscreen<<"\t Exiting program...\n";
		exit(1);
	}


	fpinput.getline(dummy_char,256);
	fpinput.getline(dummy_char,256);
	fpinput.getline(dummy_char,256);
	
	//number of particles
	fpinput >> psystem.nparticles;
	fpinput.getline(dummy_char,256);
	fpinput.getline(dummy_char,256);


	//BOX: only for orthogonal boxes right now!
	fpinput >> dummy;
	fpinput >> psystem.box[0];
	
	fpinput >> dummy;
	fpinput >> psystem.box[1];

	fpinput >> dummy;
	fpinput >> psystem.box[2];
	fpinput.getline(dummy_char,256);

	
	for(i=0;i<3;i++){
		psystem.box_2[i] = psystem.box[i]/2.0;
	}
	psystem.box_relative[0] = 1.0;
	psystem.box_relative[1] = psystem.box[1]/psystem.box[0];
	psystem.box_relative[2] = psystem.box[2]/psystem.box[0];
	//calculate density
	psystem.volume = psystem.box[0]*psystem.box[1]*psystem.box[2];
	psystem.rho = psystem.nparticles/psystem.volume;
	

	if(psystem.nparticles > MAX_PART){
		fpscreen<<"\n\n\t !ERROR! Reading input structure: no of particles larger than MAX_PART\n";
		fpscreen<<"\t MAX_PART = "<<MAX_PART<<"\t nparticles = "<<psystem.nparticles<<endl;
		fpscreen<<"\tExiting program...\n";
		exit(1);
	}

	fpinput.getline(dummy_char,256);
	
	for(i=0;i<psystem.nparticles;i++){
		fpinput >> dummy;
		fpinput >> dummy;
		fpinput >> pparticle[i].mass;
		for(j=0;j<3;j++){
			fpinput >> pparticle[i].r[j];
		}
		for(j=0;j<3;j++){
			fpinput >> pparticle[i].v[j];
		}
	}
}


