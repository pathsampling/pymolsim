/*
 * draw3D.c - OpenGL interface 
 */
#include "draw3D.h"
//#include "header.h"


//exampler for drawing various things
void draw_stuff(){
    /* make a delay for a short period. */
    //sleep(4);
    
    return;
}

//JR: stuff that was in Lizhe's main programm
//static void Idle() {  //originally a static function, not sure how this works now
void Idle() {
	double box[3];
	double alpha;
//    if (iBlock % N_mem == 0) {
//        fpscreen << "Block cycle:" << iBlock << "\n";
//    }
//    cycle_iter();
//    convert_coordinates();
    	//glutPostRedisplay();
//    iBlock++;
    //ending Monte Carlo simulation
//    if (iBlock > Nblock) {
//        fpscreen << "\n\nMonte Carlo sampling finished!";
        //ending sampling scheme
        //store information
//        infoStore();
//        exit(EXIT_SUCCESS);
//    }
	fpscreen<<"in Idle\n";
	fpscreen<<"calling main progam...\n\n";
	main_program();
	exit(0);
	
	alpha=-20.0;
	int i;
	for(i=0;i<10000;i++){
		fpscreen<<i<<" - Idle"<<endl;
		box[0] = 20.0;
		box[1] = 10.0;
		box[2] = 80.0;
		//box[2] = 40.0 + (double) i*0.20;
		alpha -= 0.20;
	       	glLoadIdentity();
    		glTranslatef( 0.0, 0.0, g_Distance );
    		//glRotatef(xRotation, 1, 0, 0);
    		//glRotatef(yRotation, 0, 1, 0);
    		//glRotatef(zRotation, 0, 0, 1);
    		//  glPolarView(g_Distance,Azimuth,Inclination,Twist);    
    		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    		//draw box
    		glMatrixMode( GL_MODELVIEW );
		//glRotated(90.0,0.0,1.0,0.0);
		glRotated(alpha,1.0,0.0,0.0);
		glRotated(90.0,0.0,1.0,0.0);
  	 	 draw_box(box,0);
    		//draw membrane
//  	  	draw_memb();
    		//Drawall, controls whether to draw all the particles
//    		if (Drawall>0) {
//    	    		for(i=0;i<Nsphere;i++)   display_sph1(i);
//    		}
		//glRotated(90.0,0.0,1.0,0.0);
    		glFlush();
    		glutSwapBuffers();
	
    	//	glutPostRedisplay();
//		Display();
		//sleep(1);
	}
	exit(0);

}

//-----------------------------------------------------update display with current system configuration-------------------------------

void update_display(System &psystem,Particle *pparticle){

	//static double alpha =  10.0;
	static double alpha =  0.0;
	int i;
	double pos[3];
	char buffer[50];
	
	glLoadIdentity();
    	glMatrixMode( GL_MODELVIEW );
    	glTranslatef( 0.0, 0.0, g_Distance );
    	//glPolarView(g_Distance,Azimuth,Inclination,Twist);    
    	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//rotate view
	glRotated(alpha,1.0,0.0,0.0);
	glRotated(90.0,0.0,1.0,0.0);
    	//draw box
  	draw_box(psystem.box,0);
    		//draw membrane
//  	  	draw_memb();
    		//Drawall, controls whether to draw all the particles
	for(i=0;i<psystem.nparticles;i++)
	{
		convert_coordinates(psystem,pparticle[i],pos);
		display_sph1(pos);
	}
	glLoadIdentity();
    	glTranslatef( 0.0, 0.0, g_Distance );
	sprintf(buffer,"time = %14.6f",psystem.time);
	print_test(buffer);
    	glFlush();
    	glutSwapBuffers();
	//alpha -= 0.20;

	return;
}


//-----------------------------------------------convert coordinates suitable for display---------------------------------------------

void convert_coordinates(System &psystem,Particle &pparticle_i,double *pos){
	int j;
	for(j=0;j<3;j++){
		pos[j] = pparticle_i.r[j] - psystem.box_2[j];
	}

}

//--------------------------------------------display a sphere------------------------------------------------------------------------
void display_sph1(double *ppos){
	double radius;
//	radius = 0.5;
	radius = 0.3;
    	glPushMatrix();
    	glTranslated(ppos[0], ppos[1], ppos[2]);
    	glScaled(radius, radius, radius);
    	//switch(pos[ind].no) {
        //	case 0:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialGreen)  ; break;
        //	case 1:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialCyan)  ; break;
        //	case 2:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialRed)   ; break;
        //	case 3:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialYellow); break;
        //	case 4:	 glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlue); break;
        //	case 5:	 glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlack); break;
    	//}
	glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialRed);
    	glutSolidSphere(1.0, Precision, Precision);
    	glPopMatrix();
}

//------------------------------------------------display function of OpenGL----------------------------------------------------------
void Display() {
//    int i;
    fpscreen<<"in Display"<<endl;
    //for(i=0;i<10;i++){
//	    fpscreen<<i<<" - Display"<<endl;
    	//glutLeaveMainLoop();
    	glLoadIdentity();
    	glTranslatef( 0.0, 0.0, g_Distance );
    	glRotatef(xRotation, 1, 0, 0);
    	glRotatef(yRotation, 0, 1, 0);
    	glRotatef(zRotation, 0, 0, 1);
    	//  glPolarView(g_Distance,Azimuth,Inclination,Twist);    
    	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    	//draw box
//  	  draw_box( Box, 0 );
    	//draw membrane
//  	  draw_memb();
    	//Drawall, controls whether to draw all the particles
//    	if (Drawall>0) {
//    	    for(i=0;i<Nsphere;i++)   display_sph1(i);
//    	}
    	glFlush();
    	glutSwapBuffers();
	//for(i=0;i<5;i++){
	//    fpscreen<<i<<" - Display"<<endl;
	//sleep(2);
	//}
    //}
}


//for drawing box
static float vert_box[17][3] = {
    {-1.0, -1.0, -1.0},
    {-1.0, -1.0, 1.0},
    { 1.0, -1.0, 1.0},
    { 1.0, -1.0, -1.0},
    {-1.0, -1.0, -1.0},
    {-1.0, 1.0, -1.0},
    { 1.0, 1.0, -1.0},
    { 1.0, -1.0, -1.0},
    { 1.0, -1.0, 1.0},
    { 1.0, 1.0, 1.0},
    { 1.0, 1.0, 1.0},
    { 1.0, 1.0, -1.0},
    { 1.0, 1.0, 1.0},
    {-1.0, 1.0, 1.0},
    {-1.0, 1.0, -1.0},
    {-1.0, 1.0, 1.0},
    {-1.0, -1.0, 1.0},
};


void draw_box(double *pbox,double boxdis) {
    int i;
    glPushMatrix();
    glScaled(0.5*pbox[0], 0.5*pbox[1], 0.5*pbox[2]);
    glTranslated(0, boxdis, 0);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlack);
    glColor3f(1.0, 0, 0);
    glBegin(GL_LINE_STRIP);
   // glRotated(90.0,0.0,1.0,0.0);
    for(i=0;i<17;i++) {
        glVertex3fv(vert_box[i]);
    }
    glEnd();
    glPopMatrix();
}


void print_test(char* text){
	glPushMatrix();
	glTranslated(-10,7,0);
	glScaled(-g_Distance/10000,-g_Distance/10000,-g_Distance/10000);
	for(char* p = text; *p; p++){
		glutStrokeCharacter(GLUT_STROKE_ROMAN,*p);
	}
	glPopMatrix();
}



