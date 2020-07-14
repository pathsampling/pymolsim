//	Ssphere display color distribution:
//	0:	Green
//	1:	Cyan
//	2:	Red
//	3:	Yellow
//	4:	Blue
//	5:	Black


#ifndef DRAW3D_H
#define DRAW3D_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>
//linux
#include <GL/glut.h>
//#include "freeglut_ext.h"
//mac
//#include <GLUT/glut.h>
#include "header.h"

using namespace std;

//vector<Atom> pos;	/* Storage of sphere positions, etc */
int     Precision=8;        /* Number of sides per sphere, etc */
int     Quit = false;       /* Quit flag */
int     Run  = false;       /* Run flag */
int     Drawall =3;         /* Draw flag*/
//virtual particle

float	Far = 1000.0;		/* The far clipping plane */
float   Boxdistance =1.25;
float	g_Distance = -80;		/* Distance, Azimuth, Inclination and Twist */
//float	g_Distance = 0.0;		/* Distance, Azimuth, Inclination and Twist */
//GLfloat	xRotation = 90,		/* define the viewing position */
//        yRotation = 180,
//        zRotation = 0;
GLfloat	xRotation = 0.0,		/* define the viewing position */
        yRotation = 0.0,
        zRotation = 0;
//char string[] = "A";

//color & font set
static float MaterialRed[]     =   { 0.7, 0.0, 0.0, 1.0 };
static float MaterialGreen[]   =   { 0.1, 0.5, 0.2, 1.0 };
static float MaterialBlue[]    =   { 0.0, 0.0, 0.7, 1.0 };
static float MaterialYellow[]  =   { 0.7, 0.7, 0.0, 1.0 };
static float MaterialCyan[]    =   { 0.0, 0.7, 0.7, 1.0 };
static float MaterialBlack[]   =   { 0.0, 0.0, 0.0, 1.0 };

static float front_shininess[] =   {60.0};
static float front_specular[]  =   { 0.7, 0.7, 0.7, 1.0 };
static float ambient[]         =   { 0.0, 0.0, 0.0, 1.0 };
static float diffuse[]         =   { 1.0, 1.0, 1.0, 1.0 };
static float position0[]       =   { 1.0, 1.0, 1.0, 0.0 };
static float position1[]       =   {-1.0, -1.0, 1.0, 0.0 };
static float lmodel_ambient[]  =   { 0.5, 0.5, 0.5, 1.0 };
static float lmodel_twoside[]  =   {GL_TRUE};





/*
//drawing all spheres
void display_sph1(int ind) {
    glPushMatrix();
    glTranslatef(pos[ind].r.x, pos[ind].r.y, pos[ind].r.z);
    glScalef(pos[ind].radius, pos[ind].radius, pos[ind].radius);
    switch(pos[ind].no) {
        case 0:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialGreen)  ; break;
        case 1:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialCyan)  ; break;
        case 2:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialRed)   ; break;
        case 3:  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialYellow); break;
        case 4:	 glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlue); break;
        case 5:	 glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialBlack); break;
    }
    glutSolidSphere(1.0, Precision, Precision);
    glPopMatrix();
}

*/


/*
//for drawing a plane;
static float vert_memb[8][3] = {
    {-1.0, -1.0, 1.0},
    {-1.0, 1.0, 1.0},
    {1.0, 1.0, 1.0},
    {1.0, -1.0, 1.0},
    {-1.0, -1.0, -1.0},
    {-1.0, 1.0, -1.0},
    {1.0, 1.0, -1.0},
    {1.0, -1.0, -1.0}
};
*/

/*
//draw a plane
void display_plane(Vect &boxl){
    int i, j, N=int(ceil(boxl.x/lthMB));
    double lth=lthMB*N;    
    glPushMatrix();
    //glTranslatef(0, 0, 0);//center of ploting
    //glScalef(plane.x, plane.y, plane.z);
    //glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialCyan);
//    glColor3f(1.0, 1.0, 1.0);
//    glBegin(GL_QUADS);
//    for(i=0;i<8;i++) {
//        glVertex3fv(vert_memb[i]);
//    }
//    glEnd();
    glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialCyan);
    glTranslatef(0.0, 0.0, 0.0);
    //glTranslatef(-lth+lthMB,-lth+(N-1)*lthMB, 0.0);
    glTranslatef(-(lth-lthMB)/2.0-lthMB,-(lth-lthMB)/2.0, 0.0);
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            glTranslatef(lthMB, 0.0, 0.0);
            glutSolidCube(lthMB);
        }
        glTranslatef(-lth,lthMB, 0.0);
    }
    glPopMatrix();
}

*/

/*
//draw the membrane
void draw_memb(void){
    display_plane(Box);
}
*/

// reponse if the size of the window is changed.
static void Reshape( int width, int height ) {
	fpscreen<<"in Reshape"<<endl;
    glViewport( 0, 0, width, height );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glFrustum( -1.0, 1.0, -1.0, 1.0, 5.0, 1000.0 );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glTranslatef( 0.0, 0.0, -15.0 );
}




static void CursorKeys(int key, int x, int y) {
	fpscreen<<"in CursorKeys"<<endl;
    switch (key) {
        case GLUT_KEY_LEFT:
            zRotation += 5;
            break;
        case GLUT_KEY_RIGHT:
            zRotation -= 5;
            break;
        case GLUT_KEY_UP:
            xRotation += 5;
            break;
        case GLUT_KEY_DOWN:
            xRotation -= 5;
            break;
        default:
            return;
    }
    glutPostRedisplay();
}

static void Key(unsigned char key, int x, int y) {
	fpscreen<<"in Key"<<endl;
    switch (key) {
        case 27:
            exit(1);
        case ',':
            yRotation += 5;
            break;
        case '.':
            yRotation -= 5;
            break;
        case 'x':
            g_Distance += 2;
            break;
        case 'z':
            g_Distance -= 2;
            break;
        case '1':
            Precision =2;
            break;
        case '2':
            Precision =4;
            break;
        case '3':
            Precision =8;
            break;
	case '4':   
            Precision=16;
	    break;
	case '5':
	    Precision=32;
	    break;
        //case 't':
        //    beta+=0.1;
        //   fpscreen<<beta<<'\n';
        //    break;
        //case 'r':
        //    beta-=0.1;
        //    fpscreen<<beta<<'\n';
        //    break;
        default:
            return;
    }
    glutPostRedisplay();
}

static void Init( void ) {
    /* setup lighting, etc */
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position0);
    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, position1);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_specular);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glHint(GL_FOG_HINT, GL_FASTEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_FASTEST);
}

#define UNSCALED  1
#define NORMALIZE 2
#define RESCALE   3
#define QUIT      4

static void ModeMenu(int entry) {
    if (entry==UNSCALED) {
        //      glDisable(GL_RESCALE_NORMAL_EXT);
        glDisable(GL_NORMALIZE);
    }
    else if (entry==NORMALIZE) {
        glEnable(GL_NORMALIZE);
        //   glDisable(GL_RESCALE_NORMAL_EXT);
    }
    else if (entry==RESCALE) {
        glDisable(GL_NORMALIZE);
        //  glEnable(GL_RESCALE_NORMAL_EXT);
    }
    else if (entry==QUIT) {
        exit(0);
    }
    glutPostRedisplay();
}

void Init_Graphics(int argc, char *argv[]) {
    
    fpscreen<<"In graphics"<<endl;
    glutInit( &argc, argv);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(800, 800);
    
    glutInitDisplayMode( GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB );
    /*        glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );*/
    
    glutCreateWindow(argv[0]);
    glClearDepth(1.0);
    glClearColor( 1.0, 1.0, 1.0, 1.0 );
    glColor3f( 1.0, 1.0, 1.0 );
    
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glFlush();
    glutSwapBuffers();
    
    fpscreen<<"going in init"<<endl;
    Init();
    
    //fpscreen<<"allocating memory"<<endl;
    //essential part, number of particles
    //fpscreen<<"done allocating memory"<<endl;
    
    //convert_coordinates();
    
    glutIdleFunc( Idle );
    glutReshapeFunc( Reshape );
    glutSpecialFunc( CursorKeys );
    glutKeyboardFunc(Key);
    glutDisplayFunc( Display );
    
    glutCreateMenu(ModeMenu);
    glutAddMenuEntry("Unscaled", UNSCALED);
    glutAddMenuEntry("Normalize", NORMALIZE);
    glutAddMenuEntry("Rescale EXT", RESCALE);
    glutAddMenuEntry("Quit", QUIT);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    //fpscreen<<"before main loop"<<endl;
    glutMainLoop();
    //glutMainLoopEvent();
    //fpscreen<<"after main loop"<<endl;
    return;
}

#endif
