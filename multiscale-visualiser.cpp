//mostly drawn from Computer Graphics by Hill and Kelley, Pearson International, 2007

#include <cstdlib>                        // standard definitions
#include <iostream>                        // C++ I/O
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>                        // C I/O (for sprintf) 
#include <cmath>                        // standard definitions
#include <vector>
#include <map>
#include <GL/glut.h>                        // GLUT
#include <GL/glu.h>                        // GLU
#include <GL/gl.h>                        // OpenGL

struct Point3D
{
    GLfloat x,y,z;
    Point3D(GLfloat _x, GLfloat _y, GLfloat _z) : x(_x),y(_y),z(_z) {};
};
struct ColorTriple
{
    float color[3];
};

//-----------------------------------------------------------------------
// Global variables
//-----------------------------------------------------------------------
double rotAngle = 0;                    // rotation angle for camera
double xpos = 6.0, ypos = 6.4, zpos = 4.2;
double zdistance = 2.0;
double cell[3]; //PBC of simulation box
int nparticles;     //number of particles in inputCoordStream
int nparticles_vis; //number of particles to be visualised, as listed in inputInfoStream
int nbonds;     //number of backbone harmonic springs, as listed in inputBondStream
int nnbs;       //number of non-bonded harmonic springs, as listed in inputBondStream
std::vector<Point3D> coordinates;
std::vector<std::vector<int> > bondsList, nonBondsList;
std::vector<int> indices_vis; //indices of particles to be included in visualisation
std::vector<std::string> parttypes; //types of particles in indices_vis
std::vector<float> sphereSizes; //sizes of particles in indices_vis
std::vector<ColorTriple> sphereColors; //colors of particles in indices_vis
//predifined maps for plotting
std::map<std::string, float> sizeMap = { {"A", 0.1}, {"C", 0.05}, {"H", 0.05}, {"O", 0.05}, {"S", 0.05}, {"N", 0.05} };
ColorTriple red = {1.0,0.0,0.0}, blue = {0.0, 0.0, 1.0}, white = {1.0, 1.0, 1.0}, cyan = {0.0, 1.0, 1.0}, grey = {0.5, 0.5, 0.5}, yellow = {1.0, 1.0, 0.0};
std::map<std::string, ColorTriple> colorMap = { {"A", grey}, {"C", cyan}, {"H", white}, {"O", red}, {"S", yellow}, {"N", blue} };
//O is red, N is blue, C is cyan and S is yellow
std::ifstream inputCoordStream, inputBondStream, inputInfoStream;

//-----------------------------------------------------------------------
// init (sets up some default OpenGL values)
//-----------------------------------------------------------------------

void init()
{
    glClearColor(1.0, 1.0, 1.0, 0);            // background color
    glClearDepth(1.0);                   // background depth value

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //gluPerspective(60, 1, 1, 1000);        // setup a perspective projection, alternatively glOrtho()
    glOrtho(-2.4,2.4,-2.4,2.4,0,50);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(                            // set up the camera
         //       0.0, 0.0, zdistance,      // eye position
         //       5.0, 3.0, 1.0,            // lookat position
         //       0.0, 1.0, 0.0);           // up direction
    xpos,ypos,zpos,3,3,3,0.5,1,0);

    glEnable(GL_DEPTH_TEST);            // enable hidden surface removal
    glEnable(GL_LIGHTING);              // enable lighting
    glEnable(GL_LIGHT0);                // enable lightsource 0
    float lpos[] = { 5, 5, 5, 0 };
    glLightfv(GL_LIGHT0, GL_POSITION, lpos);
    glShadeModel(GL_SMOOTH);            // smooth shading
}

void drawLine(Point3D point1,Point3D point2)
{
    glBegin(GL_LINES);
        glVertex3f(point1.x,point1.y,point1.z);
        glVertex3f(point2.x,point2.y,point2.z);
    glEnd();
}

//-----------------------------------------------------------------------
// readCoordinates: reads coordinated from file in GROMACS gro format
//-----------------------------------------------------------------------
void readCoordinates()
{
   std::string line;
   std::getline(inputCoordStream, line); //skip comments
   inputCoordStream >> nparticles;
   inputCoordStream.ignore(10000,'\n');
   for (int i=0;i<nparticles;i++)
   {
      std::getline(inputCoordStream, line);
      Point3D point(std::stof(line.substr(20,8)),std::stof(line.substr(28,8)),std::stof(line.substr(36,8)));
      coordinates.push_back(point);
   }
   inputCoordStream >> cell[0] >> cell[1] >> cell[2];
   inputCoordStream.ignore(10000,'\n');
}

//-----------------------------------------------------------------------
// readParticleInfo: read info about which particles to display and how
//-----------------------------------------------------------------------
void readParticleInfo()
{
   std::string line;
   std::getline(inputInfoStream, line); //skip comments
   inputInfoStream >> nparticles_vis;
   inputInfoStream.ignore(10000,'\n');
   if (nparticles_vis > nparticles) 
   {
      std::cout<<"Error!"<<nparticles_vis<<nparticles<<std::endl; //TODO
   }
   for (int i=0;i<nparticles_vis;i++)
   {
      int index;
      std::string type;
      inputInfoStream >> index;
      indices_vis.push_back(index);
      inputInfoStream >> type;
      type = type.substr(0,1);
      parttypes.push_back(type);
      sphereSizes.push_back(sizeMap[type]);
      sphereColors.push_back(colorMap[type]); 
   }
}
//-----------------------------------------------------------------------
// readBondPairs: reads the indices and distances of bonds in the ENM
//-----------------------------------------------------------------------
void readBondPairs()
{
   std::string line;
   //read bonds
   std::getline(inputBondStream, line); //skip comments
   inputBondStream >> nbonds;
   inputBondStream.ignore(10000,'\n');
   std::cout<<"Reading "<<nbonds<<" backbone bonds"<<std::endl;
   for (int i=0;i<nbonds;i++)
   {
      std::vector<int> indices;
      int temp;
      inputBondStream >> temp;
      indices.push_back(temp);
      inputBondStream >> temp;
      indices.push_back(temp);
      bondsList.push_back(indices);
      inputBondStream.ignore(10000,'\n');
   }
   //read non-bonded connections
   std::getline(inputBondStream, line); //skip comments
   std::getline(inputBondStream, line); //skip comments
   inputBondStream >> nnbs;
   inputBondStream.ignore(10000,'\n');
   std::cout<<"Reading "<<nnbs<<" non-bonded connections"<<std::endl;
   for (int i=0;i<nnbs;i++)
   {
      std::vector<int> indices;
      int temp;
      inputBondStream >> temp;
      indices.push_back(temp);
      inputBondStream >> temp;
      indices.push_back(temp);
      nonBondsList.push_back(indices);
      inputBondStream.ignore(10000,'\n');
   }

}

//-----------------------------------------------------------------------
// display callback function
//        This is called each time application needs to redraw itself.
//        Most of the opengl work is done through this function.
//-----------------------------------------------------------------------

void display()
{
    glClear(
        GL_COLOR_BUFFER_BIT |            // clear the frame buffer (color)
        GL_DEPTH_BUFFER_BIT);            // clear the depth buffer (depths)

    glPushMatrix();                        // pushes the current matrix stack down by one, duplicating the current matrix

    //glTranslated(-5,-3,-1);
    //glRotated(rotAngle, 0, 1, 0);        // rotate by rotAngle about y-axis
    //glTranslated(5,3,1);

    //gluLookAt(                            // set up the camera
    //            xpos, ypos, zpos,      // eye position
    //            3.0, 3.0, 3.0,            // lookat position
    //            0.5, 1.0, 0.0);           // up direction

    glEnable(GL_COLOR_MATERIAL);        // specify object color

    //draw particles

    for(int i=0; i<nparticles_vis; i++)
    {
        int j = indices_vis[i];
        glColor3f(sphereColors[i].color[0],sphereColors[i].color[1],sphereColors[i].color[2]); 
        double x=coordinates[j].x;
        double y=coordinates[j].y;
        double z=coordinates[j].z;
        glPushMatrix();
        glTranslated(x,y,z);
        //glScaled(0.5,0.5,0.5);
        glutSolidSphere(sphereSizes[i],25,25);
        glPopMatrix();
    }

    //draw bonds between calpha particles in ENM
    glLineWidth(4.0);

    glColor3f(0.0, 0.0, 1.0);        // blue
    for(int i=0; i<nbonds; i++)
    {
        drawLine(coordinates[bondsList[i][0]-1],coordinates[bondsList[i][1]-1]);
    }

    glColor3f(0.0, 1.0, 0.0);        // green
    for(int i=0; i<nnbs; i++)
    {
        drawLine(coordinates[nonBondsList[i][0]-1],coordinates[nonBondsList[i][1]-1]);
    }

    //draw bonds between atomistic particles
    //if distance<bondCutoff

    glPopMatrix();                        // restore the modelview matrix
    glFlush();                            // force OpenGL to render now

    glutSwapBuffers();                    // make the image visible
}

//-----------------------------------------------------------------------
// keyboard callback function
//        This is called whenever a keyboard key is hit.
//-----------------------------------------------------------------------

void keyboard(unsigned char k, int x, int y)
{
    switch (k)
    {
    case 'a':
        rotAngle += 5;                    // increase rotation by 5 degrees
        break;
    case 'l':
        rotAngle -= 5;                    // decrease rotation by 5 degrees
        break;
    case 'x':
        xpos += 1.0;
        break;
    case 'y':
        ypos += 1.0;
        break;
    case 'z':
        zpos -= 1.0;
        break;
    case 'b':
        xpos -= 1.0;
        break;
    case 'c':
        ypos -= 1.0;
        break;
    case 'd':
        zpos += 1.0;
        break;
    case 'q':
        exit(0);                        // exit
    }

    glutPostRedisplay();                // redraw the image now
}

//-----------------------------------------------------------------------
//    usage
//-----------------------------------------------------------------------

void usage()
{
    std::cout << "\n\
-----------------------------------------------------------------------\n\
  CMSC 427 Sample Program.\n\
  Usage: ./enm-calpha-cbeta grofile bondlistfile\n\
  grofile is a standard GROMACS-format coordinate file\n\
  bondlistfile contains list of bonds in ENM (for format see README)\n\
  Inputs:\n\
    a:                Rotate counterclockwise\n\
    l:                Rotate clockwise\n\
    q:                Quit\n\
  You may need to place the cursor over the graphics window for\n\
  keyboard input to be processed.\n\
-----------------------------------------------------------------------\n";
    std::cout.flush();
}

//-----------------------------------------------------------------------
// main program
//        Where everything begins.
//-----------------------------------------------------------------------

int main(int argc, char **argv)
{
    usage();                             // explain how to use
    glutInit(&argc, argv);
    glutInitDisplayMode(                 // initialize GLUT
                GLUT_DOUBLE |            // use double buffering
                GLUT_DEPTH |             // request memory for z-buffer
                GLUT_RGB );              // set RGB color mode

    glutInitWindowSize(1640,1480);
    glutCreateWindow("Elastic Network Model");    

    //read coordinates of all particles
    char* coordFileName;
    coordFileName = argv[1];
    inputCoordStream.open(coordFileName,std::ifstream::in); 
    readCoordinates();
    inputCoordStream.close();

    //read info about which particles to display and how
    char* infoFileName;
    infoFileName = argv[2];
    inputInfoStream.open(infoFileName,std::ifstream::in); 
    readParticleInfo();
    inputInfoStream.close();

    //read bond pairs
    char* bondFileName;
    bondFileName = argv[3];
    inputBondStream.open(bondFileName,std::ifstream::in); 
    readBondPairs();
    inputBondStream.close();

    glutDisplayFunc(display);            // call display() to redraw window
    glutKeyboardFunc(keyboard);          // call keyboard() when key is hit

    init();                     

    glutMainLoop();              
    return 0; 
}
