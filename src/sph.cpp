#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string> 

//include header file for glfw library so that we can use OpenGL
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "lodepng.h"


#include "sim.cpp"
#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define PI 3.14159265 // Should be used from mathlib

using namespace std;

//****************************************************
// Global Variables
//****************************************************
GLfloat translation[3] = {0.0f, 0.0f, 0.0f};
bool auto_strech = false;
int Width_global = 800;
int Height_global = 800;
int Z_buffer_bit_depth = 128;

float cameraLeft = -5;
float cameraRight = 25;
float cameraBottom = -5;
float cameraTop = 25;
float lastRun = 0;

float ratioWin = 1;

Simulator* sim;

inline float sqr(float x) { return x*x; }


//****************************************************
// Simple init function
//****************************************************
void initializeRendering()
{
    glfwInit();
}


//****************************************************
// Keyboard inputs. Add things to match the spec! 
//****************************************************
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    switch (key) {
            
        case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
        case GLFW_KEY_Q: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
        case GLFW_KEY_LEFT :
            if (action && mods == GLFW_MOD_SHIFT) translation[0] -= 0.001f * Width_global; break;
        case GLFW_KEY_RIGHT:
            if (action && mods == GLFW_MOD_SHIFT) translation[0] += 0.001f * Width_global; break;
        case GLFW_KEY_UP   :
            if (action && mods == GLFW_MOD_SHIFT) translation[1] += 0.001f * Height_global; break;
        case GLFW_KEY_DOWN :
            if (action && mods == GLFW_MOD_SHIFT) translation[1] -= 0.001f * Height_global; break;
        case GLFW_KEY_F:
            if (action && mods == GLFW_MOD_SHIFT) auto_strech = !auto_strech; break;
        case GLFW_KEY_SPACE: break;
        case GLFW_KEY_EQUAL:
            if (action) {
                cameraLeft += 0.1;
                cameraRight -= 0.1;
                cameraBottom += 0.1;
                cameraTop -= 0.1;
            }
            break;
        case GLFW_KEY_MINUS:
            if (action) {
                cameraLeft -= 0.1;
                cameraRight += 0.1;
                cameraBottom -= 0.1;
                cameraTop += 0.1;
            }
            break;
        case GLFW_KEY_S:
            if (action) {
                sim->stepSPH();
            }
            break;
        default: break;
    }
    
}
const float DEG2RAD = 3.14159/180;
 
void drawCircle(float radius, float x, float y, float z) {
    glBegin(GL_TRIANGLE_FAN);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glVertex3f(x, y, z);
    for (int i=0; i <= 360; i++) {
        float degInRad = i*DEG2RAD;
        glVertex3f(cos(degInRad)*radius + x,sin(degInRad)*radius + y, z);
    }

    glEnd();
}
void drawSim() {
    vector<Particle*> particles = sim->particles;
    for (int i = 0; i < particles.size(); i++) {
        if (particles[i]->ghost == false) {
            drawCircle(1, particles[i]->pos.x, particles[i]->pos.z, particles[i]->pos.y);
        } /*else {
            glPointSize(10.0f);
            glBegin(GL_POINTS);
            glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
            glVertex3f(particles[i]->pos.x, particles[i]->pos.z, particles[i]->pos.y);
            glEnd();
        }*/
    }
}

int iteration = 0;
//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void display( GLFWwindow* window )
{
    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f ); //clear background screen to black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                // clear the color buffer (sets everything to black)

    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(cameraLeft * ratioWin, cameraRight * ratioWin, cameraBottom, cameraTop, 5, -5);
    glMatrixMode(GL_MODELVIEW);                  // indicate we are specifying camera transformations
    glLoadIdentity();                            // make sure transformation is "zero'd"

    //----------------------- code to draw objects --------------------------
    glPushMatrix();
    glTranslatef(translation[0], translation[1], translation[2]);
    if (iteration < 400) {
        string filename = "images/" + to_string(iteration) + ".png";
        //printf("%f\n", delta);
        sim->stepSPH(0.1);
        drawSim();
        int width = Width_global;
        int height = Height_global;
        int w = width;
        int h = height;
        vector<unsigned char> data(width*height*4);
        glReadPixels(0,0,width,height,GL_BGR,GL_UNSIGNED_BYTE,&data[0]);
        vector<unsigned char> imagedata;
        imagedata.resize(w * h * 4);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < w; x++) {
                int index = y * width * 3 + x * 3;
                imagedata[4 * w * (height-y - 1) + 4 * x + 0] = data[index+2];
                imagedata[4 * w * (height-y - 1) + 4 * x + 1] = data[index+1];
                imagedata[4 * w * (height-y - 1) + 4 * x + 2] = data[index];
                imagedata[4 * w * (height-y - 1) + 4 * x + 3] = 255;
            }
        }
        unsigned error = lodepng::encode(filename.c_str(), imagedata, width, height);
        if (error) {
            printf("There was an error saving to %s.\n", filename.c_str());
        } else {
            printf("Wrote image to %s!\n", filename.c_str());
        }
        iteration++;
    } else {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    glPopMatrix();
    
    glfwSwapBuffers(window);

    // note: check out glPolygonMode and glShadeModel 
    // for wireframe and shading commands
    
}

//****************************************************
// function that is called when window is resized
//***************************************************
void size_callback(GLFWwindow* window, int width, int height)
{
    // Get the pixel coordinate of the window
    // it returns the size, in pixels, of the framebuffer of the specified window
    glfwGetFramebufferSize(window, &Width_global, &Height_global);
    
    glViewport(0, 0, Width_global, Height_global);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    ratioWin = Width_global * 1.0f/Height_global;
    glOrtho(cameraLeft * ratioWin, cameraRight * ratioWin, cameraBottom, cameraTop, 5, -5);
    
    display(window);
}


//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {

    //Simulator Setup
    sim = new Simulator(60, 1, 60, 20, 1, 20);

    //This initializes glfw
    initializeRendering();
    
    GLFWwindow* window = glfwCreateWindow( Width_global, Height_global, "Water sim", NULL, NULL );
    if ( !window )
    {
        cerr << "Error on window creating" << endl;
        glfwTerminate();
        return -1;
    }

    const GLFWvidmode * mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    if ( !mode )
    {
        cerr << "Error on getting monitor" << endl;
        glfwTerminate();
        return -1;
    }
    
    glfwMakeContextCurrent( window );
    
    // Get the pixel coordinate of the window
    // it returns the size, in pixels, of the framebuffer of the specified window
    glfwGetFramebufferSize(window, &Width_global, &Height_global);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	glOrtho(cameraLeft, cameraRight, cameraBottom, cameraTop, 5, -5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

    glEnable(GL_DEPTH_TEST);	// enable z-buffering
    glDepthFunc(GL_LESS);

    glfwSetWindowTitle(window, "CS184");
    glfwSetWindowSizeCallback(window, size_callback);
    glfwSetKeyCallback(window, key_callback);

    while( !glfwWindowShouldClose( window ) ) // infinite loop to draw object again and again
    {   // because once object is draw then window is terminated
        display( window );
        
        if (auto_strech){
            glfwSetWindowSize(window, mode->width, mode->height);
            glfwSetWindowPos(window, 0, 0);
        }
        
        glfwPollEvents();
        
    }

    return 0;
}