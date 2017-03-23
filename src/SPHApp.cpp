#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <memory>
#include <chrono>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <omp.h>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Particle.hpp"
#include "HashTable.hpp"
#include "Shader.hpp"

#include "Settings.hpp"
#include "Kernels.hpp"
#include "BoundingBox.hpp"
#include "SPHSolver.hpp"
#include "WriteFunctions.hpp"


////void generateParticles();
void init();

Shader particleShader;
Shader forceShader;
GLuint VAO,VBO_p,VBO_c,VBO_f;
glm::mat4 cameraView;

//ScreenPoints sCloud;
glm::vec2* sPosition;
glm::vec3* sColor;
glm::vec2* sForce;

//glm::vec3* fColor;

SPHSolver SPHsolver;

float calcScale(glm::vec2 Fmax)
{
    static float pixelWidth = RenderSettings::displayBox.dx()/float(RenderSettings::width);
    float scale = 50.0f*pixelWidth/float(glm::length(Fmax)+1.e-5);
    return scale;
}

glm::dvec2 calcMaxForce(glm::vec2* force)
{
    glm::vec2 Fmax(0.0);
    float     len2Fmax=0.0;
    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        if (glm::length2(force[iPart])>len2Fmax)
        {
            Fmax = force[iPart];
            len2Fmax = glm::length2(Fmax);
        }
    }
    return Fmax;
}

//void keyboard(unsigned char ch, int x, int y)
void keyboard(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    switch(key)
    {
    case 27:
    case 'q':
    case 'Q':
        if (action == GLFW_PRESS) glfwSetWindowShouldClose(window, GL_TRUE);
        break;
    case ' ':
    {
        if (action == GLFW_PRESS)
        {
            if (RenderSettings::displayRender==RenderSettings::SIMPLE)
            {
                RenderSettings::displayRender=RenderSettings::PRESSFORCES;
                glfwSetWindowTitle(window, "Visualization of pressure forces");
            }
            else if (RenderSettings::displayRender==RenderSettings::PRESSFORCES)
            {
                RenderSettings::displayRender=RenderSettings::VISCFORCES;
                glfwSetWindowTitle(window, "Visualization of viscous forces");
            }
            else if (RenderSettings::displayRender==RenderSettings::VISCFORCES)
            {
                RenderSettings::displayRender=RenderSettings::OTHERFORCES;
                glfwSetWindowTitle(window, "Visualization of other forces");
            }
            else if (RenderSettings::displayRender==RenderSettings::OTHERFORCES)
            {
                RenderSettings::displayRender=RenderSettings::ALLFORCES;
                glfwSetWindowTitle(window, "Visualization of total forces");
            }
            else if (RenderSettings::displayRender==RenderSettings::ALLFORCES)
            {
                RenderSettings::displayRender=RenderSettings::SIMPLE;
                glfwSetWindowTitle(window, "Force visualization - disabled");
            }
        }
        break;
    }
    default:{}
    }
}

void renderParticles(GLFWwindow* window)
{
    //for (int i=0;i<stepsPerFrame;i++)
    static std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    static std::chrono::time_point<std::chrono::system_clock> now   = std::chrono::system_clock::now();
    static std::chrono::duration<double> deltaTime =std::chrono::duration<double>(0.0);

    now = std::chrono::system_clock::now();
    deltaTime = now - start;

    while (deltaTime.count()<16.5e-3) 
    {
        SPHsolver.step(window);
        now = std::chrono::system_clock::now();
        deltaTime = now - start;
    }
    start = std::chrono::system_clock::now();

    //update screen positions
    std::vector<Particle>& cloud = SPHsolver.cloud();
    #pragma omp parallel for
    for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
    {
        if (!cloud[iPart].active) continue;
        Particle& iParticle = cloud[iPart];
        //sCloud.position[iPart] = iParticle.position;
        sPosition[iPart] = iParticle.position;
        float velMag = float(glm::length(iParticle.velocity));
        float scale = 0.2f;
        //sCloud.color[iPart] = glm::vec3(
        sColor[iPart] = glm::vec3
        (
            scale*velMag,
            scale*velMag, 
            scale*velMag*0.3f+0.7
        );
    }

    particleShader.use();
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_p);
    glBufferSubData(GL_ARRAY_BUFFER, 0, SPHSettings::NParticles*sizeof(glm::vec2), sPosition);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_c);
    glBufferSubData(GL_ARRAY_BUFFER, 0, SPHSettings::NParticles*sizeof(glm::vec3), sColor);

    float particleSize = 0.75f*float(BoundaryConditions::bndBox.dx()*SPHSettings::initDx)*float(RenderSettings::width);
    GLint cameraViewLoc = glGetUniformLocation(particleShader.program(), "cameraView");
    glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView));
    GLint pointSizeLoc = glGetUniformLocation(particleShader.program(), "pointSize");
    glUniform1f(pointSizeLoc, particleSize);
    glDrawArrays(GL_POINTS, 0, SPHSettings::NParticles);
    glBindVertexArray(0);

    if (RenderSettings::displayRender!=RenderSettings::SIMPLE)
    {
        forceShader.use();

        cameraViewLoc = glGetUniformLocation(forceShader.program(), "cameraView");
        glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView));

        GLint colorLoc = glGetUniformLocation(forceShader.program(), "Color");

        if (RenderSettings::displayRender==RenderSettings::PRESSFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
            {
                if (!cloud[iPart].active) continue;
                Particle& iParticle = cloud[iPart];
                sForce[iPart] = iParticle.Fpress;

            }

            glUniform3f(colorLoc, 1.0f, 1.0f, 1.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::VISCFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
            {
                if (!cloud[iPart].active) continue;
                Particle& iParticle = cloud[iPart];
                sForce[iPart] = iParticle.Fvisc;
            }

            glUniform3f(colorLoc, 0.0f, 1.0f, 1.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::OTHERFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
            {
                if (!cloud[iPart].active) continue;
                Particle& iParticle = cloud[iPart];
                sForce[iPart] = iParticle.Fother;
            }

            glUniform3f(colorLoc, 1.0f, 0.0f, 1.0f);
        }

        else if (RenderSettings::displayRender==RenderSettings::ALLFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
            {
                if (!cloud[iPart].active) continue;
                Particle& iParticle = cloud[iPart];
                sForce[iPart] = iParticle.Ftot;
            }

            glUniform3f(colorLoc, 1.0f, 1.0f, 0.0f);
        }

        glm::vec2 Fmax=calcMaxForce(sForce);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO_f);
        glBufferSubData(GL_ARRAY_BUFFER, 0, SPHSettings::NParticles*sizeof(glm::vec2), sForce);

        GLint scaleLoc = glGetUniformLocation(forceShader.program(), "scale");
        glUniform1f(scaleLoc, calcScale(Fmax));
        glDrawArrays(GL_POINTS, 0, SPHSettings::NParticles);
        glBindVertexArray(0);
    }
    else if (RenderSettings::displayRender==RenderSettings::OTHERFORCES)
    {
        #pragma omp parallel for
        for (int iPart = 0;iPart<SPHSettings::NParticles;iPart++) 
        {
            if (!cloud[iPart].active) continue;
            Particle& iParticle = cloud[iPart];
            sForce[iPart] = iParticle.Fother;
        }

        forceShader.use();
        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO_f);
        glBufferSubData(GL_ARRAY_BUFFER, 0, SPHSettings::NParticles*sizeof(glm::vec2), sForce);

        cameraViewLoc = glGetUniformLocation(forceShader.program(), "cameraView");
        glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView));
        GLint colorLoc = glGetUniformLocation(forceShader.program(), "Color");
        glUniform3f(colorLoc, 1.0f, 0.0f, 0.0f);
        GLint scaleLoc = glGetUniformLocation(forceShader.program(), "scale");
        glUniform1f(scaleLoc, 0.00001f);
        glDrawArrays(GL_POINTS, 0, SPHSettings::NParticles);
        glBindVertexArray(0);
    }
}

int main( int argc, char** argv )
{
    init();

	// Init GLFW
	glfwInit();
	// Set all the required options for GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

	// retrive primary monitor
	//GLFWmonitor* primary = glfwGetPrimaryMonitor();

	// Create a GLFWwindow object that we can use for GLFW's functions
    std::cout<<"Creating window"<<std::endl;
	GLFWwindow* window = glfwCreateWindow(RenderSettings::width, RenderSettings::height, "SPH", nullptr, nullptr);    
	if (window == nullptr)
	{
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		exit(1);
	}
	glfwMakeContextCurrent(window);
	// Set the required callback functions
    std::cout<<"Setting callback functions"<<std::endl;
	//glfwSetMouseButtonCallback(window, mouseButton_callback);
	glfwSetKeyCallback(window, keyboard);

	// Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
    std::cout<<"Initializing GLEW"<<std::endl;
    glewExperimental = GL_TRUE;
    // Initialize GLEW to setup the OpenGL Function pointers
    if (glewInit() != GLEW_OK)
    {
    	std::cerr << "Failed to initialize GLEW" << std::endl;
    	exit(1);
    }    

    //create shader programs
    particleShader.compileVertexShader("src/shaders/pointShader.vs");
    particleShader.compileFragmentShader("src/shaders/pointShader.fs");
    particleShader.linkShaders();

    forceShader.compileVertexShader("src/shaders/forceShader.vs");
    forceShader.compileGeometryShader("src/shaders/forceShader.gs");
    forceShader.compileFragmentShader("src/shaders/forceShader.fs");
    forceShader.linkShaders();

    // Create Vertex Buffer and Array Objects
    // Create Vertex Buffer
    glGenBuffers(1, &VBO_p);
    glGenBuffers(1, &VBO_c);
    glGenBuffers(1, &VBO_f);
    // Create Array for triangle
    glGenVertexArrays(1, &VAO);

    // Bind Vertex Array Object
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_p);
    glBufferData(GL_ARRAY_BUFFER, SPHSettings::NParticles*(sizeof(glm::vec2)), sPosition, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_c);
    glBufferData(GL_ARRAY_BUFFER, SPHSettings::NParticles*(sizeof(glm::vec3)), sColor, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (GLvoid*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_f);
    glBufferData(GL_ARRAY_BUFFER, SPHSettings::NParticles*(sizeof(glm::vec2)), sForce, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (GLvoid*)0);
    glEnableVertexAttribArray(2);

	// Define the viewport dimensions
	glViewport(0, 0, RenderSettings::width, RenderSettings::height);

    glEnable(GL_PROGRAM_POINT_SIZE);

    // set clear color
    glClearColor( 0.1f, 0.1f, 0.1f, 1.0f );

    // enable vsync
    glfwSwapInterval(1);

    //define camera matrix
    glm::mat4 viewMatrix = glm::lookAt(glm::vec3(0.0f,0.0f, 1.0f),
                                       glm::vec3(0.0f,0.0f, 0.0f),
                                       glm::vec3(0.0f,1.0f, 0.0f));
    glm::mat4 projMatrix = glm::ortho( RenderSettings::displayBox.minX(), 
                                       RenderSettings::displayBox.maxX(), 
                                       RenderSettings::displayBox.minY(), 
                                       RenderSettings::displayBox.maxY(), 
                                       0.f, 10.f);
    cameraView = projMatrix*viewMatrix;

    std::cout<<"starting game loop "<<std::endl;
	// Game loop
	while (!glfwWindowShouldClose(window))
	{
		// Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
		glfwPollEvents();

		// Render
		// Clear the colorbuffer
        glClear( GL_COLOR_BUFFER_BIT );

        //render image
	    renderParticles(window);

	    glfwSwapBuffers(window);
	}

    glDeleteBuffers(1, &VBO_p );
    glDeleteBuffers(1, &VBO_c );
    glDeleteBuffers(1, &VBO_f );
    glDeleteVertexArrays(1, &VAO );

	std::cout<< "destroying window" <<std::endl;
	glfwDestroyWindow(window);

	std::cout<< "terminating glfw" <<std::endl;
	glfwTerminate();

	std::cout<< "deallocating arrays" <<std::endl;
    //delete [] sCloud;
      delete [] sPosition;
      delete [] sColor;
      delete [] sForce;
    //delete [] fColor;
}


void init()
{
    Kernel::SmoothingLength::setSmoothingLength(SPHSettings::initDx*2.);

    SPHsolver.init();

  //sCloud.position.resize(SPHSettings::NParticles);
  //sCloud.color.resize(SPHSettings::NParticles);
  //sCloud.force.resize(SPHSettings::NParticles);

    sPosition = new glm::vec2[SPHSettings::NParticles];
    sColor    = new glm::vec3[SPHSettings::NParticles];
    sForce    = new glm::vec2[SPHSettings::NParticles];

    std::vector<Particle>& cloud = SPHsolver.cloud();
    for (int i=0;i<SPHSettings::NParticles;i++)
    {
        //sCloud.position[i] = cloud[i].position;
        sPosition[i] = cloud[i].position;
        sColor[i]    = glm::vec3(0.0f);
        sForce[i]    = glm::vec2(0.0f);
    }
}

