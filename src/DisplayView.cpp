
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "DisplayView.hpp"
#include "Settings.hpp"

//********************************************************************************
void DisplayView::keyboard(GLFWwindow* window, int key, int scancode, int action, int mode)
//********************************************************************************
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
                RenderSettings::displayRender=RenderSettings::SURFFORCES;
                glfwSetWindowTitle(window, "Visualization of surface forces");
            }
            else if (RenderSettings::displayRender==RenderSettings::SURFFORCES)
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

//********************************************************************************
void DisplayView::WindowManager::init(std::vector<Particle>& cloud)
//********************************************************************************
{

    sPosition_.resize(SPHSettings::NParticles);
    sColor_.resize(SPHSettings::NParticles);
    sForce_.resize(SPHSettings::NParticles);

    for (int i=0;i<SPHSettings::NParticles;i++)
    {
        //sCloud.position[i] = cloud[i].position;
        sPosition_[i] = cloud[i].position;
        sColor_[i]    = glm::vec3(1.0f);
        sForce_[i]    = glm::vec2(0.0f);
    }

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
	window_ = glfwCreateWindow(RenderSettings::width, RenderSettings::height, "SPH", nullptr, nullptr);    
	if (window_ == nullptr)
	{
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		exit(1);
	}
	glfwMakeContextCurrent(window_);
	// Set the required callback functions
    std::cout<<"Setting callback functions"<<std::endl;
	//glfwSetMouseButtonCallback(window, mouseButton_callback);
	glfwSetKeyCallback(window_, DisplayView::keyboard);

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
    particleShader_.compileShaderPart("src/shaders/pointShader.vs",GL_VERTEX_SHADER);
    particleShader_.compileShaderPart("src/shaders/pointShader.fs",GL_FRAGMENT_SHADER);
    particleShader_.linkProgram();

    forceShader_.compileShaderPart("src/shaders/forceShader.vs",GL_VERTEX_SHADER);
    forceShader_.compileShaderPart("src/shaders/forceShader.gs",GL_GEOMETRY_SHADER);
    forceShader_.compileShaderPart("src/shaders/forceShader.fs",GL_FRAGMENT_SHADER);
    forceShader_.linkProgram();

    // Create Vertex Buffer and Array Objects
    // Create Vertex Buffer
    glGenBuffers(1, &VBO_p_);
    glGenBuffers(1, &VBO_c_);
    glGenBuffers(1, &VBO_f_);
    // Create Array for triangle
    glGenVertexArrays(1, &VAO_);

    // Bind Vertex Array Object
    glBindVertexArray(VAO_);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_p_);
    glBufferData(GL_ARRAY_BUFFER, SPHSettings::NParticles*(sizeof(glm::vec2)), &sPosition_[0], GL_STREAM_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_c_);
    glBufferData(GL_ARRAY_BUFFER, SPHSettings::NParticles*(sizeof(glm::vec3)), &sColor_[0], GL_STREAM_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (GLvoid*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_f_);
    glBufferData(GL_ARRAY_BUFFER, SPHSettings::NParticles*(sizeof(glm::vec2)), &sForce_[0], GL_STREAM_DRAW);
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
    cameraView_ = projMatrix*viewMatrix;
}

//********************************************************************************
bool DisplayView::WindowManager::renderParticles(std::vector<Particle>& cloud, const int NParticles)
//********************************************************************************
{
    // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
    glfwPollEvents();

    // Clear the colorbuffer
    glClear( GL_COLOR_BUFFER_BIT );
    
    //update screen positions
    #pragma omp parallel for
    for (int iPart = 0;iPart<NParticles;iPart++) 
    {
        sPosition_[iPart] = cloud[iPart].position;

        float velMag = float(glm::length(cloud[iPart].velocity));
        float scale = 0.2f;
        //blue to white
        sColor_[iPart] = glm::vec3
        (
            scale*velMag,
            scale*velMag, 
            scale*velMag*0.3f+0.7
        );
    }

    particleShader_.use();
    glBindVertexArray(VAO_);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_p_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, NParticles*sizeof(glm::vec2), &sPosition_[0]);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_c_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, NParticles*sizeof(glm::vec3), &sColor_[0]);

    float particleSize = 0.75f*float(BoundaryConditions::bndBox.dx()*SPHSettings::initDx)*float(RenderSettings::width);

    GLint cameraViewLoc = glGetUniformLocation(particleShader_.program(), "cameraView");
    glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView_));

    GLint pointSizeLoc = glGetUniformLocation(particleShader_.program(), "pointSize");
    glUniform1f(pointSizeLoc, particleSize);

    glDrawArrays(GL_POINTS, 0, SPHSettings::NParticles);
    glBindVertexArray(0);

    if (RenderSettings::displayRender!=RenderSettings::SIMPLE)
    {
        forceShader_.use();

        cameraViewLoc = glGetUniformLocation(forceShader_.program(), "cameraView");
        glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView_));

        GLint colorLoc = glGetUniformLocation(forceShader_.program(), "Color");

        if (RenderSettings::displayRender==RenderSettings::PRESSFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fpress;
            }

            glUniform3f(colorLoc, 1.0f, 1.0f, 1.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::VISCFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fvisc;
            }

            glUniform3f(colorLoc, 0.0f, 1.0f, 1.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::SURFFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fsurf;
              //sForce_[iPart] = cloud[iPart].normal;
            }

            glUniform3f(colorLoc, 1.0f, 0.0f, 1.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::OTHERFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fother;
            }

            glUniform3f(colorLoc, 1.0f, 0.0f, 1.0f);
        }

        else if (RenderSettings::displayRender==RenderSettings::ALLFORCES)
        {
            #pragma omp parallel for
            for (int iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Ftot;
            }

            glUniform3f(colorLoc, 1.0f, 1.0f, 0.0f);
        }

        glm::vec2 Fmax=calcMaxForce(sForce_);

        glBindVertexArray(VAO_);

        glBindBuffer(GL_ARRAY_BUFFER, VBO_f_);
        glBufferSubData(GL_ARRAY_BUFFER, 0, NParticles*sizeof(glm::vec2), &sForce_[0]);

        GLint scaleLoc = glGetUniformLocation(forceShader_.program(), "scale");
        glUniform1f(scaleLoc, calcScale(Fmax));
        glDrawArrays(GL_POINTS, 0, SPHSettings::NParticles);
        glBindVertexArray(0);
    }
    
    glfwSwapBuffers(window_);
    return !glfwWindowShouldClose(window_);
}

//********************************************************************************
float DisplayView::WindowManager::calcScale(glm::vec2 Fmax)
//********************************************************************************
{
    static float pixelWidth = RenderSettings::displayBox.dx()/float(RenderSettings::width);
    float scale = 50.0f*pixelWidth/float(glm::length(Fmax)+1.e-5);
    return scale;
}

//********************************************************************************
glm::dvec2 DisplayView::WindowManager::calcMaxForce(std::vector<glm::vec2>& force)
//********************************************************************************
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
