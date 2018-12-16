
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
#include "HardCodedShaders.hpp"
#include "Simulation/Settings.hpp"

#define GLLogging 0

//********************************************************************************
void DisplayView::keyboard(GLFWwindow* window, int key, int /*scancode*/, int action, int /*mode*/)
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
            if (RenderSettings::displayRender==RenderSettings::INDEX)
            {
                RenderSettings::displayRender=RenderSettings::SIMPLE;
                glfwSetWindowTitle(window, "Velocity-based Visualization");
            }
            else if (RenderSettings::displayRender==RenderSettings::SIMPLE)
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
                RenderSettings::displayRender=RenderSettings::INDEX;
                glfwSetWindowTitle(window, "Index-based visualization");
            }
        }
        break;
    }
    default:{}
    }
}
//
//********************************************************************************
DisplayView::WindowManager::~WindowManager()
//********************************************************************************
{
    glDeleteBuffers(1, &VBO_p_ );
    glDeleteBuffers(1, &VBO_c_ );
    glDeleteBuffers(1, &VBO_f_ );
    glDeleteVertexArrays(1, &VAO_ );
    
#if GLLogging
    std::cout<< "destroying window" <<std::endl;
#endif
    glfwDestroyWindow(window_);
    
#if GLLogging
    std::cout<< "terminating glfw" <<std::endl;
#endif
    glfwTerminate();
}


//********************************************************************************
void DisplayView::WindowManager::init()
//********************************************************************************
{

    sPosition_.resize(SPHSettings::NParticles);
    sColor_.resize(SPHSettings::NParticles);
    sForce_.resize(SPHSettings::NParticles);

    for (unsigned i=0;i<SPHSettings::NParticles;i++)
    {
        sPosition_[i] = glm::vec2(0.);
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
#if GLLogging
    std::cout<<"Creating window"<<std::endl;
#endif
	window_ = glfwCreateWindow(int(RenderSettings::width), int(RenderSettings::height), "SPH", nullptr, nullptr);    
	if (window_ == nullptr)
	{
		std::cerr << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		exit(1);
	}
	glfwMakeContextCurrent(window_);
	// Set the required callback functions
#if GLLogging
    std::cout<<"Setting callback functions"<<std::endl;
#endif
	//glfwSetMouseButtonCallback(window, mouseButton_callback);
	glfwSetKeyCallback(window_, DisplayView::keyboard);

	// Set this to true so GLEW knows to use a modern approach to retrieving function pointers and extensions
#if GLLogging
    std::cout<<"Initializing GLEW"<<std::endl;
#endif
    glewExperimental = GL_TRUE;
    // Initialize GLEW to setup the OpenGL Function pointers
    if (glewInit() != GLEW_OK)
    {
    	std::cerr << "Failed to initialize GLEW" << std::endl;
    	exit(1);
    }    

    //create shader programs
    particleShader_.compileShaderPart(HardCodedShader::pointShader_vs,GL_VERTEX_SHADER);
    particleShader_.compileShaderPart(HardCodedShader::pointShader_fs,GL_FRAGMENT_SHADER);
    particleShader_.linkProgram();

    forceShader_.compileShaderPart(HardCodedShader::forceShader_vs,GL_VERTEX_SHADER);
    forceShader_.compileShaderPart(HardCodedShader::forceShader_gs,GL_GEOMETRY_SHADER);
    forceShader_.compileShaderPart(HardCodedShader::forceShader_fs,GL_FRAGMENT_SHADER);
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
    glBufferData(GL_ARRAY_BUFFER, sPosition_.size()*(sizeof(glm::vec2)), sPosition_.data(), GL_STREAM_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_c_);
    glBufferData(GL_ARRAY_BUFFER, sColor_.size()*(sizeof(glm::vec3)), sColor_.data(), GL_STREAM_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (GLvoid*)0);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_f_);
    glBufferData(GL_ARRAY_BUFFER, sForce_.size()*(sizeof(glm::vec2)), sForce_.data(), GL_STREAM_DRAW);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (GLvoid*)0);
    glEnableVertexAttribArray(2);

	// Define the viewport dimensions
	glViewport(0, 0, int(RenderSettings::width), int(RenderSettings::height));

    glEnable(GL_PROGRAM_POINT_SIZE);

    // set clear color
    glClearColor( 0.2f, 0.2f, 0.0f, 1.0f );

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
bool DisplayView::WindowManager::renderParticles(const ParticleCloud& cloud)
//********************************************************************************
{
    // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
    glfwPollEvents();

    // Clear the colorbuffer
    glClear( GL_COLOR_BUFFER_BIT );
    
    const auto& particlePos = cloud.get<Attr::ePosition>();
    const auto& particleVel = cloud.get<Attr::eVelocity>();

    //update screen positions
    const size_t NParticles = cloud.size();
    #pragma omp parallel for
    for (size_t iPart = 0;iPart<NParticles;iPart++) 
    {
        sPosition_[iPart] = particlePos[iPart];

        if (RenderSettings::displayRender==RenderSettings::INDEX)
        {
            float clr = float(iPart)/float(NParticles-1);
            sColor_[iPart] = glm::vec3
            (
                clr,
                clr,
                clr
            );
        }
        else if (RenderSettings::displayRender==RenderSettings::SIMPLE)
        {
            float velMag = float(glm::length(particleVel[iPart]));
            float scale = 0.2f;
            //blue to white
            sColor_[iPart] = glm::vec3
            (
                scale*velMag,
                scale*velMag, 
                scale*velMag*0.3f+0.7
            );
        }
        else 
        {
            sColor_[iPart] = glm::vec3(0.0f);
        }
    }

    particleShader_.use();
    glBindVertexArray(VAO_);

    glBindBuffer(GL_ARRAY_BUFFER, VBO_p_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, NParticles*sizeof(glm::vec2), sPosition_.data());
    glBindBuffer(GL_ARRAY_BUFFER, VBO_c_);
    glBufferSubData(GL_ARRAY_BUFFER, 0, NParticles*sizeof(glm::vec3), sColor_.data());

    float particleSize = float((SPHSettings::initDx*double(RenderSettings::width)))/RenderSettings::displayBox.dx();

    GLint cameraViewLoc = glGetUniformLocation(particleShader_.program(), "cameraView");
    glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView_));

    GLint pointSizeLoc = glGetUniformLocation(particleShader_.program(), "pointSize");
    glUniform1f(pointSizeLoc, particleSize);

    glDrawArrays(GL_POINTS, 0, int(NParticles));
    glBindVertexArray(0);

    if (RenderSettings::displayRender>RenderSettings::SIMPLE)
    {
        forceShader_.use();

        cameraViewLoc = glGetUniformLocation(forceShader_.program(), "cameraView");
        glUniformMatrix4fv(cameraViewLoc, 1, GL_FALSE, glm::value_ptr(cameraView_));

        GLint colorLoc = glGetUniformLocation(forceShader_.program(), "Color");

        if (RenderSettings::displayRender==RenderSettings::PRESSFORCES)
        {
            const auto& particlePress = cloud.get<Attr::ePressForce>();
            #pragma omp parallel for
            for (size_t iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = particlePress[iPart];
            }

            glUniform3f(colorLoc, 1.0f, 0.0f, 0.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::VISCFORCES)
        {
            #pragma omp parallel for
            for (size_t iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fvisc;
            }

            glUniform3f(colorLoc, 0.0f, 1.0f, 0.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::SURFFORCES)
        {
            #pragma omp parallel for
            for (size_t iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fsurf;
            }

            glUniform3f(colorLoc, 0.0f, 0.0f, 1.0f);
        }
        else if (RenderSettings::displayRender==RenderSettings::OTHERFORCES)
        {
            #pragma omp parallel for
            for (size_t iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Fother;
            }

            glUniform3f(colorLoc, 1.0f, 1.0f, 0.0f);
        }

        else if (RenderSettings::displayRender==RenderSettings::ALLFORCES)
        {
            #pragma omp parallel for
            for (size_t iPart = 0;iPart<NParticles;iPart++) 
            {
                sForce_[iPart] = cloud[iPart].Ftot;
            }

            glUniform3f(colorLoc, 1.0f, 1.0f, 1.0f);
        }

        glm::vec2 Fmax=calcMaxForce(sForce_);

        glBindVertexArray(VAO_);

        glBindBuffer(GL_ARRAY_BUFFER, VBO_f_);
        glBufferSubData(GL_ARRAY_BUFFER, 0, NParticles*sizeof(glm::vec2), sForce_.data());

        GLint scaleLoc = glGetUniformLocation(forceShader_.program(), "scale");
        glUniform1f(scaleLoc, calcScale(Fmax));
        glDrawArrays(GL_POINTS, 0, int(NParticles));
        glBindVertexArray(0);
    }
    
    glfwSwapBuffers(window_);
    renderTimer_.start();
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
    for (const auto& vec : force) 
    {
        if (glm::length2(vec)>len2Fmax)
        {
            Fmax = vec;
            len2Fmax = glm::length2(Fmax);
        }
    }
    return Fmax;
}
