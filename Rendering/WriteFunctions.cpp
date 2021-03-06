
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/


#include <fstream>
#include <png.h>

#include "WriteFunctions.hpp"

#include "Simulation/ParticleCloud.hpp"
#include "Simulation/Settings.hpp"

//********************************************************************************
void writeRAWfile(std::string fileName, const ParticleCloud& cloud)
//********************************************************************************
{
    static int writeID = 0;
    if      (writeID<10   ) fileName+="0000";
    else if (writeID<100  ) fileName+="000";
    else if (writeID<1000 ) fileName+="00";
    else if (writeID<10000) fileName+="0";
    fileName+=std::to_string(writeID);
    std::string fileNameRAW = fileName+".raw";
    std::string fileNamePNG = fileName+".png";

    std::cout<<"#writing file "<<fileName<<std::endl;
    std::ofstream outfile (fileNameRAW.c_str(),std::ofstream::binary);
    outfile<<"#ID x y u v dens Fpress-x,y Fvisc-x,y Fother-x,y\n";
    
    const unsigned NParticles = unsigned(cloud.size());
    for (unsigned i=0;i<NParticles;i++)
    {
        Particle iParticle = cloud.particle(i);
        outfile<<std::scientific;
        outfile<<i<<"\t"                    
               <<iParticle.get<Attr::ePosition>().x<<"\t" 
               <<iParticle.get<Attr::ePosition>().y<<"\t" 
               <<iParticle.get<Attr::eVelocity>().x<<"\t" 
               <<iParticle.get<Attr::eVelocity>().y<<"\t" 
               <<iParticle.get<Attr::eDensity>()<<"\t" 
               <<iParticle.get<Attr::ePressForce>().x<<"\t" 
               <<iParticle.get<Attr::ePressForce>().y<<"\t" 
               <<iParticle.get<Attr::eViscForce>().x<<"\t" 
               <<iParticle.get<Attr::eViscForce>().y<<"\t" 
               <<iParticle.get<Attr::eOtherForce>().x<<"\t" 
               <<iParticle.get<Attr::eOtherForce>().y<<"\t" 
               <<"\n";                   
    }
    outfile.close();
    writeID++;
}

//********************************************************************************
glm::ivec3 floatToIntColor(glm::vec3 fCol, int res)
//********************************************************************************
{
    return glm::ivec3(fCol*float(res-1));
}

//********************************************************************************
void writePNG(std::string fileName, const std::vector<glm::vec3>& color)
//********************************************************************************
{
    fileName += ".png";
    std::cout<<"#writing file "<<fileName<<std::endl;

    FILE* outfile = fopen(fileName.c_str(), "wb");
    if (!outfile) 
    {
        std::cerr<<"ERROR::writePNG::couldn't open file '"<<fileName<<"'"<<std::endl;
        exit(1);
    }

    png_structp png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) 
    {
        std::cerr<<"ERROR::writePNG::couldn't create write struct"<<std::endl;
        exit(1);
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
    {
        png_destroy_write_struct (&png_ptr, (png_infopp)NULL);
        std::cerr<<"ERROR::writePNG::couldn't create info struct"<<std::endl;
        exit(1);
    }

    if (setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(outfile);
        std::cerr<<"ERROR::writePNG::stjmp smth"<<std::endl;
        exit(1);
    }

    png_init_io(png_ptr, outfile);
    // Write header (8 bit colour depth)
    png_set_IHDR(png_ptr, info_ptr, unsigned(RenderSettings::width), unsigned(RenderSettings::height),
        8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);

    // Allocate memory for one row (3 bytes per pixel - RGB)
    png_bytep row = NULL;
    row = (png_bytep) malloc(3 * RenderSettings::width * sizeof(png_byte));

    // Write image data
    for (unsigned j=0 ; j<RenderSettings::height ; j++) 
    {
        for (unsigned i=0 ; i<RenderSettings::width ; i++) 
        {
            glm::ivec3 iColor = floatToIntColor(color[j*RenderSettings::width+i],256);
            row[i*3+0] = (unsigned char)(std::min(iColor.r,255));
            row[i*3+1] = (unsigned char)(std::min(iColor.g,255));
            row[i*3+2] = (unsigned char)(std::min(iColor.b,255));
        }
        png_write_row(png_ptr, row);
    }

    // End write
    png_write_end(png_ptr, NULL);

    if (outfile != NULL) fclose(outfile);
    if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
    if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    if (row != NULL) free(row);
}

//********************************************************************************
void renderImage(std::string fileName, const ParticleCloud& cloud, const HashTable& neibhs)
//********************************************************************************
{
    static int writeID = 0;

    if      (writeID<10   ) fileName+="0000";
    else if (writeID<100  ) fileName+="000";
    else if (writeID<1000 ) fileName+="00";
    else if (writeID<10000) fileName+="0";
    fileName+=std::to_string(writeID);

    const double denomHEIGTH = 1./double(RenderSettings::height-1);
    const double denomWIDTH  = 1./double(RenderSettings::width-1);
    std::vector<glm::vec3> pixelColor(RenderSettings::width*RenderSettings::height);

    glm::dvec2 pos(0.0);
    glm::uvec2 gridPos(0);
    double potential = 0.;
    double potentialThres = 1./(0.4*SPHSettings::initDx*SPHSettings::initDx);
    glm::dvec2 vel(0.0);
    bool inFluid = false;
    double mult = RenderSettings::height/RenderSettings::width;
    //metaball apprach
    if (RenderSettings::fileRender==RenderSettings::METABALL)
    {
        const auto& particlePos = cloud.get<Attr::ePosition>();
        const auto& particleVel = cloud.get<Attr::eVelocity>();
        for (unsigned j=0;j<RenderSettings::height;j++)
        {
            pos.y = mult*BoundaryConditions::bndBox.maxX()-double(j)*denomHEIGTH*(mult*BoundaryConditions::bndBox.dx());
            for (unsigned i=0;i<RenderSettings::width;i++)
            {
                unsigned iNode = j*RenderSettings::width+i;
                pixelColor[iNode] = glm::vec3(0.0f);
                potential = 0.;
                inFluid = false;

                pos.x   = BoundaryConditions::bndBox.minX()+double(i)*denomWIDTH*(BoundaryConditions::bndBox.dx());
                gridPos =  neibhs.findGridPos(pos);
                std::vector<unsigned> neiParts = neibhs.neiParticlesFor(gridPos);
                for (unsigned nei : neiParts)
                {
                    double dist = glm::length(pos-particlePos[nei]);
                    if (dist>1.e-10) 
                    {
                        double potCon = 1./(dist*dist);
                        potential+= potCon;
                        vel += particleVel[nei]*potCon;
                    }
                    inFluid = (potential>=potentialThres);
                    if (inFluid) break;
                }
                if (inFluid) 
                {
                    vel/=potential;
                    float velMag = float(glm::length(vel));
                    float scale = 0.2f;
                    pixelColor[iNode] = glm::vec3(
                                    scale*velMag,
                                    scale*velMag, 
                                    scale*velMag*0.3f+0.7
                                    );
                }
            }
        }
    }
    //point shader approach
    else if (RenderSettings::fileRender==RenderSettings::DISCRETE)
    {
        const auto& particlePos = cloud.get<Attr::ePosition>();
        const auto& particleVel = cloud.get<Attr::eVelocity>();
        for (unsigned j=0;j<RenderSettings::height;j++)
        {
            pos.y = mult*BoundaryConditions::bndBox.maxX()-double(j)*denomHEIGTH*(mult*BoundaryConditions::bndBox.dx());
            for (unsigned i=0;i<RenderSettings::width;i++)
            {
                unsigned iNode = j*RenderSettings::width+i;
                pixelColor[iNode] = glm::vec3(0.0f);
                potential = 0.;
                inFluid = false;

                pos.x   = BoundaryConditions::bndBox.minX()+double(i)*denomWIDTH*(BoundaryConditions::bndBox.dx());
                gridPos =  neibhs.findGridPos(pos);
                std::vector<unsigned> neiParts = neibhs.neiParticlesFor(gridPos);
                for (unsigned nei : neiParts)
                {
                    double dist = glm::length(pos-particlePos[nei]);
                    if (dist<SPHSettings::initDx*0.4) inFluid=true;
                    if (inFluid) 
                    {
                        float velMag = float(glm::length(particleVel[nei]));
                        float scale = 0.2f;
                        pixelColor[iNode] = glm::vec3(
                                        scale*velMag,
                                        scale*velMag, 
                                        scale*velMag*0.3f+0.7
                                        );
                        break;
                    }
                }
                if (inFluid) break;
            }
        }
    }
    writePNG(fileName, pixelColor);
    writeID++;
}
