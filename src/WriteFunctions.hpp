
#ifndef WRITEFUNCTIONS_H
#define WRITEFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#include "Particle.hpp"
#include "HashTable.hpp"

void writeGNUfile(std::string fileName, std::vector<Particle>& cloud);
void renderImage(std::string fileName, std::vector<Particle>& cloud, HashTable& neibhs);
void writePNG(std::string fileName, std::vector<glm::vec3>& color);
glm::ivec3 floatToIntColor(glm::vec3 fCol, int res);

#endif
