#ifndef KERNELS_H
#define KERNELS_H

#include <glm/glm.hpp>
namespace Kernel
{
    class SmoothingLength
    {
    public:
        static double h  ;
        static double h2 ;
        static double dh ;
        static double dh4;
        static double dh8;

        static void setSmoothingLength(const double);
    };

    namespace poly6
    {
        double         W(glm::dvec2& xi, glm::dvec2& xj);
        glm::dvec2 gradW(glm::dvec2& xi, glm::dvec2& xj);
        double     laplW(glm::dvec2& xi, glm::dvec2& xj);
    }
    namespace spiky
    {
        glm::dvec2 gradW(glm::dvec2& xi, glm::dvec2& xj);
    }
    namespace visc
    {
        double laplW(glm::dvec2& xi, glm::dvec2& xj);
    }
}

#endif
