
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    BoundingBox
 
Description
    AA bounding box

SourceFiles
    -

\************************************************************************/

#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include <glm/glm.hpp>

//created for M-> {glm::dvec2 or glm::vec2}
template <class M>
class BoundingBox
{
protected:
    M boxMinPos_;
    M boxMaxPos_;
    M delta_;

    using CoordinateType = typename M::value_type;
public:
    BoundingBox():
    boxMinPos_(0.0),
    boxMaxPos_(0.0),
    delta_(0.0)
    {}

    BoundingBox(const M& minPos, const M& maxPos):
    boxMinPos_(minPos),
    boxMaxPos_(maxPos),
    delta_(maxPos-minPos)
    {}

    void setBoundingBox(const M& minPos, const M& maxPos)
    {
        boxMinPos_ = minPos;
        boxMaxPos_ = maxPos;
        delta_     = maxPos-minPos;
    }
    
    const M& minPos() const { return boxMinPos_;}
    const M& maxPos() const { return boxMaxPos_;}
    const M& delta () const { return delta_;}

    const CoordinateType& minX() const { return boxMinPos_.x; }
    const CoordinateType& minY() const { return boxMinPos_.y; }
    const CoordinateType& maxX() const { return boxMaxPos_.x; }
    const CoordinateType& maxY() const { return boxMaxPos_.y; }
    const CoordinateType& dx()   const { return delta_.x; }
    const CoordinateType& dy()   const { return delta_.y; }
};

#endif

