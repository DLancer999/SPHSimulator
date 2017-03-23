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
    
    const M& minPos() const {return boxMinPos_;}
    const M& maxPos() const {return boxMaxPos_;}
    const M& delta () const {return delta_;}

    const typename M::value_type& minX() const { return boxMinPos_.x; }
    const typename M::value_type& minY() const { return boxMinPos_.y; }
    const typename M::value_type& maxX() const { return boxMaxPos_.x; }
    const typename M::value_type& maxY() const { return boxMaxPos_.y; }
    const typename M::value_type& dx()   const { return delta_.x; }
    const typename M::value_type& dy()   const { return delta_.y; }
};

#endif

