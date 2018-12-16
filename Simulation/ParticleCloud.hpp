
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    ParticleCloud
 
\************************************************************************/

#ifndef PARTICLECLOUD_H
#define PARTICLECLOUD_H

#include <vector>

#include "Particle.hpp"

//this is an itermediate class used to migrate from old system to new.
//once the migration is done this should be deleted.
class LesserParticle
{
public:
    glm::dvec2 position;
    glm::dvec2 velocity;
    glm::dvec2 normal;
    glm::dvec2 Fpress;
    glm::dvec2 Fvisc;
    glm::dvec2 Fsurf;
    glm::dvec2 Fother;
    glm::dvec2 Ftot;
    glm::ivec2 gridPos; //in background grid
    double mass;
    double density;
    double ddensity;    //1./density
    double densityErr;
    double pressure;
    std::vector<Neigbhor> nei; //list of neighbours

public:
    LesserParticle():
    position(-100.0),
    velocity(0.0),
    normal(0.0),
    Fpress(0.0),
    Fvisc(0.0),
    Fsurf(0.0),
    Fother(0.0),
    Ftot(0.0),
    gridPos(0),
    mass(0.0),
    density(0.0),
    ddensity(0.0),
    densityErr(0.0),
    pressure(0.0),
    nei()
    {}

    LesserParticle(const Particle& p):
    position   (p.position),
    velocity   (p.velocity),
    normal     (p.normal),
    Fpress     (p.Fpress),
    Fvisc      (p.Fvisc),
    Fsurf      (p.Fsurf),
    Fother     (p.Fother),
    Ftot       (p.Ftot),
    gridPos    (p.gridPos),
    mass       (p.mass),
    density    (p.density),
    ddensity   (p.ddensity),
    densityErr (p.densityErr),
    pressure   (p.pressure),
    nei        (p.nei)
    {}
};


class ParticleCloud
{
public:
    using ParticleVector = std::vector<LesserParticle>;
    using iterator = ParticleVector::iterator;
    using const_iterator = ParticleVector::const_iterator;

          iterator  begin()       { return _cloud.begin(); }
    const_iterator  begin() const { return _cloud.begin(); }
    const_iterator cbegin() const { return _cloud.cbegin(); }
          iterator  end()       { return _cloud.end(); }
    const_iterator  end() const { return _cloud.end(); }
    const_iterator cend() const { return _cloud.cend(); }

          LesserParticle& front()       { return _cloud.front(); }
    const LesserParticle& front() const { return _cloud.front(); }
          LesserParticle& back()       { return _cloud.back(); }
    const LesserParticle& back() const { return _cloud.back(); }

    ParticleCloud(): _cloud() {}
    ParticleCloud(const ParticleCloud&) = default;
    ParticleCloud(ParticleCloud&&)      = default;

    ParticleCloud(size_t s): _cloud(s) {}

    ~ParticleCloud() = default;

    ParticleCloud& operator=(const ParticleCloud&) = default;
    ParticleCloud& operator=(ParticleCloud&&)      = default;

    LesserParticle& operator[](size_t i) { return _cloud[i]; }
    const LesserParticle& operator[](size_t i) const { return _cloud[i]; }

    void reserve(size_t s) { _cloud.reserve(s); }

    void push_back(const Particle& p) {
      _cloud.push_back(LesserParticle{p});
    }

    void push_back(const LesserParticle& p) {
      _cloud.push_back(p);
    }

    size_t size() const { return _cloud.size(); }

    bool empty() const { return _cloud.empty(); }

private:
    ParticleVector _cloud;
};

#endif
