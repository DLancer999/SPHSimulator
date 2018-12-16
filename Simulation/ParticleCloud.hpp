
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
#include <tuple>

#include "Particle.hpp"

//this is an itermediate class used to migrate from old system to new.
//once the migration is done this should be deleted.
class LesserParticle
{
public:
    //glm::dvec2 position;
    //glm::dvec2 velocity;
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
    //position(-100.0),
    //velocity(0.0),
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
    //position   (p.position),
    //velocity   (p.velocity),
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

//order must be in aggrement with data tuple
enum class Attr : size_t {
  ePosition = 0,
  eVelocity,
  nAttr
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

    ParticleCloud(): _cloud(), _data() {}
    ParticleCloud(const ParticleCloud&) = default;
    ParticleCloud(ParticleCloud&&)      = default;

    ~ParticleCloud() = default;

    ParticleCloud& operator=(const ParticleCloud&) = default;
    ParticleCloud& operator=(ParticleCloud&&)      = default;

    LesserParticle& operator[](size_t i) { return _cloud[i]; }
    const LesserParticle& operator[](size_t i) const { return _cloud[i]; }

    template <Attr attr>
    decltype(auto) get() { return std::get<static_cast<size_t>(attr)>(_data); }
    template <Attr attr>
    decltype(auto) get() const { return std::get<static_cast<size_t>(attr)>(_data); }

    void reserve(size_t s) {
      _cloud.reserve(s);
      get<Attr::ePosition>().reserve(s);
      get<Attr::eVelocity>().reserve(s);
    }

    void push_back(const Particle& p) {
      _cloud.push_back(LesserParticle{p});
      get<Attr::ePosition>().push_back(p.position);
      get<Attr::eVelocity>().push_back(p.velocity);
    }

    //void push_back(const LesserParticle& p) {
    //  _cloud.push_back(p);
    //}
    
    Particle particle(size_t i) const
    {
      return Particle {
        get<Attr::ePosition>()[i],
        get<Attr::eVelocity>()[i],
        _cloud[i].normal,
        _cloud[i].Fpress,
        _cloud[i].Fvisc,
        _cloud[i].Fsurf,
        _cloud[i].Fother,
        _cloud[i].Ftot,
        _cloud[i].gridPos,
        _cloud[i].mass,
        _cloud[i].density,
        _cloud[i].ddensity,
        _cloud[i].densityErr,
        _cloud[i].pressure,
        _cloud[i].nei
      };
    }

    size_t size() const { return _cloud.size(); }

    bool empty() const { return _cloud.empty(); }

    ParticleVector& getCloud() { return _cloud; }
    const ParticleVector& getCloud() const { return _cloud; }

private:
    ParticleVector _cloud;
    std::tuple<
      std::vector<glm::dvec2> //position
    , std::vector<glm::dvec2> //velocity
    //, std::vector<glm::dvec2> //normal
    //, std::vector<glm::dvec2> //Fpress
    //, std::vector<glm::dvec2> //Fvisc
    //, std::vector<glm::dvec2> //Fsurf
    //, std::vector<glm::dvec2> //Fother
    //, std::vector<glm::dvec2> //Ftot
    //, std::vector<glm::ivec2> //gridPos in background grid
    //, std::vector<double>     //mass
    //, std::vector<double>     //density
    //, std::vector<double>     //ddensity
    //, std::vector<double>     //densityErr
    //, std::vector<double>     //pressure
    //, std::vector<std::vector<Neigbhor>> // list of neighbours
    > _data;
};

#endif
