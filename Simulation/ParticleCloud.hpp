
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

class ParticleCloud
{
public:
    using ParticleVector = std::vector<Particle>;
    using iterator = ParticleVector::iterator;
    using const_iterator = ParticleVector::const_iterator;

          iterator  begin()       { return _cloud.begin(); }
    const_iterator  begin() const { return _cloud.begin(); }
    const_iterator cbegin() const { return _cloud.cbegin(); }
          iterator  end()       { return _cloud.end(); }
    const_iterator  end() const { return _cloud.end(); }
    const_iterator cend() const { return _cloud.cend(); }

          Particle& front()       { return _cloud.front(); }
    const Particle& front() const { return _cloud.front(); }
          Particle& back()       { return _cloud.back(); }
    const Particle& back() const { return _cloud.back(); }

    ParticleCloud(): _cloud() {}
    ParticleCloud(const ParticleCloud&) = default;
    ParticleCloud(ParticleCloud&&)      = default;

    ParticleCloud(size_t s): _cloud(s) {}

    ~ParticleCloud() = default;

    ParticleCloud& operator=(const ParticleCloud&) = default;
    ParticleCloud& operator=(ParticleCloud&&)      = default;

    Particle& operator[](size_t i) { return _cloud[i]; }
    const Particle& operator[](size_t i) const { return _cloud[i]; }

    void reserve(size_t s) { _cloud.reserve(s); }

    void push_back(const Particle& p) { _cloud.push_back(p); }

    size_t size() const { return _cloud.size(); }

private:
    ParticleVector _cloud;
};

#endif
