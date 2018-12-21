
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

namespace detail {
  template<class Tuple, class Functor, std::size_t... Is>
  void applyOnTuple(Tuple& t, const Functor& f, std::index_sequence<Is...>)
  {
    (f(std::get<Is>(t)),...);
  }
}

//order must be in aggrement with data tuple
enum class Attr : size_t {
  ePosition = 0,
  eVelocity,
  eNormal,
  ePressForce,
  eViscForce,
  eSurfForce,
  eOtherForce,
  eTotalForce,
  eMass,
  eDensity,
  eDDensity,
  eDensErr,
  ePressure,
  eNei,
  nAttr
};

class ParticleCloud
{
private:
    using DataType = std::tuple<
      std::vector<glm::dvec2> //position
    , std::vector<glm::dvec2> //velocity
    , std::vector<glm::dvec2> //normal
    , std::vector<glm::dvec2> //Fpress
    , std::vector<glm::dvec2> //Fvisc
    , std::vector<glm::dvec2> //Fsurf
    , std::vector<glm::dvec2> //Fother
    , std::vector<glm::dvec2> //Ftot
    , std::vector<double    > //mass
    , std::vector<double    > //density
    , std::vector<double    > //ddensity
    , std::vector<double    > //densityErr
    , std::vector<double    > //pressure
    , std::vector<std::vector<Neigbhor>> // list of neighbours
    >;
    static constexpr size_t nAttributes = std::tuple_size<DataType>::value;
public:

    ParticleCloud(): _data() {}
    ParticleCloud(const ParticleCloud&) = default;
    ParticleCloud(ParticleCloud&&)      = default;

    ~ParticleCloud() = default;

    ParticleCloud& operator=(const ParticleCloud&) = default;
    ParticleCloud& operator=(ParticleCloud&&)      = default;

    template <Attr attr>
    decltype(auto) get() { return std::get<static_cast<size_t>(attr)>(_data); }
    template <Attr attr>
    decltype(auto) get() const { return std::get<static_cast<size_t>(attr)>(_data); }

    void reserve(size_t s) {
      auto reserveFunctor = [s](auto& v){ v.reserve(s); };
      applyFunctor(reserveFunctor);
    }

    void push_back(const Particle& p) {
      get<Attr::ePosition>().push_back(p.position);
      get<Attr::eVelocity>().push_back(p.velocity);
      get<Attr::eNormal>().push_back(p.normal);
      get<Attr::ePressForce>().push_back(p.Fpress);
      get<Attr::eViscForce>().push_back(p.Fvisc);
      get<Attr::eSurfForce>().push_back(p.Fvisc);
      get<Attr::eOtherForce>().push_back(p.Fother);
      get<Attr::eTotalForce>().push_back(p.Ftot);
      get<Attr::eMass>().push_back(p.mass);
      get<Attr::eDensity>().push_back(p.density);
      get<Attr::eDDensity>().push_back(p.ddensity);
      get<Attr::eDensErr>().push_back(p.densityErr);
      get<Attr::ePressure>().push_back(p.pressure);
      get<Attr::eNei>().push_back(p.nei);
    }
    
    Particle particle(size_t i) const
    {
      return Particle {
        get<Attr::ePosition   >()[i],
        get<Attr::eVelocity   >()[i],
        get<Attr::eNormal     >()[i],
        get<Attr::ePressForce >()[i],
        get<Attr::eViscForce  >()[i],
        get<Attr::eSurfForce  >()[i],
        get<Attr::eOtherForce >()[i],
        get<Attr::eTotalForce >()[i],
        get<Attr::eMass       >()[i],
        get<Attr::eDensity    >()[i],
        get<Attr::eDDensity   >()[i],
        get<Attr::eDensErr    >()[i],
        get<Attr::ePressure   >()[i],
        get<Attr::eNei        >()[i],
      };
    }

    size_t size() const { return get<Attr::ePosition>().size(); }

    bool empty() const { return get<Attr::ePosition>().empty(); }

    template <class Functor>
    void applyFunctor(const Functor& f) {
      detail::applyOnTuple(_data, f, std::make_index_sequence<nAttributes>());
    }
    template <class Functor>
    void applyFunctor(const Functor& f) const {
      detail::applyOnTuple(_data, f, std::make_index_sequence<nAttributes>());
    }

private:
    DataType _data; 
};

#endif
