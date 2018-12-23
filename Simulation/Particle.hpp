
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Particle 
 
Description
    All parameters required to set-up simulation

SourceFiles
    -

\************************************************************************/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <glm/glm.hpp>
#include <type_traits>
#include <utility>

#include "ParticleAttributes.hpp"

namespace detail
{
  //template <size_t ... DataIndexes>
  template <typename DataType>
  class ParticleT
  {
  public:
    ParticleT():_data{}{}
    ParticleT(const ParticleT&) = default;
    ParticleT(ParticleT&&)      = default;
    ParticleT(const DataType& data)
    :_data{data}
    {}
    ParticleT(DataType&& data)
    :_data{std::move(data)}
    {}

    ~ParticleT() = default;
    ParticleT& operator=(const ParticleT&) = default;
    ParticleT& operator=(ParticleT&&     ) = default;

    template <Attr attr>
    decltype(auto) get() { return std::get<static_cast<size_t>(attr)>(_data); }
    template <Attr attr>
    decltype(auto) get() const { return std::get<static_cast<size_t>(attr)>(_data); }
  private:
    DataType _data;
  };
}

using Particle = detail::ParticleT<decltype(
    AttrUtil::expandToType(std::make_index_sequence<size_t(Attr::nAttr)>())
)>;

#endif
