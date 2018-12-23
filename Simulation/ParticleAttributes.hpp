
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Description
    Particle Attributes and related boilerplate code

SourceFiles
    -

\************************************************************************/

#ifndef PARTICLEATTRIBUTES_H
#define PARTICLEATTRIBUTES_H

#include <glm/glm.hpp>
#include <vector>
#include <tuple>
#include <type_traits>
#include <utility>

#include "Neighbour.hpp"

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

using ParticleDataTypes = std::tuple<
  glm::dvec2,           //position
  glm::dvec2,           //velocity
  glm::dvec2,           //normal
  glm::dvec2,           //Fpress
  glm::dvec2,           //Fvisc
  glm::dvec2,           //Fsurf
  glm::dvec2,           //Fother
  glm::dvec2,           //Ftot
  double,               //mass
  double,               //density
  double,               //ddensity
  double,               //densityErr
  double,               //pressure
  std::vector<Neigbhor> // nei
>;
  
namespace AttrUtil
{
  template <typename TupleType, size_t I>
  using GenericTupleI = std::remove_reference_t <
    decltype (
      std::get<I> (
        std::declval<TupleType>()
      )
    )
  >;

  template <size_t I>
  using AttributeTypeI = GenericTupleI<ParticleDataTypes, I>;

  template <size_t I>
  using AttributeVectorTypeI = std::vector<AttributeTypeI<I>>;

  template<std::size_t... Is>
  constexpr decltype(auto) expandToType (std::index_sequence<Is...>) {
    return std::tuple< AttributeTypeI<Is>... >{ AttributeTypeI<Is>{}... };
  }

  template<std::size_t... Is>
  constexpr decltype(auto) expandToVectorType (std::index_sequence<Is...>) {
    return std::tuple< AttributeVectorTypeI<Is>... >{ AttributeVectorTypeI<Is>{}... };
  }

  template <typename VectorTupleType, size_t I>
  using TypeAtVectorTypeI = typename GenericTupleI<VectorTupleType, I>::value_type;

  template<typename VectorTupleType, size_t... Is>
  constexpr decltype(auto) fromVectorTypeToSimpleType_impl (std::index_sequence<Is...>) {
    return std::tuple< TypeAtVectorTypeI<VectorTupleType, Is>... >{ TypeAtVectorTypeI<VectorTupleType, Is>{}... };
  }

  template<typename VectorTupleType>
  constexpr decltype(auto) fromVectorTypeToSimpleType () {
    return fromVectorTypeToSimpleType_impl<VectorTupleType>(std::make_index_sequence<std::tuple_size<VectorTupleType>::value>());
  }

}

#endif
