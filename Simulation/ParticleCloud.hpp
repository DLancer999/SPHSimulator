
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

#include "ParticleAttributes.hpp"
#include "Particle.hpp"

namespace detail {
  template<class Tuple, class Functor, std::size_t... Is>
  void applyOnTuple(Tuple& t, const Functor& f, std::index_sequence<Is...>)
  {
    (f(std::get<Is>(t)),...);
  }

  template <typename ParticleVectorData>
  class ParticleCloudT
  {
  private:
      static constexpr size_t nAttributes = std::tuple_size<ParticleVectorData>::value;
  public:

      ParticleCloudT(): _data() {}
      ParticleCloudT(const ParticleCloudT&) = default;
      ParticleCloudT(ParticleCloudT&&)      = default;

      ~ParticleCloudT() = default;

      ParticleCloudT& operator=(const ParticleCloudT&) = default;
      ParticleCloudT& operator=(ParticleCloudT&&)      = default;

      template <Attr attr>
      decltype(auto) get() { return std::get<static_cast<size_t>(attr)>(_data); }
      template <Attr attr>
      decltype(auto) get() const { return std::get<static_cast<size_t>(attr)>(_data); }

      void reserve(size_t s) {
        auto reserveFunctor = [s](auto& v){ v.reserve(s); };
        applyFunctor(reserveFunctor);
      }

      void push_back(const Particle& p) {
        push_back_impl<0>(p);
      }

      decltype(auto) particle(size_t i) const {
        return particle_impl( i,
          std::make_index_sequence<static_cast<size_t>(nAttributes)>()
        );
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
      ParticleVectorData _data; 

      template <size_t ... Is>
      decltype(auto) particle_impl(size_t i, std::index_sequence<Is...>) const
      {
        using RetType = decltype(
          AttrUtil::fromVectorTypeToSimpleType<ParticleVectorData>()
        );
        return RetType { std::make_tuple( get<static_cast<Attr>(Is)>()[i]...)};
      }

      template <size_t Is>
      void push_back_impl(const Particle& p) {
        if constexpr ( Is < nAttributes ) {
          get<static_cast<Attr>(Is)>().push_back(p.get<static_cast<Attr>(Is)>());
          push_back_impl<Is+1>(p);
        }
      }

  };
}

//TODO:: current implementation only works for this definitions of ParticleCloud and Particle.
//However, in theory we should be able to have Particle/ParticleCloud instances that support a different
//subset of the attributes. Current implementation hints at this possibility, but we are not there yet.
using ParticleCloud = detail::ParticleCloudT<decltype(
    AttrUtil::expandToVectorType(std::make_index_sequence<size_t(Attr::nAttr)>())
)>;

#endif
