
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Function
    ZValue
 
Description
    Utility function that computes the z-curve value of a (x,y) point
    This algorithm works for values up to half the capacity of Int.

SourceFiles
    -

\************************************************************************/

#ifndef ZValue_H
#define ZValue_H

#include <stdlib.h>
#include <type_traits>

namespace Util
{
  template <typename Int>
  Int ZValue(Int x, Int y) 
  {
      static_assert(std::is_integral<Int>::value);

      constexpr size_t nBits = sizeof(Int)*8;
      constexpr size_t halfBits = nBits/2;

      Int result = 0;
      size_t activeBit = 1;
      size_t bitOffset = 0;
      for(size_t ii = 0; ii < halfBits; ++ii)
      {
          result |= (y & activeBit) << bitOffset;
          result |= (x & activeBit) << (bitOffset+1);
          ++bitOffset;
          activeBit<<=1;
      }
      return result;
  }
}

#endif
