
/*************************************************************************\
License
    Copyright (c) 2017 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

\************************************************************************/

#include "Reorderer.hpp"

#include "Kernels.hpp"
#include "Statistics/Statistics.hpp"
#include "Util/ZValue.hpp"

#include <boost/range/irange.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/indexed.hpp>

using boost::adaptors::transformed;
using boost::adaptors::indexed;

namespace
{
  template <typename RandomAccessCont, typename RandomAccessIndexes>
  void indexedReordering (RandomAccessCont& c, const RandomAccessIndexes& indexes)
  {
    auto indexedValuesR = indexes | transformed([&c]( auto i ){ return c[i]; });

    RandomAccessCont tempCont(std::size(indexes));
    boost::copy(indexedValuesR, std::begin(tempCont));
    boost::copy(tempCont, std::begin(c));
  }
}

//********************************************************************************
void Reorderer::reorderCloud(std::vector<Particle>& cloud, const unsigned NParticles)
//********************************************************************************
{
    static auto reorTimerID = Statistics::createTimer("Reorderer::reorderTimer");
    Statistics::TimerGuard reorderTimerGuard(reorTimerID);

    glm::dvec2 minPos = cloud.front().position;
    for (unsigned i=1; i<NParticles; ++i) {
      const auto& iPos = cloud[i].position;
      if (minPos.x > iPos.x)
        minPos.x = iPos.x;
      if (minPos.y > iPos.y)
        minPos.y = iPos.y;
    }

    auto getPos = [&cloud](unsigned i) { return cloud[i].position; };
    auto calcZIndex = [&minPos](const glm::dvec2& pos) {
      const glm::uvec2 gridIndex = floor((pos-minPos)*Kernel::SmoothingLength::dh);
      return Util::ZValue(gridIndex.x, gridIndex.y);
    };

    auto indexesR = boost::irange(0u, NParticles);
    auto zValR = indexesR | transformed(getPos) | transformed(calcZIndex);

    std::vector<unsigned> zVal(zValR.begin(), zValR.end());

    std::vector<unsigned> indexes(indexesR.begin(), indexesR.end());
    std::sort(indexes.begin(),indexes.end(), [&zVal](const unsigned& i, const unsigned& j) {
      return zVal[i]<zVal[j];
    });

    ::indexedReordering(cloud, indexes);
}



