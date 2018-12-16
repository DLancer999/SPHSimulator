
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

#include <numeric>
#include <limits>

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
    c = std::move(tempCont);
  }
}

//********************************************************************************
void Reorderer::reorderCloud(ParticleCloud& cloud)
//********************************************************************************
{
    static auto reorTimerID = Statistics::createTimer("Reorderer::reorderTimer");
    Statistics::TimerGuard reorderTimerGuard(reorTimerID);

    //if (cloud.empty())
    //  return;

    const auto& particlePos = cloud.get<Attr::ePosition>();

    glm::dvec2 minPos = glm::dvec2(std::numeric_limits<double>::max());
    for (const auto& iPos : particlePos ) {
      if (minPos.x > iPos.x)
        minPos.x = iPos.x;
      if (minPos.y > iPos.y)
        minPos.y = iPos.y;
    }

    auto calcZIndex = [&minPos](const glm::dvec2& pos) {
      const glm::uvec2 gridIndex = floor((pos-minPos)*Kernel::SmoothingLength::dh);
      return Util::ZValue(gridIndex.x, gridIndex.y);
    };

    auto zValR = particlePos | transformed(calcZIndex);

    std::vector<unsigned> zVal(zValR.begin(), zValR.end());

    std::vector<unsigned> indexes(cloud.size());
    std::iota(indexes.begin(), indexes.end(), 0);

    std::sort(indexes.begin(),indexes.end(), [&zVal](const unsigned& i, const unsigned& j) {
      return zVal[i]<zVal[j];
    });

    auto indexedReordered = [&indexes](auto& v){ ::indexedReordering(v, indexes); };
    cloud.applyFunctor(indexedReordered);
}



