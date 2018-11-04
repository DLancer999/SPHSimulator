
/*************************************************************************\
License
    Copyright (c) 2018 Kavvadias Ioannis.
    
    This file is part of SPHSimulator.
    
    Licensed under the MIT License. See LICENSE file in the project root for 
    full license information.  

Class
    Matrix
 
Description
    This is an elegant vector<vector<T>>

SourceFiles
    -

\************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

namespace Util
{

template <typename T>
class Matrix
{
public:
  using iterator       = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator       begin()        { return _data.begin(); }
  const_iterator begin()  const { return _data.begin(); }
  const_iterator cbegin() const { return _data.cbegin(); }
  iterator       end()        { return _data.end(); }
  const_iterator end()  const { return _data.end(); }
  const_iterator cend() const { return _data.cend(); }

  Matrix() :_data() ,_sizeX(0) ,_sizeY(0) {};
  Matrix(const Matrix&) = default;
  Matrix(Matrix&&)      = default;
  ~Matrix() = default;
  Matrix& operator=(const Matrix&) = default;
  Matrix& operator=(Matrix&&)      = default;

  void resize(size_t sizeX, size_t sizeY)
  {
    _sizeX = sizeX;
    _sizeY = sizeY;
    _data.resize(sizeX*sizeY, T{});
  }
  void clear()
  {
    _data.clear();
  }

  T& operator()(size_t i, size_t j) {
    return _data[i*_sizeY + j];
  }

  const T& operator()(size_t i, size_t j) const {
    return _data[i*_sizeY + j];
  }

private:
  std::vector<T> _data;
  size_t _sizeX;
  size_t _sizeY;
};

}

#endif
