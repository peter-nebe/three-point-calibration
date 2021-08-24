/*
 * three-point-calibration
 * Copyright (c) 2020 Peter Nebe (mail@peter-nebe.dev)
 *
 * This file is part of three-point-calibration.
 *
 * three-point-calibration is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * three-point-calibration is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with three-point-calibration.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_

#include "boost/qvm/vec.hpp"
#include "boost/qvm/mat.hpp"
#include <array>

namespace threePointCalibration
{

typedef double Coordinate_t;

template<std::size_t ndims>
using Vector_ = boost::qvm::vec<Coordinate_t, ndims>;

template<std::size_t ndims>
using Matrix_ = boost::qvm::mat<Coordinate_t, ndims, ndims>;

template<std::size_t ndims>
struct Point_;

template<>
struct Point_<2> : Vector_<2>
{
  const Coordinate_t &x() const { return a[0]; }
  const Coordinate_t &y() const { return a[1]; }
};

template<>
struct Point_<3> : Vector_<3>
{
  const Coordinate_t &x() const { return a[0]; }
  const Coordinate_t &y() const { return a[1]; }
  const Coordinate_t &z() const { return a[2]; }
};

template<std::size_t ndims, std::size_t npoints>
using Points_ = std::array<Point_<ndims>, npoints>;

template<std::size_t ndims>
using ReferencePoints_ = Points_<ndims, ndims>;

template<std::size_t ndims>
class Transformation_
{
public:
  using Vec = Vector_<ndims>;
  using Mat = Matrix_<ndims>;
  using RefPoints = ReferencePoints_<ndims>;

  Transformation_(const RefPoints &rp, const RefPoints &rpMapping);
  Transformation_(const RefPoints &triangleInPlane);

  Vec transform(const Vec &x) const
  {
    return _a * x + _b;
  }

  Vec transformInv(const Vec &x) const
  {
    return _aInv * (x - _b);
  }

  Point_<ndims> operator()(const Point_<ndims> &p) const
  {
    return transform(p);
  }

private:
  Mat _a, _aInv;
  Vec _b;
};

using Point = Point_<3>;
using Point2 = Point_<2>;
using Transformation = Transformation_<3>;
using Transformation2D = Transformation_<2>;

} /* namespace threePointCalibration */

namespace boost { namespace qvm
{
  template<std::size_t ndims>
  struct vec_traits<threePointCalibration::Point_<ndims>> : vec_traits<threePointCalibration::Vector_<ndims>>
  {
  };
}}

#endif /* TRANSFORMATION_H_ */
