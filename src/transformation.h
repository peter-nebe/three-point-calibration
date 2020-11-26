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
typedef boost::qvm::vec<Coordinate_t, 3> Vec_t;
typedef boost::qvm::mat<Coordinate_t, 3, 3> Mat_t;

struct Point_t : Vec_t
{
  const Coordinate_t &x() const { return a[0]; }
  const Coordinate_t &y() const { return a[1]; }
  const Coordinate_t &z() const { return a[2]; }
};

template<std::size_t N>
using Points_t = std::array<Point_t, N>;
using RefPoints_t = Points_t<3>;

class Transformation
{
public:
  Transformation(const RefPoints_t &refPts, const RefPoints_t &refPtsT);
  Vec_t transform(const Vec_t &x) const;
  Point_t operator()(const Point_t &p) const;

private:
  Transformation(const Vec_t &origin, const Vec_t &baseX, const Vec_t &baseY);
  Transformation(const Vec_t &origin, const RefPoints_t &refPts, const RefPoints_t &refPtsT);
  static Vec_t calcOrigin(const RefPoints_t &refPts, const RefPoints_t &refPtsT);
  Vec_t transformInv(const Vec_t &x) const;

  Mat_t a;
  Vec_t b;
};

} /* namespace threePointCalibration */

namespace boost { namespace qvm
{
  template<>
  struct vec_traits<threePointCalibration::Point_t> : vec_traits<threePointCalibration::Vec_t>
  {
  };
}}

#endif /* TRANSFORMATION_H_ */
