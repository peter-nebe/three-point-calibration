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

#include "boost/qvm/mat.hpp"
#include "boost/qvm/vec_operations.hpp"
#include "boost/qvm/vec_traits_defaults.hpp"
#include "boost/qvm/deduce_vec.hpp"
#include <array>

namespace threePointCalibration
{

typedef double Coordinate_t;

template<int Dim>
struct Point_;
using Point2 = Point_<2>;
using Point3 = Point_<3>;

template<>
struct Point_<2>
{
  Coordinate_t x, y;
};

template<>
struct Point_<3>
{
  Coordinate_t x, y, z;
};

template<int Dim>
using Vector_ = Point_<Dim>;

template<int Dim>
using Matrix_ = boost::qvm::mat<Coordinate_t, Dim, Dim>;

template<int Dim>
using ReferencePoints_ = std::array<Point_<Dim>, Dim>;

template<int Dim>
class Transformation_
{
public:
  using RefPoints = ReferencePoints_<Dim>;
  using Point = Point_<Dim>;
  using Mat = Matrix_<Dim>;
  using Vec = Vector_<Dim>;

  Transformation_(const RefPoints &rp, const RefPoints &rpMapping);
  Transformation_(const RefPoints &triangleInPlane);

  Point transform(const Point &x) const
  {
    using boost::qvm::operator+;
    return _a * x + _b;
  }

  Point transformInv(const Point &x) const
  {
    using boost::qvm::operator-;
    return _aInv * (x - _b);
  }

  Point operator()(const Point &p) const
  {
    return transform(p);
  }

private:
  Mat _a, _aInv;
  Vec _b;
}; // class Transformation_

using Transformation = Transformation_<3>;
using Transformation2D = Transformation_<2>;

} // namespace threePointCalibration

namespace boost { namespace qvm
{

namespace tpc = threePointCalibration;

template<int Dim>
struct vec_traits<tpc::Point_<Dim>> : vec_traits_defaults<tpc::Point_<Dim>, tpc::Coordinate_t, Dim>
{
  template<int I>
  static tpc::Coordinate_t &write_element(tpc::Point_<Dim>&);
};
template<> template<>
inline tpc::Coordinate_t &vec_traits<tpc::Point2>::write_element<0>(tpc::Point2 &p)
{
  return p.x;
}
template<> template<>
inline tpc::Coordinate_t &vec_traits<tpc::Point2>::write_element<1>(tpc::Point2 &p)
{
  return p.y;
}
template<> template<>
inline tpc::Coordinate_t &vec_traits<tpc::Point3>::write_element<0>(tpc::Point3 &p)
{
  return p.x;
}
template<> template<>
inline tpc::Coordinate_t &vec_traits<tpc::Point3>::write_element<1>(tpc::Point3 &p)
{
  return p.y;
}
template<> template<>
inline tpc::Coordinate_t &vec_traits<tpc::Point3>::write_element<2>(tpc::Point3 &p)
{
  return p.z;
}

template<int Dim>
struct deduce_vec2<tpc::Matrix_<Dim>, tpc::Point_<Dim>, Dim>
{
  typedef tpc::Point_<Dim> type;
};

}} // namespace boost::qvm

#endif // TRANSFORMATION_H_
