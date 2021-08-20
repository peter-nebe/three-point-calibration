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

#include "transformation.h"
#include "boost/qvm/vec_operations.hpp"
#include "boost/qvm/vec_mat_operations.hpp"
#include "boost/qvm/mat_operations.hpp"
#include "boost/qvm/map_mat_vec.hpp"
#include "boost/qvm/map_mat_mat.hpp"
namespace qvm = boost::qvm;

namespace threePointCalibration
{

namespace
{

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp = 2)
{
  const auto abs = std::fabs(x - y);

  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return abs <= std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp
      // unless the result is subnormal
      || abs < std::numeric_limits<T>::min();
}

bool isMagOne(const Vec_t &v)
{
  return almost_equal(qvm::mag(v), 1.0);
}

bool isDotZero(const Vec_t &a, const Vec_t &b)
{
  return std::fabs(qvm::dot(a, b)) < 1e-15;
}

} // namespace

Transformation::Transformation(const RefPoints_t &triangleInPlane)
{
  const Vec_t p = triangleInPlane[0];
  const Vec_t u = triangleInPlane[1] - p;
  const Vec_t v = triangleInPlane[2] - p;

  const Vec_t n0 = qvm::normalized(qvm::cross(u, v));
  assert(isMagOne(n0));

  const Vec_t baseZ = -n0;
  assert(isMagOne(baseZ));

  // dot(baseY, baseZ) = 0
  // byx bzx + byy bzy + byz bzz = 0
  const Coordinate_t byx = 0;
  const Coordinate_t byz = 1;
  // byy bzy + bzz = 0
  // byy bzy = -bzz
  // byy = -bzz / bzy
  const Coordinate_t byy = -baseZ.a[2] / baseZ.a[1];
  const Vec_t baseY = qvm::normalized(Vec_t{ byx, byy, byz });
  assert(isMagOne(baseY));
  assert(isDotZero(baseY, baseZ));

  const Vec_t baseX = qvm::cross(baseY, baseZ);
  assert(isMagOne(baseX));
  assert(isDotZero(baseX, baseZ));
  assert(isDotZero(baseX, baseY));

  // set rotation matrix
  qvm::col<0>(_a) = baseX;
  qvm::col<1>(_a) = baseY;
  qvm::col<2>(_a) = baseZ;
  assert(almost_equal(qvm::determinant(_a), 1.0));

  // rotation matrix: inverse = transpose
  _aInv = qvm::transposed(_a);
  assert(almost_equal(qvm::determinant(_aInv), 1.0));

  const Coordinate_t distFromOrigin = qvm::dot(p, n0);
  assert(distFromOrigin > 0);

  // set translation vector
  _b = n0 * distFromOrigin;
}

Transformation::Transformation(const RefPoints_t &refPts, const RefPoints_t &refPtsT)
: Transformation(calcOrigin(refPts, refPtsT), refPts, refPtsT)
{
}

Transformation::Transformation(const Vec_t &origin, const RefPoints_t &refPts, const RefPoints_t &refPtsT)
{
  Mat_t ref;
  qvm::col<0>(ref) = refPts[0] - origin;
  qvm::col<1>(ref) = refPts[1] - origin;
  qvm::col<2>(ref) = refPts[2] - origin;

  Mat_t refT;
  qvm::col<0>(refT) = refPtsT[0];
  qvm::col<1>(refT) = refPtsT[1];
  qvm::col<2>(refT) = refPtsT[2];

  _a = ref * qvm::inverse(refT);
  _aInv = qvm::inverse(_a);
  _b = origin;
}

Transformation::Transformation(const Vec_t &origin, const Vec_t &baseX, const Vec_t &baseY)
{
  const Vec_t x = baseX - origin;
  const Vec_t y = baseY - origin;

  qvm::col<0>(_a) = x;
  qvm::col<1>(_a) = y;
  qvm::col<2>(_a) = qvm::cross(x, y);

  _aInv = qvm::inverse(_a);
  _b = origin;
}

Vec_t Transformation::calcOrigin(const RefPoints_t &refPts, const RefPoints_t &refPtsT)
{
  const Transformation refT(refPtsT[0], refPtsT[1], refPtsT[2]);
  const Vec_t originT = refT.transformInv({});

  const Transformation ref(refPts[0], refPts[1], refPts[2]);
  const Vec_t origin = ref.transform(originT);

  return origin;
}

Vec_t Transformation::transform(const Vec_t &x) const
{
  return _a * x + _b;
}

Vec_t Transformation::transformInv(const Vec_t &x) const
{
  return _aInv * (x - _b);
}

Point_t Transformation::operator()(const Point_t &p) const
{
  return transform(p);
}

} /* namespace threePointCalibration */
