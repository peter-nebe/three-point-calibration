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
#include <tuple>
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

template<typename VecTp>
bool isMagOne(const VecTp &v)
{
  return almost_equal(qvm::mag(v), 1.0);
}

template<typename VecTp>
bool isDotZero(const VecTp &a, const VecTp &b)
{
  return std::fabs(qvm::dot(a, b)) < 1e-15;
}

template<typename MatTp>
bool isDetOne(const MatTp &m)
{
  return almost_equal(qvm::determinant(m), 1.0);
}

Matrix_<2> makeRotation(const ReferencePoints_<2> &rp, const ReferencePoints_<2> &rpMapping)
{
  auto normalizedVector = [](const ReferencePoints_<2> &rp)
                          {
                            const Vector_<2> v = qvm::normalized(rp[1] - rp[0]);
                            return std::make_tuple(v.a[0], v.a[1]);
                          };

  const auto [dx, dy] = normalizedVector(rp);
  const auto [dxm, dym] = normalizedVector(rpMapping);
  const auto xBaseX = dx * dxm + dy * dym;
  const auto xBaseY = dy * dxm - dx * dym;

  const Vector_<2> xBase{ xBaseX, xBaseY };
  const Vector_<2> yBase{ -xBaseY, xBaseX };
  assert(isMagOne(xBase));
  assert(isMagOne(yBase));
  assert(isDotZero(xBase, yBase));

  Matrix_<2> rot;
  qvm::col<0>(rot) = xBase;
  qvm::col<1>(rot) = yBase;
  assert(isDetOne(rot));

  return rot;
}

} // namespace

template<>
Transformation_<2>::Transformation_(const RefPoints &rp, const RefPoints &rpMapping)
{
  // set rotation matrix
  _a = makeRotation(rp, rpMapping);

  // rotation matrix: inverse = transpose
  _aInv = qvm::transposed(_a);
  assert(isDetOne(_aInv));

  // set translation vector
  _b = rp.front() - _a * rpMapping.front();
}

template<>
Transformation_<3>::Transformation_(const RefPoints &triangleInPlane)
{
  // calculate orientation of the plane in mapping coordinates

  const Vec p = triangleInPlane[0];
  const Vec u = triangleInPlane[1] - p;
  const Vec v = triangleInPlane[2] - p;

  const Vec n0 = qvm::normalized(qvm::cross(u, v));
  assert(isMagOne(n0));
  const Vec zBase = -n0;

  // dot(yBase, zBase) = 0
  // ybx zbx + yby zby + ybz zbz = 0
  const Coordinate_t ybx = 0;
  const Coordinate_t ybz = 1;
  // yby zby + zbz = 0
  // yby zby = -zbz
  // yby = -zbz / zby
  const Coordinate_t yby = -zBase.a[2] / zBase.a[1];

  const Vec yBase = qvm::normalized(Vec{ ybx, yby, ybz });
  const Vec xBase = qvm::cross(yBase, zBase);
  assert(isMagOne(xBase));
  assert(isMagOne(yBase));
  assert(isMagOne(zBase));
  assert(isDotZero(xBase, yBase));
  assert(isDotZero(yBase, zBase));
  assert(isDotZero(zBase, xBase));

  Mat rot;
  qvm::col<0>(rot) = xBase;
  qvm::col<1>(rot) = yBase;
  qvm::col<2>(rot) = zBase;

  // set inverse rotation matrix
  _aInv = rot;

  // rotation matrix: inverse = transpose
  _a = qvm::transposed(_aInv);
  assert(isDetOne(_a));
  assert(isDetOne(_aInv));

  const Coordinate_t distFromOrigin = qvm::dot(p, n0);
  assert(distFromOrigin > 0);

  // set translation vector
  _b = Vec{ 0, 0, distFromOrigin };
}

template<>
Transformation_<3>::Transformation_(const RefPoints &rp, const RefPoints &rpMapping)
{
  auto makeLinearTransformation = [](const RefPoints &threePoints)
                                  {
                                    const Vec p = threePoints[0];
                                    const Vec u = threePoints[1] - p;
                                    const Vec v = threePoints[2] - p;

                                    Mat m;
                                    qvm::col<0>(m) = u;
                                    qvm::col<1>(m) = v;
                                    qvm::col<2>(m) = qvm::cross(u, v);

                                    return m;
                                  };
  const Mat m = makeLinearTransformation(rp);
  const Mat mm = makeLinearTransformation(rpMapping);

  _a = m * qvm::inverse(mm);
  _aInv = qvm::inverse(_a);

  // set translation vector
  _b = rp.front() - _a * rpMapping.front();
}

} /* namespace threePointCalibration */
