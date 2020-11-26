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
namespace qvm = boost::qvm;

namespace threePointCalibration
{

Transformation::Transformation(const RefPoints_t &refPts, const RefPoints_t &refPtsT)
: Transformation(calcOrigin(refPts, refPtsT), refPts, refPtsT)
{
}

Vec_t Transformation::transform(const Vec_t &x) const
{
  return a * x + b;
}

Point_t Transformation::operator()(const Point_t &p) const
{
  return transform(p);
}

Transformation::Transformation(const Vec_t &origin, const Vec_t &baseX, const Vec_t &baseY)
{
  const Vec_t x = baseX - origin;
  const Vec_t y = baseY - origin;

  qvm::col<0>(a) = x;
  qvm::col<1>(a) = y;
  qvm::col<2>(a) = qvm::cross(x, y);

  b = origin;
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

  a = ref * qvm::inverse(refT);
  b = origin;
}

Vec_t Transformation::calcOrigin(const RefPoints_t &refPts, const RefPoints_t &refPtsT)
{
  const Transformation refT(refPtsT[0], refPtsT[1], refPtsT[2]);
  const Vec_t originT = refT.transformInv({});

  const Transformation ref(refPts[0], refPts[1], refPts[2]);
  const Vec_t origin = ref.transform(originT);

  return origin;
}

Vec_t Transformation::transformInv(const Vec_t &x) const
{
  return qvm::inverse(a) * (x - b);
}

} /* namespace threePointCalibration */
