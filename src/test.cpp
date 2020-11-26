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
#include "boost/geometry/geometries/point_xyz.hpp"
#include "boost/geometry/geometries/point_xy.hpp"
#include "boost/geometry/strategies/transform/matrix_transformers.hpp"
#include "boost/geometry/algorithms/transform.hpp"
#include "boost/qvm/vec_operations.hpp"
#include <functional>
#include <algorithm>
#include <iostream>
#include <iomanip>
using namespace threePointCalibration;
using namespace std;
namespace geo = boost::geometry;

typedef geo::model::d3::point_xyz<Coordinate_t> GeoPnt_t;
typedef geo::model::d2::point_xy<Coordinate_t> Geo2dPnt_t;

struct Map2d
{
protected:
  Coordinate_t unmapped = 0;
};

struct MapXY : Map2d
{
  Geo2dPnt_t map(const Point_t &p)
  {
    unmapped = p.z();
    return { p.x(), p.y() };
  }
  Point_t unmap(const Geo2dPnt_t &p)
  {
    return { p.x(), p.y(), unmapped };
  }
};

struct MapXZ : Map2d
{
  Geo2dPnt_t map(const Point_t &p)
  {
    unmapped = p.y();
    return { p.x(), p.z() };
  }
  Point_t unmap(const Geo2dPnt_t &p)
  {
    return { p.x(), unmapped, p.y() };
  }
};

struct MapYZ : Map2d
{
  Geo2dPnt_t map(const Point_t &p)
  {
    unmapped = p.x();
    return { p.z(), p.y() };
  }
  Point_t unmap(const Geo2dPnt_t &p)
  {
    return { unmapped, p.y(), p.x() };
  }
};

ostream &operator<<(ostream &os, const Point_t &p)
{
  os << fixed << setprecision(3)
     << setw(8) << p.x() << ", "
     << setw(8) << p.y() << ", "
     << setw(8) << p.z();
  return os;
}

template<size_t N>
ostream &operator<<(ostream &os, const Points_t<N> &points)
{
  for(const Point_t &pnt : points)
    os << pnt << endl;
  return os;
}

template<size_t N>
void transformInPlace(Points_t<N> &points, const function<Point_t(const Point_t&)> &transformPoint)
{
  ranges::transform(points, points.begin(), transformPoint);
}

class Test
{
  const RefPoints_t refPoints
  {{
    { -400, 100, 1000 },
    {  400, 100, 1000 },
    { -400, 900, 1000 }
  }};

  using Cuboid_t = Points_t<8>;
  const Cuboid_t cube
  {{ // center 0, 500, 500; edge length 600
    { -300, 200, 200 },
    {  300, 200, 200 },
    { -300, 200, 800 },
    {  300, 200, 800 },
    { -300, 800, 200 },
    {  300, 800, 200 },
    { -300, 800, 800 },
    {  300, 800, 800 }
  }};

  RefPoints_t refPointsT = refPoints;
  Cuboid_t cubeT = cube;

  void translateAll(Coordinate_t x, Coordinate_t y, Coordinate_t z)
  {
    const geo::strategy::transform::translate_transformer<Coordinate_t, 3, 3> transl(x, y, z);
    const auto geoTransform = [&transl](const Point_t &point)
                              {
                                const GeoPnt_t gp(point.x(), point.y(), point.z());
                                GeoPnt_t gpT;
                                geo::transform(gp, gpT, transl);
                                return Point_t{ gpT.x(), gpT.y(), gpT.z() };
                              };
    transformInPlace(refPointsT, geoTransform);
    transformInPlace(cubeT, geoTransform);
  }

  template<typename Map2d>
  void rotateAll(const Point_t &rotPoint, Map2d map2d, Coordinate_t degrees)
  {
    const geo::strategy::transform::rotate_transformer<geo::degree, Coordinate_t, 2, 2> rot(-degrees);
    const auto geoTransform = [&](const Point_t &point)
                              {
                                const Geo2dPnt_t gp = map2d.map(point - rotPoint);
                                Geo2dPnt_t gpT;
                                geo::transform(gp, gpT, rot);
                                return map2d.unmap(gpT) + rotPoint;
                              };
    transformInPlace(refPointsT, geoTransform);
    transformInPlace(cubeT, geoTransform);
  }

  void calcTransformedPoints()
  {
    translateAll(0, -((refPoints[0].y() + refPoints[2].y()) / 2), 0);

    const Point_t rotPoint{ 0, 0, refPoints[0].z() };
    const Coordinate_t degreesX = 50;
    rotateAll(rotPoint, MapYZ(), degreesX);

    const Coordinate_t degreesY = 35;
    rotateAll(rotPoint, MapXZ(), degreesY);

    const Coordinate_t degreesZ = 20;
    rotateAll(rotPoint, MapXY(), degreesZ);

    translateAll(-123.4, 345.6, -98.7);
  }

public:
  void test()
  {
    calcTransformedPoints();

    cout << "transformed reference points:" << endl;
    cout << refPointsT << endl;
    cout << "transformed object points:" << endl;
    cout << cubeT << endl;

    const Transformation trans(refPoints, refPointsT);
    transformInPlace(cubeT, trans);

    cout << "back-transformed object points:" << endl;
    cout << cubeT;
  }
};

int main()
{
  Test t;
  t.test();

  return 0;
}
