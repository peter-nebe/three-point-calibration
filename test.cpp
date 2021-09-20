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
#include <functional>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numbers>
using namespace std;
using namespace threePointCalibration;
namespace geo = boost::geometry;

typedef geo::model::d3::point_xyz<Coordinate_t> GeoPnt3_t;
typedef geo::model::d2::point_xy<Coordinate_t> GeoPnt2_t;

struct Map2d
{
protected:
  Coordinate_t unmapped = 0;
};

struct MapXY : Map2d
{
  GeoPnt2_t map(const Point3 &p)
  {
    unmapped = p.z;
    return { p.x, p.y };
  }
  Point3 unmap(const GeoPnt2_t &p)
  {
    return { p.x(), p.y(), unmapped };
  }
};

struct MapXZ : Map2d
{
  GeoPnt2_t map(const Point3 &p)
  {
    unmapped = p.y;
    return { p.x, p.z };
  }
  Point3 unmap(const GeoPnt2_t &p)
  {
    return { p.x(), unmapped, p.y() };
  }
};

struct MapYZ : Map2d
{
  GeoPnt2_t map(const Point3 &p)
  {
    unmapped = p.x;
    return { p.z, p.y };
  }
  Point3 unmap(const GeoPnt2_t &p)
  {
    return { unmapped, p.y(), p.x() };
  }
};

ostream &operator<<(ostream &os, const Point3 &p)
{
  os << fixed << setprecision(3)
     << setw(8) << p.x << ", "
     << setw(8) << p.y << ", "
     << setw(8) << p.z;
  return os;
}

template<size_t N>
ostream &operator<<(ostream &os, const array<Point3, N> &points)
{
  for(const Point3 &p : points)
    os << p << endl;
  return os;
}

template<typename Points>
void transformInPlace(Points &points, const function<Point3(const Point3&)> &transformPoint)
{
  ranges::transform(points, points.begin(), transformPoint);
}

class Test
{
protected:
  using RefPoints = ReferencePoints_<3>;

  const RefPoints _refPoints
  {{
    // x: to the right
    // y: height
    // z: forward
    { -400, 100, 1000 },
    {  400, 100, 1000 },
    { -400, 900, 1000 }
  }};

  RefPoints _refPointsMapping = _refPoints;

  array<Point3, 8> _cube
  {{
    // cube: center 0, 500, 500; edge length 600
    { -300, 200, 200 },
    {  300, 200, 200 },
    { -300, 200, 800 },
    {  300, 200, 800 },
    { -300, 800, 200 },
    {  300, 800, 200 },
    { -300, 800, 800 },
    {  300, 800, 800 }
  }};

  void translateAll(Coordinate_t x, Coordinate_t y, Coordinate_t z = 0)
  {
    const geo::strategy::transform::translate_transformer<Coordinate_t, 3, 3> transl(x, y, z);
    const auto geoTransform = [&transl](const Point3 &point)
                              {
                                const GeoPnt3_t gp(point.x, point.y, point.z);
                                GeoPnt3_t gpT;
                                geo::transform(gp, gpT, transl);

                                return Point3{ gpT.x(), gpT.y(), gpT.z() };
                              };
    transformInPlace(_refPointsMapping, geoTransform);
    transformInPlace(_cube, geoTransform);
  }

  template<typename Map2d>
  void rotateAll(const Point3 &rotPoint, Map2d map2d, Coordinate_t degrees)
  {
    const geo::strategy::transform::rotate_transformer<geo::degree, Coordinate_t, 2, 2> rot(-degrees);
    const auto geoTransform = [&](const Point3 &point)
                              {
                                using boost::qvm::operator-;
                                using boost::qvm::operator+;

                                const GeoPnt2_t gp = map2d.map(point - rotPoint);
                                GeoPnt2_t gpT;
                                geo::transform(gp, gpT, rot);

                                return map2d.unmap(gpT) + rotPoint;
                              };
    transformInPlace(_refPointsMapping, geoTransform);
    transformInPlace(_cube, geoTransform);
  }

private:
  // transform into camera coordinates
  virtual void calcMappedPoints()
  {
    translateAll(0, -((_refPoints[0].y + _refPoints[2].y) / 2), 0);

    const Point3 rotPoint{ 0, 0, _refPoints[0].z };
    const Coordinate_t degreesX = 50;
    rotateAll(rotPoint, MapYZ(), degreesX);

    const Coordinate_t degreesY = 35;
    rotateAll(rotPoint, MapXZ(), degreesY);

    const Coordinate_t degreesZ = 20;
    rotateAll(rotPoint, MapXY(), degreesZ);

    translateAll(-123.4, 345.6, -98.7);
  }

  // transform into world coordinates
  virtual void transformToWorld()
  {
    const Transformation trans(_refPoints, _refPointsMapping);
    transformInPlace(_cube, trans);
  }

public:
  void test()
  {
    calcMappedPoints();

    cout << "mapping of the reference points:" << endl;
    cout << _refPointsMapping << endl;
    cout << "mapping of the object points:" << endl;
    cout << _cube << endl;

    transformToWorld();

    cout << "object points transformed into world coordinates:" << endl;
    cout << _cube;
  }
}; // class Test

class TestSimple : public Test
{
  const RefPoints _triangleOnTheFloor
  {{
    // x: to the right
    // y: forward
    // z: height
    { -500, 1000, 0 },
    {  500, 1000, 0 },
    {  400,  100, 0 }
  }};
  const Point3 _objPoint = _cube.front();

  void calcMappedPoints() override
  {
    translateAll(0, 0, -1000);

    const Point3 rotPoint{};
    const Coordinate_t degreesX = -140;
    rotateAll(rotPoint, MapYZ(), degreesX);

    const Coordinate_t degreesY = 35;
    rotateAll(rotPoint, MapXZ(), degreesY);

    const Coordinate_t degreesZ = 20;
    rotateAll(rotPoint, MapXY(), degreesZ);
  }

  void transformToWorld() override
  {
    // Use only three points in camera coordinates. They determine the xy plane
    // in world coordinates. The position of the points within the plane is
    // irrelevant, they just need to form a large triangle. The world coordinate
    // origin is determined by the foot of the perpendicular from the camera to
    // the xy plane.
    const Transformation trans(_refPointsMapping);
    transformInPlace(_cube, trans);

    // The cube object has now been transformed into world coordinates. The
    // z-coordinates of all points have their original value (see above).
    // However, the xy values appear rotated. Because the direction of the
    // y-axis within the xy-plane is determined by the viewing direction of
    // the camera. We now correct this rotation. The angle is determined using
    // the first object point.
    const Point3 objPointMapping = _cube.front();
    auto atan = [](const Point3 &p)
                {
                  return std::atan2(p.y, p.x);
                };
    const Coordinate_t angle = (atan(objPointMapping) - atan(_objPoint)) * 180 / numbers::pi;

    rotateAll(Point3{}, MapXY(), -angle);
  }

public:
  TestSimple()
  {
    _refPointsMapping = _triangleOnTheFloor;
  }
}; // class TestSimple

class Test2D : public Test
{
  const RefPoints _refPoints
  {{
    { -500, 1000 },
    {  500, 1000 }
  }};

  void calcMappedPoints() override
  {
    translateAll(700, -800);

    const Point3 rotPoint{};
    const Coordinate_t degreesZ = 50;
    rotateAll(rotPoint, MapXY(), degreesZ);

    translateAll(-123, 45);
  }

  void transformToWorld() override
  {
    struct P2 : Point2
    {
      P2(const Point3 &p) : Point2{ p.x, p.y } {}
    };
    struct P3 : Point3
    {
      P3(const Point2 &p, Coordinate_t z) : Point3{ p.x, p.y, z } {}
    };
    const Transformation2D trans({
                                   P2(_refPoints[0]),
                                   P2(_refPoints[1])
                                 },
                                 {
                                   P2(_refPointsMapping[0]),
                                   P2(_refPointsMapping[1])
                                 });
    transformInPlace(_cube, [&trans](const Point3 &p)
                            {
                              return P3(trans(P2(p)), p.z);
                            });
  }

public:
  Test2D()
  {
    _refPointsMapping = _refPoints;
  }
}; // class Test2D

int main()
{
  Test().test();

  cout << endl << "---simplified calibration with three points on the floor---" << endl;
  TestSimple().test();

  cout << endl << "---2D calibration with two points---" << endl;
  Test2D().test();

  return 0;
}
