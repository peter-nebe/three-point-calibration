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
#include <numbers>
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

template<size_t N>
void transformInvInPlace(Points_t<N> &points, const Transformation &trans)
{
  ranges::transform(points, points.begin(), [&trans](const Point_t &p){ return trans.transformInv(p); });
}

class Test
{
  const RefPoints_t _refPoints
  {{
    // x: to the right
    // y: height
    // z: forward
    { -400, 100, 1000 },
    {  400, 100, 1000 },
    { -400, 900, 1000 }
  }};

protected:
  RefPoints_t _refPointsT = _refPoints;

  Points_t<8> _cubeT
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
    transformInPlace(_refPointsT, geoTransform);
    transformInPlace(_cubeT, geoTransform);
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
    transformInPlace(_refPointsT, geoTransform);
    transformInPlace(_cubeT, geoTransform);
  }

private:
  // transform to camera coordinates
  virtual void calcTransformedPoints()
  {
    translateAll(0, -((_refPoints[0].y() + _refPoints[2].y()) / 2), 0);

    const Point_t rotPoint{ 0, 0, _refPoints[0].z() };
    const Coordinate_t degreesX = 50;
    rotateAll(rotPoint, MapYZ(), degreesX);

    const Coordinate_t degreesY = 35;
    rotateAll(rotPoint, MapXZ(), degreesY);

    const Coordinate_t degreesZ = 20;
    rotateAll(rotPoint, MapXY(), degreesZ);

    translateAll(-123.4, 345.6, -98.7);
  }

  // transform to world coordinates
  virtual void transformBack()
  {
    const Transformation trans(_refPoints, _refPointsT);
    transformInPlace(_cubeT, trans);
  }

public:
  void test()
  {
    calcTransformedPoints();

    cout << "transformed reference points:" << endl;
    cout << _refPointsT << endl;
    cout << "transformed object points:" << endl;
    cout << _cubeT << endl;

    transformBack();

    cout << "back-transformed object points:" << endl;
    cout << _cubeT;
  }
}; // class Test

class TestSimple : public Test
{
  const RefPoints_t _triangleOnTheFloor
  {{
    // x: to the right
    // y: forward
    // z: height
    { -500, 1000, 0 },
    {  500, 1000, 0 },
    {  400,  100, 0 }
  }};
  const Point_t _objPoint = _cubeT.front();

  void calcTransformedPoints() override
  {
    translateAll(0, 0, -1000);

    const Point_t rotPoint{};
    const Coordinate_t degreesX = -140;
    rotateAll(rotPoint, MapYZ(), degreesX);

    const Coordinate_t degreesY = 35;
    rotateAll(rotPoint, MapXZ(), degreesY);

    const Coordinate_t degreesZ = 20;
    rotateAll(rotPoint, MapXY(), degreesZ);
  }

  void transformBack() override
  {
    // Use only three points in camera coordinates. They determine the xy plane
    // in world coordinates. The position of the points within the plane is
    // irrelevant, they just need to form a large triangle. The world coordinate
    // origin is determined by the foot of the perpendicular from the camera to
    // the xy plane.
    const Transformation trans(_refPointsT);
    transformInvInPlace(_cubeT, trans);

    // The cube object has now been transformed into world coordinates. The
    // z-coordinates of all points have their original value (see above).
    // However, the xy values appear rotated. Because the direction of the
    // y-axis within the xy-plane is determined by the viewing direction of
    // the camera. We now correct this rotation. The angle is determined using
    // the first object point.
    const Point_t objPointT = _cubeT.front();
    const Coordinate_t angle = (atan2(objPointT.y(), objPointT.x()) - atan2(_objPoint.y(), _objPoint.x())) * 180 / numbers::pi;
    rotateAll(Point_t{}, MapXY(), -angle);
  }

public:
  TestSimple()
  {
    _refPointsT = _triangleOnTheFloor;
  }
}; // class TestSimple

int main()
{
  Test().test();

  cout << endl << "---simplified calibration with triangle on the floor---" << endl;
  TestSimple().test();

  return 0;
}
