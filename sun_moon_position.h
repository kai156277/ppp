#ifndef SUN_MOON_POSITION_H
#define SUN_MOON_POSITION_H

#include "time_system.h"
#include "coordinate_system.h"
#include <istream>
#include <cmath>

using namespace std;
class sun_moon_position
{
public:
    sun_moon_position();
    static XYZ_coordinate sunPostion(const julian_Time &Julian);
    static XYZ_coordinate moonPostion(const julian_Time &Julian);
    static const double min_julian;
    static const double max_julian;
private:
    static XYZ_coordinate getSunPostionCIS(const julian_Time &Julian);
    static XYZ_coordinate getMoonPostionCIS(const julian_Time &Julian);
    static XYZ_coordinate CIStoCTS(const XYZ_coordinate &posCIS,const julian_Time &Juliday);
};

#endif // SUN_MOON_POSITION_H
