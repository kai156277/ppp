#ifndef TIDE_CORRECTION_H
#define TIDE_CORRECTION_H

#include "ppp_data.h"
#include "coordinate_system.h"
#include "time_system.h"
#include "ocean_date.h"
#include "erp_data.h"
#include <Eigen/Eigen>

using namespace Eigen;

class tide_correction
{
public:
    tide_correction();
    static Vector3d solidTide(const XYZ_coordinate &receiver,const XYZ_coordinate &sunPos,const XYZ_coordinate &moonPos,GPS_Time t);
    static Vector3d oceanTide(const ocean &ocean_data, GPS_Time t);
    static Vector3d poleTide(const erp &erp_data, XYZ_coordinate receiver, GPS_Time t);

private:
    static Vector3d ST1DIU(const Vector3d &XSTA, const Vector3d &XSUN, const Vector3d &XMON, double F2SUN, double F2MON);
    static Vector3d ST1SEM(const Vector3d &XSTA, const Vector3d &XSUN, const Vector3d &XMON, double F2SUN, double F2MON);
    static Vector3d ST1L1(const Vector3d &XSTA, const Vector3d &XSUN, const Vector3d &XMON, double F2SUN, double F2MON);
    static Vector3d STEP2DIU(const Vector3d &XSTA, double FHR, double T);
    static Vector3d ST2LON(const Vector3d &XSTA, double T);
    static Vector3d ST2DIU(const Vector3d &XSTA, double FHR, double T);
    static MatrixXd getArg(GPS_Time t);
    static const double Pi;

};

#endif // TIDE_CORRECTION_H
