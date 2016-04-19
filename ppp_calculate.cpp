#include "ppp_calculate.h"
#include "gc_gpss.h"

ppp_calculate::ppp_calculate()
{

}

void ppp_calculate::ppp_coordinate(const o_file_date &ofile, const sp3_file &sp3file, ppp_file &ppp)
{
    GC_GPSS time;
    for(int i = 0; i<ofile.satellite_file.size(); i++)
    {
        ppp_epoch epoch;
        epoch.year   = ofile.satellite_file[i].year;
        epoch.month  = ofile.satellite_file[i].month;
        epoch.day    = ofile.satellite_file[i].day;
        epoch.hour   = ofile.satellite_file[i].hour;
        epoch.minute = ofile.satellite_file[i].minute;
        epoch.second = ofile.satellite_file[i].second;
        epoch.sate_number = ofile.satellite_file[i].number_of_satellite;
        time.setGC(epoch.year,epoch.month,epoch.day,epoch.hour,epoch.minute,epoch.second);
        time.GCtoGPS();
        epoch.GPSS = time.GPSS;
        time.clear();
    }
}

