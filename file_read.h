#ifndef FILE_READ_H
#define FILE_READ_H

#include<QString>

#include"o_date.h"
#include"sp3_date.h"
#include"clock_date.h"
#include"snx_date.h"
#include"antmod_date.h"

class file_read
{
public:
    file_read();
    void ppp_o_read(const QString &file_path,o_file_date &obs);
    void ppp_sp3_read(const QString &file_path,sp3_file &sp3);
    void ppp_clock_read(const QString &file_path,clock_file &clock);
    void ppp_snx_read(const QString &file_path,snx_date &snx,QString mark_name);
    void ppp_ant_read(const QString &file_path,antmod_file &ant);
private:
    void phase_matching(const QStringList &match_list, system_signal &sys_list);//选择观测信号
    int each_phase_matching(const QStringList & reference_list,const QStringList & obs_descriptor);//匹配最精确的观测信号
};

#endif // FILE_READ_H
