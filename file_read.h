#ifndef FILE_READ_H
#define FILE_READ_H

#include<QString>

#include"o_date.h"
#include"sp3_date.h"

class file_read
{
public:
    file_read();
    void ppp_o_read(QString file_path,o_file_date &obs);
    void ppp_sp3_read(QString file_path,sp3_file &sp3);
private:
    void phase_matching(const QVector<sys_record> &match_list, system_signal &sys_list);
    int each_phase_matching(const QStringList & reference_list,const QStringList & obs_descriptor);
};

#endif // FILE_READ_H
