#ifndef PARAMETER_H
#define PARAMETER_H

#include <QString>

struct Parameter
{
    QString working_directory;
    QString file_dialog_directory;
    QString bin_path;

    int txt_size;
    int log_size;
};

extern Parameter *par;

#endif // PARAMETER_H
