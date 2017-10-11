#ifndef DIALOGGSC_H
#define DIALOGGSC_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogGSC;
}

class DialogGSC : public QDialog
{
    Q_OBJECT

public:

    explicit DialogGSC(QWidget *parent = 0);

    ~DialogGSC();

    QString getProg() const;

    QStringList getArgs() const;

public slots:

    void apply();

private slots:

    void on_pushButtonVcf_clicked();

private:

    Ui::DialogGSC *ui;
};

#endif // DIALOGGSC_H
