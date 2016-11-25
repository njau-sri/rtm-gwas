#ifndef DIALOGDATA_H
#define DIALOGDATA_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogData;
}

class DialogData : public QDialog
{
    Q_OBJECT

public:
    explicit DialogData(QWidget *parent = 0);

    ~DialogData();

    QStringList arguments() const;

public slots:
    void apply();

private slots:
    void on_pushButtonGenotype_clicked();

    void on_pushButtonLocus_clicked();

    void on_pushButtonIndividual_clicked();

private:
    Ui::DialogData *ui;
};

#endif // DIALOGDATA_H
