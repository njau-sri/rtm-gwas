#ifndef DIALOGLDB_H
#define DIALOGLDB_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogLDB;
}

class DialogLDB : public QDialog
{
    Q_OBJECT

public:
    explicit DialogLDB(QWidget *parent = 0);

    ~DialogLDB();

    QStringList arguments() const;

public slots:
    void apply();

private slots:
    void on_pushButtonGenotype_clicked();

    void on_pushButtonBlock_clicked();

private:
    Ui::DialogLDB *ui;
};

#endif // DIALOGLDB_H
