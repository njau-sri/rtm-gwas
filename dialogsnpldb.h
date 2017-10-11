#ifndef DIALOGSNPLDB_H
#define DIALOGSNPLDB_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogSNPLDB;
}

class DialogSNPLDB : public QDialog
{
    Q_OBJECT

public:

    explicit DialogSNPLDB(QWidget *parent = 0);

    ~DialogSNPLDB();

    QString getProg() const;

    QStringList getArgs() const;

public slots:

    void apply();

private slots:

    void on_pushButtonVcf_clicked();

    void on_pushButtonBlock_clicked();

private:

    Ui::DialogSNPLDB *ui;
};

#endif // DIALOGSNPLDB_H
