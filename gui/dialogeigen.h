#ifndef DIALOGEIGEN_H
#define DIALOGEIGEN_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogEigen;
}

class DialogEigen : public QDialog
{
    Q_OBJECT

public:
    explicit DialogEigen(QWidget *parent = 0);

    ~DialogEigen();

    QStringList arguments() const;

public slots:
    void apply();

private slots:
    void on_pushButtonGenotype_clicked();

private:
    Ui::DialogEigen *ui;
};

#endif // DIALOGEIGEN_H
