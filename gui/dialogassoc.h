#ifndef DIALOGASSOC_H
#define DIALOGASSOC_H

#include <QDialog>
#include <QStringList>

namespace Ui {
class DialogAssoc;
}

class DialogAssoc : public QDialog
{
    Q_OBJECT

public:
    explicit DialogAssoc(QWidget *parent = 0);

    ~DialogAssoc();

    QStringList arguments() const;

public slots:
    void apply();

private slots:
    void on_pushButtonGenotype_clicked();

    void on_pushButtonPhenotype_clicked();

    void on_pushButtonCovariate_clicked();

    void on_pushButtonKinship_clicked();

private:
    Ui::DialogAssoc *ui;
};

#endif // DIALOGASSOC_H
