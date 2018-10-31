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

    QString getProg() const;

    QStringList getArgs() const;

public slots:

    void apply();

private slots:

    void on_pushButtonVCF_clicked();

    void on_pushButtonPheno_clicked();

    void on_pushButtonCovar_clicked();

private:

    Ui::DialogAssoc *ui;
};

#endif // DIALOGASSOC_H
