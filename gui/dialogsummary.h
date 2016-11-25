#ifndef DIALOGSUMMARY_H
#define DIALOGSUMMARY_H

#include <QDialog>

namespace Ui {
class DialogSummary;
}

class DialogSummary : public QDialog
{
    Q_OBJECT

public:
    explicit DialogSummary(QWidget *parent = 0);

    ~DialogSummary();

    QStringList arguments() const;

public slots:
    void apply();

private slots:
    void on_pushButtonGenotype_clicked();

private:
    Ui::DialogSummary *ui;
};

#endif // DIALOGSUMMARY_H
