#ifndef DIALOGLD_H
#define DIALOGLD_H

#include <QDialog>
#include <QStringList>

class QSpinBox;
class QLineEdit;
class QDoubleSpinBox;

class DialogLD : public QDialog
{
    Q_OBJECT

public:
    DialogLD(QWidget *parent = nullptr);
    QString program() const;
    QStringList arguments() const;

private:
    void setup_layout();
    void setup_group_box_data();
    void setup_group_box_option();
    void setup_button_box();

private slots:
    void get_vcf();
    void get_loc();
    void restore_defaults();

private:
    QLineEdit *vcf_;
    QLineEdit *loc_;
    QLineEdit *out_;
    QSpinBox *maxdist_;
    QDoubleSpinBox *minr2_;
};

#endif // DIALOGLD_H
