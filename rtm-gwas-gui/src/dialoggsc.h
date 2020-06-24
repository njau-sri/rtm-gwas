#ifndef DIALOGGSC_H
#define DIALOGGSC_H

#include <QDialog>
#include <QStringList>

class QSpinBox;
class QLineEdit;

class DialogGSC : public QDialog
{
    Q_OBJECT

public:
    DialogGSC(QWidget *parent = nullptr);
    QString program() const;
    QStringList arguments() const;

private:
    void setup_layout();
    void setup_group_box_data();
    void setup_group_box_option();
    void setup_button_box();

private slots:
    void get_vcf();
    void get_grm();
    void restore_defaults();

private:
    QLineEdit *vcf_;
    QLineEdit *grm_;
    QLineEdit *out_;
    QSpinBox *top_;
    QSpinBox *thread_;
};

#endif // DIALOGGSC_H
