#ifndef DIALOGASSOC_H
#define DIALOGASSOC_H

#include <QDialog>
#include <QStringList>

class QSpinBox;
class QComboBox;
class QLineEdit;

class DialogAssoc : public QDialog
{
    Q_OBJECT

public:
    DialogAssoc(QWidget *parent = nullptr);
    QString program() const;
    QStringList arguments() const;

private:
    void setup_layout();
    void setup_group_box_data();
    void setup_group_box_option();
    void setup_button_box();

private slots:
    void get_vcf();
    void get_pheno();
    void get_covar();
    void restore_defaults();

private:
    QLineEdit *vcf_;
    QLineEdit *pheno_;
    QLineEdit *covar_;
    QLineEdit *out_;
    QLineEdit *alpha_;
    QLineEdit *preselect_;
    QLineEdit *rsq_;
    QComboBox *mtc_;
    QSpinBox *thread_;
};

#endif // DIALOGASSOC_H
