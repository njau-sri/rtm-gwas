#ifndef DIALOGSNPLDB_H
#define DIALOGSNPLDB_H

#include <QDialog>
#include <QStringList>

class QSpinBox;
class QCheckBox;
class QLineEdit;
class QDoubleSpinBox;

class DialogSNPLDB : public QDialog
{
    Q_OBJECT

public:
    DialogSNPLDB(QWidget *parent = nullptr);
    QString program() const;
    QStringList arguments() const;

private:
    void setup_layout();
    void setup_group_box_data();
    void setup_group_box_option();
    void setup_button_box();

private slots:
    void get_vcf();
    void get_gene();
    void get_block();
    void restore_defaults();
    void set_nam_state(int state);

private:
    QLineEdit *vcf_;
    QLineEdit *gene_;
    QLineEdit *block_;
    QLineEdit *out_;
    QLineEdit *maf_;
    QSpinBox *maxlen_;
    QSpinBox *llim_;
    QSpinBox *ulim_;
    QSpinBox *recomb_;
    QDoubleSpinBox *inform_;
    QSpinBox *thread_;
    QLineEdit *identity_;
    QCheckBox *natchk_;
    QCheckBox *rilchk_;
    QCheckBox *namchk_;
    QLineEdit *nam_;
};

#endif // DIALOGSNPLDB_H
