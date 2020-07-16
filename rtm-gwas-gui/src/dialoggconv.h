#ifndef DIALOGGCONV_H
#define DIALOGGCONV_H

#include <QDialog>
#include <QStringList>

class QCheckBox;
class QComboBox;
class QLineEdit;

class DialogGConv : public QDialog
{
    Q_OBJECT

public:
    DialogGConv(QWidget *parent = nullptr);
    QString program() const;
    QStringList arguments() const;

private:
    void setup_layout();
    void setup_group_box();
    void setup_button_box();

private slots:
    void get_input();
    void restore_defaults();

private:
    QLineEdit *in_file_;
    QLineEdit *out_file_;
    QComboBox *in_format_;
    QComboBox *out_format_;
    QCheckBox *sort_;
};

#endif // DIALOGGCONV_H
