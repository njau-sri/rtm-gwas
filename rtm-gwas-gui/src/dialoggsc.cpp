#include "dialoggsc.h"

#include <QLabel>
#include <QSpinBox>
#include <QDateTime>
#include <QGroupBox>
#include <QLineEdit>
#include <QFileDialog>
#include <QFormLayout>
#include <QGridLayout>
#include <QPushButton>
#include <QSpacerItem>
#include <QVBoxLayout>
#include <QDialogButtonBox>

#include "parameter.h"

namespace {

QString current_datetime_string()
{
    return QDateTime::currentDateTime().toString("yyMMddhhmmsszzz");
}

} // namespace

DialogGSC::DialogGSC(QWidget *parent) : QDialog(parent)
{
    setWindowTitle("GSC - RTM-GWAS");

    setup_layout();
    setup_group_box_data();
    setup_group_box_option();
    setup_button_box();

    restore_defaults();
}

QString DialogGSC::program() const
{
    return QDir(par->bin_path).absoluteFilePath("rtm-gwas-gsc");
}

QStringList DialogGSC::arguments() const
{
    QStringList args;

    if (!vcf_->text().isEmpty())
        args << "--vcf" << vcf_->text();

    if (!grm_->text().isEmpty())
        args << "--grm" << grm_->text();

    args << "--top" << top_->cleanText();
    args << "--thread" << thread_->cleanText();

    QString prefix = out_->text();
    if (prefix.isEmpty())
        prefix = "gsc.out." + current_datetime_string();
    args << "--out" << QDir(par->working_directory).absoluteFilePath(prefix);

    return args;
}

void DialogGSC::setup_layout()
{
    setLayout(new QVBoxLayout);
}

void DialogGSC::setup_group_box_data()
{
    QGroupBox *box = new QGroupBox("Data", this);
    layout()->addWidget(box);

    QGridLayout *grid = new QGridLayout;
    box->setLayout(grid);

    QLabel *vcfl = new QLabel("Genotype data:", box);
    vcf_ = new QLineEdit(box);
    vcf_->setPlaceholderText("required, VCF format");
    QPushButton *vcfb = new QPushButton("Browse...", box);
    grid->addWidget(vcfl, 0, 0);
    grid->addWidget(vcf_, 0, 1);
    grid->addWidget(vcfb, 0, 2);
    connect(vcfb, SIGNAL(clicked()), this, SLOT(get_vcf()));

    QLabel *grml = new QLabel("Genetic relationship matrix:", box);
    grm_ = new QLineEdit(box);
    grm_->setPlaceholderText("optional");
    QPushButton *grmb = new QPushButton("Browse...", box);
    grid->addWidget(grml, 1, 0);
    grid->addWidget(grm_, 1, 1);
    grid->addWidget(grmb, 1, 2);
    connect(grmb, SIGNAL(clicked()), this, SLOT(get_grm()));

    QLabel *outl = new QLabel("Output file prefix:", box);
    out_ = new QLineEdit(box);
    out_->setPlaceholderText("required");
    out_->setText("gsc.out." + current_datetime_string());
    grid->addWidget(outl, 2, 0);
    grid->addWidget(out_, 2, 1);
}

void DialogGSC::setup_group_box_option()
{
    QGroupBox *box = new QGroupBox("Option", this);
    layout()->addWidget(box);

    QFormLayout *form = new QFormLayout;
    box->setLayout(form);

    top_ = new QSpinBox(box);
    top_->setRange(1, 9999999);
    form->addRow("Number of eigenvectors:", top_);

    thread_ = new QSpinBox(box);
    form->addRow("Number of threads:", thread_);
}

void DialogGSC::setup_button_box()
{
    QSpacerItem *spacer = new QSpacerItem(1, 1, QSizePolicy::Minimum, QSizePolicy::Expanding);
    layout()->addItem(spacer);

    QDialogButtonBox *btn = new QDialogButtonBox(this);
    btn->setOrientation(Qt::Horizontal);
    btn->setStandardButtons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel | QDialogButtonBox::Reset);
    layout()->addWidget(btn);

    QPushButton *btn_reset = btn->button(QDialogButtonBox::Reset);

    connect(btn, SIGNAL(accepted()), this, SLOT(accept()));
    connect(btn, SIGNAL(rejected()), this, SLOT(reject()));
    connect(btn_reset, SIGNAL(clicked()), this, SLOT(restore_defaults()));
}

void DialogGSC::get_vcf()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        vcf_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogGSC::get_grm()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        grm_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogGSC::restore_defaults()
{
    vcf_->clear();
    grm_->clear();
    out_->setText("gsc.out." + current_datetime_string());
    top_->setValue(10);
    thread_->setValue(0);
}
