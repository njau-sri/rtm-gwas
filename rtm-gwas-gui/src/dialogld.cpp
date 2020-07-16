#include "dialogld.h"

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
#include <QDoubleSpinBox>
#include <QDialogButtonBox>

#include "parameter.h"

namespace {

QString current_datetime_string()
{
    return QDateTime::currentDateTime().toString("yyMMddhhmmsszzz");
}

} // namespace

DialogLD::DialogLD(QWidget *parent) : QDialog(parent)
{
    setWindowTitle("LD - RTM-GWAS");

    setup_layout();
    setup_group_box_data();
    setup_group_box_option();
    setup_button_box();

    restore_defaults();
}

QString DialogLD::program() const
{
    return QDir(par->bin_path).absoluteFilePath("rtm-gwas-ld");
}

QStringList DialogLD::arguments() const
{
    QStringList args;

    if (!vcf_->text().isEmpty())
        args << "--vcf" << vcf_->text();

    if (!loc_->text().isEmpty())
        args << "--with-loc" << loc_->text();

    args << "--maxdist" << maxdist_->cleanText();

    args << "--min-r2" << minr2_->cleanText();

    QString prefix = out_->text();
    if (prefix.isEmpty())
        prefix = "gsc.out." + current_datetime_string();
    args << "--out" << QDir(par->working_directory).absoluteFilePath(prefix);

    return args;
}

void DialogLD::setup_layout()
{
    setLayout(new QVBoxLayout);
}

void DialogLD::setup_group_box_data()
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

    QLabel *grml = new QLabel("Locus list file:", box);
    loc_ = new QLineEdit(box);
    loc_->setPlaceholderText("optional");
    QPushButton *grmb = new QPushButton("Browse...", box);
    grid->addWidget(grml, 1, 0);
    grid->addWidget(loc_, 1, 1);
    grid->addWidget(grmb, 1, 2);
    connect(grmb, SIGNAL(clicked()), this, SLOT(get_loc()));

    QLabel *outl = new QLabel("Output file prefix:", box);
    out_ = new QLineEdit(box);
    out_->setPlaceholderText("required");
    out_->setText("ld.out." + current_datetime_string());
    grid->addWidget(outl, 2, 0);
    grid->addWidget(out_, 2, 1);
}

void DialogLD::setup_group_box_option()
{
    QGroupBox *box = new QGroupBox("Option", this);
    layout()->addWidget(box);

    QFormLayout *form = new QFormLayout;
    box->setLayout(form);

    maxdist_ = new QSpinBox(box);
    maxdist_->setRange(2, 99999999);
    maxdist_->setSingleStep(10000);
    form->addRow("Maximum inter-locus distance:", maxdist_);

    minr2_ = new QDoubleSpinBox(box);
    minr2_->setRange(0.0, 1.0);
    minr2_->setSingleStep(0.05);
    form->addRow("Minimum r2 threshold (locus list):", minr2_);
}

void DialogLD::setup_button_box()
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

void DialogLD::get_vcf()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        vcf_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogLD::get_loc()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        loc_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogLD::restore_defaults()
{
    vcf_->clear();
    loc_->clear();
    out_->setText("ld.out." + current_datetime_string());
    maxdist_->setValue(500000);
    minr2_->setValue(0.5);
}
