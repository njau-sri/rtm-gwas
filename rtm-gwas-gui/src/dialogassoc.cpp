#include "dialogassoc.h"

#include <QLabel>
#include <QSpinBox>
#include <QComboBox>
#include <QDateTime>
#include <QGroupBox>
#include <QLineEdit>
#include <QCheckBox>
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

DialogAssoc::DialogAssoc(QWidget *parent) : QDialog(parent)
{
    setWindowTitle("Association - RTM-GWAS");

    setup_layout();
    setup_group_box_data();
    setup_group_box_option();
    setup_button_box();

    restore_defaults();
}

QString DialogAssoc::program() const
{
    return QDir(par->bin_path).absoluteFilePath("rtm-gwas-assoc");
}

QStringList DialogAssoc::arguments() const
{
    QStringList args;

    if (!vcf_->text().isEmpty())
        args << "--vcf" << vcf_->text();

    if (!pheno_->text().isEmpty())
        args << "--pheno" << pheno_->text();

    if (!covar_->text().isEmpty())
        args << "--covar" << covar_->text();

    if (!alpha_->text().isEmpty())
        args << "--alpha" << alpha_->text();

    if (!preselect_->text().isEmpty())
        args << "--preselect" << preselect_->text();

    if (!rsq_->text().isEmpty())
        args << "--rsq" << rsq_->text();

    if (mtc_->currentIndex() != 0)
        args << "--mtc" << mtc_->currentText();

    args << "--thread" << thread_->cleanText();

    if (nogxe_->isChecked())
        args << "--no-gxe";

    QString prefix = out_->text();
    if (prefix.isEmpty())
        prefix = "assoc.out." + current_datetime_string();
    args << "--out" << QDir(par->working_directory).absoluteFilePath(prefix);

    return args;
}

void DialogAssoc::setup_layout()
{
    setLayout(new QVBoxLayout);
}

void DialogAssoc::setup_group_box_data()
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

    QLabel *phenol = new QLabel("Phenotype data:", box);
    pheno_ = new QLineEdit(box);
    pheno_->setPlaceholderText("required");
    QPushButton *phenob = new QPushButton("Browse...", box);
    grid->addWidget(phenol, 1, 0);
    grid->addWidget(pheno_, 1, 1);
    grid->addWidget(phenob, 1, 2);
    connect(phenob, SIGNAL(clicked()), this, SLOT(get_pheno()));

    QLabel *covarl = new QLabel("Covariate data:", box);
    covar_ = new QLineEdit(box);
    covar_->setPlaceholderText("optional");
    QPushButton *covarb = new QPushButton("Browse...", box);
    grid->addWidget(covarl, 2, 0);
    grid->addWidget(covar_, 2, 1);
    grid->addWidget(covarb, 2, 2);
    connect(covarb, SIGNAL(clicked()), this, SLOT(get_covar()));

    QLabel *outl = new QLabel("Output file prefix:", box);
    out_ = new QLineEdit(box);
    out_->setPlaceholderText("required");
    out_->setText("assoc.out." + current_datetime_string());
    grid->addWidget(outl, 3, 0);
    grid->addWidget(out_, 3, 1);
}

void DialogAssoc::setup_group_box_option()
{
    QGroupBox *box = new QGroupBox("Option", this);
    layout()->addWidget(box);

    QFormLayout *form = new QFormLayout;
    box->setLayout(form);

    alpha_ = new QLineEdit(box);
    alpha_->setPlaceholderText("0.01");
    alpha_->setToolTip("value range (0, 1)");
    form->addRow("Significance level:", alpha_);

    preselect_ = new QLineEdit(box);
    preselect_->setPlaceholderText("0.05");
    preselect_->setToolTip("value range (0, 1)");
    form->addRow("Pre-selection threshold:", preselect_);

    rsq_ = new QLineEdit(box);
    rsq_->setPlaceholderText("0.95");
    rsq_->setToolTip("value range [0, 1]");
    form->addRow("Maximum model r-square:", rsq_);

    mtc_ = new QComboBox(box);
    mtc_->addItems(QStringList() << "None" << "BON" << "FDR");
    form->addRow("Multiple testing correction:", mtc_);

    thread_ = new QSpinBox(box);
    form->addRow("Number of threads:", thread_);

    nogxe_ = new QCheckBox("Genotype-environment interaction", box);
    form->addRow("", nogxe_);
}

void DialogAssoc::setup_button_box()
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

void DialogAssoc::get_vcf()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        vcf_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogAssoc::get_pheno()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        pheno_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogAssoc::get_covar()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        covar_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogAssoc::restore_defaults()
{
    vcf_->clear();
    pheno_->clear();
    covar_->clear();
    out_->setText("assoc.out." + current_datetime_string());
    alpha_->setText("0.01");
    preselect_->setText("0.05");
    rsq_->setText("0.95");
    mtc_->setCurrentIndex(0);
    thread_->setValue(0);
    nogxe_->setChecked(false);
}
