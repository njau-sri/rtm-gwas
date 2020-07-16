#include "dialogsnpldb.h"

#include <QLabel>
#include <QSpinBox>
#include <QCheckBox>
#include <QDateTime>
#include <QGroupBox>
#include <QLineEdit>
#include <QFileDialog>
#include <QFormLayout>
#include <QGridLayout>
#include <QPushButton>
#include <QSpacerItem>
#include <QVBoxLayout>
#include <QButtonGroup>
#include <QDoubleSpinBox>
#include <QDialogButtonBox>

#include "parameter.h"

namespace {

QString current_datetime_string()
{
    return QDateTime::currentDateTime().toString("yyMMddhhmmsszzz");
}

} // namespace

DialogSNPLDB::DialogSNPLDB(QWidget *parent) : QDialog(parent)
{
    setWindowTitle("SNPLDB - RTM-GWAS");

    setup_layout();
    setup_group_box_data();
    setup_group_box_option();
    setup_button_box();

    restore_defaults();
}

QString DialogSNPLDB::program() const
{
    return QDir(par->bin_path).absoluteFilePath("rtm-gwas-snpldb");
}

QStringList DialogSNPLDB::arguments() const
{
    QStringList args;

    if (!vcf_->text().isEmpty())
        args << "--vcf" << vcf_->text();

    if (!gene_->text().isEmpty())
        args << "--gene" << gene_->text();

    if (!block_->text().isEmpty())
        args << "--block" << block_->text();

    if (!maf_->text().isEmpty())
        args << "--maf" << maf_->text();

    args << "--maxlen" << maxlen_->cleanText();
    args << "--llim" << llim_->cleanText();
    args << "--ulim" << ulim_->cleanText();
    args << "--recomb" << recomb_->cleanText();
    args << "--inform" << inform_->cleanText();
    args << "--thread" << thread_->cleanText();

    if (!identity_->text().isEmpty())
        args << "--identity" << identity_->text();

    if (rilchk_->isChecked())
        args << "--ril";
    else if (namchk_->isChecked()) {
        if (nam_->text().isEmpty())
            args << "--nam";
        else
            args << "--nam" << nam_->text();
    }

    QString prefix = out_->text();
    if (prefix.isEmpty())
        prefix = "snpldb.out." + current_datetime_string();
    args << "--out" << QDir(par->working_directory).absoluteFilePath(prefix);

    return args;
}

void DialogSNPLDB::setup_layout()
{
    setLayout(new QVBoxLayout);
}

void DialogSNPLDB::setup_group_box_data()
{
    QGroupBox *box = new QGroupBox("Data", this);
    layout()->addWidget(box);

    QGridLayout *grid = new QGridLayout;
    box->setLayout(grid);

    QLabel *vcfl = new QLabel("SNP genotype data:", box);
    vcf_ = new QLineEdit(box);
    vcf_->setPlaceholderText("required, VCF foramt");
    QPushButton *vcfb = new QPushButton("Browse...", box);
    grid->addWidget(vcfl, 0, 0);
    grid->addWidget(vcf_, 0, 1);
    grid->addWidget(vcfb, 0, 2);
    connect(vcfb, SIGNAL(clicked()), this, SLOT(get_vcf()));

    QLabel *genel = new QLabel("Gene coordinate:", box);
    gene_ = new QLineEdit(box);
    gene_->setPlaceholderText("optional");
    QPushButton *geneb = new QPushButton("Browse...", box);
    grid->addWidget(genel, 1, 0);
    grid->addWidget(gene_, 1, 1);
    grid->addWidget(geneb, 1, 2);
    connect(geneb, SIGNAL(clicked()), this, SLOT(get_gene()));

    QLabel *blockl = new QLabel("Reference haplotype block:", box);
    block_ = new QLineEdit(box);
    block_->setPlaceholderText("optional");
    QPushButton *blockb = new QPushButton("Browse...", box);
    grid->addWidget(blockl, 2, 0);
    grid->addWidget(block_, 2, 1);
    grid->addWidget(blockb, 2, 2);
    connect(blockb, SIGNAL(clicked()), this, SLOT(get_block()));

    QLabel *outl = new QLabel("Output file prefix:", box);
    out_ = new QLineEdit(box);
    out_->setPlaceholderText("required");
    out_->setText("snpldb.out." + current_datetime_string());
    grid->addWidget(outl, 3, 0);
    grid->addWidget(out_, 3, 1);
}

void DialogSNPLDB::setup_group_box_option()
{
    QHBoxLayout *midlayout = new QHBoxLayout;

    QGroupBox *box1 = new QGroupBox("Option", this);
    midlayout->addWidget(box1);

    QFormLayout *form = new QFormLayout;
    box1->setLayout(form);

    maf_ = new QLineEdit(box1);
    maf_->setPlaceholderText("0.01");
    maf_->setToolTip("value range [0, 0.5]");
    form->addRow("Minimum haplotype frequency:", maf_);

    maxlen_ = new QSpinBox(box1);
    maxlen_->setRange(2, 99999999);
    maxlen_->setSingleStep(10000);
    form->addRow("Maximum block length:", maxlen_);

    llim_ = new QSpinBox(box1);
    llim_->setRange(0, 100);
    form->addRow("Strong LD lower limit:", llim_);

    ulim_ = new QSpinBox(box1);
    ulim_->setRange(0, 100);
    form->addRow("Strong LD upper limit:", ulim_);

    recomb_ = new QSpinBox(box1);
    recomb_->setRange(0, 100);
    form->addRow("Strong recombination threshold:", recomb_);

    inform_ = new QDoubleSpinBox(box1);
    inform_->setRange(0.0, 1.0);
    inform_->setSingleStep(0.01);
    form->addRow("Informative strong LD threshold:", inform_);

    thread_ = new QSpinBox(box1);
    form->addRow("Number of threads:", thread_);

    identity_ = new QLineEdit(box1);
    identity_->setPlaceholderText("0");
    identity_->setToolTip("value range [0, 1]");
    form->addRow("Haplotype identity:", identity_);

    QGroupBox *box2 = new QGroupBox("Population", this);
    midlayout->addWidget(box2);

    QVBoxLayout *vbox = new QVBoxLayout;
    box2->setLayout(vbox);

    natchk_ = new QCheckBox("Natural/Germplasm", box2);
    rilchk_ = new QCheckBox("Recombinant inbred line (RIL)", box2);
    namchk_ = new QCheckBox("Nested association mapping (NAM)", box2);

    QButtonGroup *btngrp = new QButtonGroup(box2);
    btngrp->addButton(natchk_);
    btngrp->addButton(rilchk_);
    btngrp->addButton(namchk_);
    btngrp->setExclusive(true);

    natchk_->setChecked(true);
    nam_ = new QLineEdit(box2);
    nam_->setPlaceholderText("family sample size, n1,n2,...");
    nam_->setDisabled(true);

    connect(namchk_, SIGNAL(stateChanged(int)), this, SLOT(set_nam_state(int)));

    vbox->addWidget(natchk_);
    vbox->addWidget(rilchk_);
    vbox->addWidget(namchk_);
    vbox->addWidget(nam_);

    QSpacerItem *spacer = new QSpacerItem(1, 1, QSizePolicy::Minimum, QSizePolicy::Expanding);
    vbox->addItem(spacer);

    static_cast<QVBoxLayout*>(layout())->addLayout(midlayout);
}

void DialogSNPLDB::setup_button_box()
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

void DialogSNPLDB::get_vcf()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        vcf_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogSNPLDB::get_gene()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        gene_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogSNPLDB::get_block()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        block_->setText(filename);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogSNPLDB::restore_defaults()
{
    vcf_->clear();
    gene_->clear();
    block_->clear();
    out_->setText("snpldb.out." + current_datetime_string());
    maf_->setText("0.01");
    maxlen_->setValue(100000);
    llim_->setValue(70);
    ulim_->setValue(98);
    recomb_->setValue(90);
    inform_->setValue(0.95);
    thread_->setValue(0);
    identity_->setText("0");
    natchk_->setChecked(true);
    nam_->clear();
}

void DialogSNPLDB::set_nam_state(int state)
{
    if (state == Qt::Checked)
        nam_->setEnabled(true);
    else if (state == Qt::Unchecked)
        nam_->setDisabled(true);
}
