#include "dialoggconv.h"

#include <QLabel>
#include <QCheckBox>
#include <QComboBox>
#include <QDateTime>
#include <QGroupBox>
#include <QLineEdit>
#include <QFileDialog>
#include <QGridLayout>
#include <QPushButton>
#include <QSpacerItem>
#include <QVBoxLayout>
#include <QButtonGroup>
#include <QDialogButtonBox>

#include "parameter.h"

namespace {

QString current_datetime_string()
{
    return QDateTime::currentDateTime().toString("yyMMddhhmmsszzz");
}

} // namespace

DialogGConv::DialogGConv(QWidget *parent) : QDialog(parent)
{
    setWindowTitle("GCONV - RTM-GWAS");

    setup_layout();
    setup_group_box();
    setup_button_box();

    restore_defaults();
}

QString DialogGConv::program() const
{
    return QDir(par->bin_path).absoluteFilePath("rtm-gwas-gconv");
}

QStringList DialogGConv::arguments() const
{
    QStringList args;

    QStringList format, suffix;
    format << "--vcf" << "--ped" << "--hmp" << "--geno";
    suffix << ".vcf" << ".ped" << ".hmp" << ".geno";

    if (!in_file_->text().isEmpty())
        args << format.at(in_format_->currentIndex()) << in_file_->text();

    if (!out_file_->text().isEmpty()) {
        int wh = out_format_->currentIndex();
        QString prefix = out_file_->text();
        if (!prefix.endsWith(suffix.at(wh)))
            prefix.append(suffix.at(wh));
        args << "--out" << QDir(par->working_directory).absoluteFilePath(prefix);
    }

    if (sort_->isChecked())
        args << "--sort";

    return args;
}

void DialogGConv::setup_layout()
{
    setLayout(new QVBoxLayout);
}

void DialogGConv::setup_group_box()
{
    QGroupBox *box = new QGroupBox("Data", this);
    layout()->addWidget(box);

    QGridLayout *grid = new QGridLayout;
    box->setLayout(grid);

    QLabel *inl = new QLabel("Input:", box);
    in_format_ = new QComboBox(box);
    in_format_->addItems(QStringList() << "VCF" << "PED/MAP" << "HapMap" << "Geno");
    in_file_ = new QLineEdit(box);
    QPushButton *inb = new QPushButton("Browse...", box);
    grid->addWidget(inl, 0, 0);
    grid->addWidget(in_format_, 0, 1);
    grid->addWidget(in_file_, 0, 2);
    grid->addWidget(inb, 0, 3);
    connect(inb, SIGNAL(clicked()), this, SLOT(get_input()));

    QLabel *outl = new QLabel("Output:", box);
    out_format_ = new QComboBox(box);
    out_format_->addItems(QStringList() << "VCF" << "PED/MAP" << "HapMap" << "Geno");
    out_file_ = new QLineEdit(box);
    grid->addWidget(outl, 1, 0);
    grid->addWidget(out_format_, 1, 1);
    grid->addWidget(out_file_, 1, 2);

    sort_ = new QCheckBox("Sort chromosome position", box);
    grid->addWidget(sort_, 2, 2);
}

void DialogGConv::setup_button_box()
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

void DialogGConv::get_input()
{
    QString filename = QFileDialog::getOpenFileName(this, QString(), par->file_dialog_directory);
    if (!filename.isEmpty()) {
        in_file_->setText(filename);
        if (filename.endsWith(".vcf"))
            in_format_->setCurrentIndex(0);
        else if (filename.endsWith(".ped") || filename.endsWith(".map"))
            in_format_->setCurrentIndex(1);
        else if (filename.endsWith(".hmp"))
            in_format_->setCurrentIndex(2);
        else if (filename.endsWith(".geno"))
            in_format_->setCurrentIndex(3);
        par->file_dialog_directory = QFileInfo(filename).absolutePath();
    }
}

void DialogGConv::restore_defaults()
{
    in_file_->clear();
    in_format_->setCurrentIndex(0);
    out_file_->setText("gconv.out." + current_datetime_string());
    out_format_->setCurrentIndex(0);
    sort_->setChecked(false);
}
