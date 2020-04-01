#include <QDir>
#include <QDateTime>
#include <QFileDialog>
#include "dialoggsc.h"
#include "ui_dialoggsc.h"
#include "parameter.h"

DialogGSC::DialogGSC(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogGSC)
{
    ui->setupUi(this);

    ui->lineEditVCF->setText(Parameter::vcf);
    ui->lineEditGRM->setText(Parameter::grm);
    ui->spinBoxOpenMP->setValue(Parameter::openmp);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));
}

DialogGSC::~DialogGSC()
{
    delete ui;
}

QString DialogGSC::getProg() const
{
    return QDir(Parameter::exe).filePath(QLatin1String("rtm-gwas-gsc"));
}

QStringList DialogGSC::getArgs() const
{
    QStringList args;

    if ( ! ui->lineEditVCF->text().isEmpty() )
        args << QLatin1String("--vcf") << ui->lineEditVCF->text();

    if ( ! ui->lineEditGRM->text().isEmpty() )
        args << QLatin1String("--grm") << ui->lineEditGRM->text();

    if ( ! ui->lineEditTop->text().isEmpty() )
        args << QLatin1String("--top") << ui->lineEditTop->text();

    if (ui->spinBoxOpenMP->value() > 0)
        args << QLatin1String("--thread") << ui->spinBoxOpenMP->text();

    QString prefix = QLatin1String("gsc.out.");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    args << QLatin1String("--out") << QDir(Parameter::work).absoluteFilePath(prefix);

    return args;
}

void DialogGSC::apply()
{
    Parameter::vcf = ui->lineEditVCF->text();
    Parameter::grm = ui->lineEditGRM->text();
    Parameter::openmp = ui->spinBoxOpenMP->value();
}

void DialogGSC::on_pushButtonVCF_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose VCF file"), Parameter::open);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditVCF->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}

void DialogGSC::on_pushButtonGRM_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose GRM file"), Parameter::open);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditGRM->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}
