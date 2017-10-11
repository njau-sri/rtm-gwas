#include <QDir>
#include <QDateTime>
#include <QFileDialog>
#include "dialogassoc.h"
#include "ui_dialogassoc.h"
#include "parameter.h"

DialogAssoc::DialogAssoc(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogAssoc)
{
    ui->setupUi(this);

    ui->lineEditVcf->setText(Parameter::geno);
    ui->lineEditPheno->setText(Parameter::pheno);
    ui->lineEditCovar->setText(Parameter::covar);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));
}

DialogAssoc::~DialogAssoc()
{
    delete ui;
}

QString DialogAssoc::getProg() const
{
    return QDir(QApplication::applicationDirPath()).filePath(QLatin1String("assoc"));
}

QStringList DialogAssoc::getArgs() const
{
    QStringList args;

    if ( ! ui->lineEditVcf->text().isEmpty() )
        args << QLatin1String("--vcf") << ui->lineEditVcf->text();

    if ( ! ui->lineEditPheno->text().isEmpty() )
        args << QLatin1String("--pheno") << ui->lineEditPheno->text();

    if ( ! ui->lineEditCovar->text().isEmpty() )
        args << QLatin1String("--covar") << ui->lineEditCovar->text();

    if ( ! ui->lineEditAlpha->text().isEmpty() )
        args << QLatin1String("--alpha") << ui->lineEditAlpha->text();

    if ( ! ui->lineEditPreselect->text().isEmpty() )
        args << QLatin1String("--preselect") << ui->lineEditPreselect->text();

    if ( ! ui->lineEditRsq->text().isEmpty() )
        args << QLatin1String("--rsq") << ui->lineEditRsq->text();

    if (ui->comboBoxMtc->currentIndex() != 0)
        args << QLatin1String("--mtc") << ui->comboBoxMtc->currentText();

    args << QLatin1String("--sstype") << ui->comboBoxStype->currentText();

    if ( ! ui->checkBoxGxe->isChecked() )
        args << QLatin1String("--no-gxe");

    QString prefix = QLatin1String("assoc.out.");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    args << QLatin1String("--out") << QDir(Parameter::work).absoluteFilePath(prefix);

    return args;
}

void DialogAssoc::apply()
{
    Parameter::geno = ui->lineEditVcf->text();
    Parameter::pheno = ui->lineEditPheno->text();
    Parameter::covar = ui->lineEditCovar->text();
}

void DialogAssoc::on_pushButtonVcf_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose VCF file"), Parameter::open);
    if (!fileName.isEmpty()) {
        ui->lineEditVcf->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}

void DialogAssoc::on_pushButtonPheno_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose phenotype file"), Parameter::open);
    if (!fileName.isEmpty()) {
        ui->lineEditPheno->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}

void DialogAssoc::on_pushButtonCovar_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose covariate file"), Parameter::open);
    if (!fileName.isEmpty()) {
        ui->lineEditCovar->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}
