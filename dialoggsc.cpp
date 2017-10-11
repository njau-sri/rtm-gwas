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

    ui->lineEditVcf->setText(Parameter::geno);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));
}

DialogGSC::~DialogGSC()
{
    delete ui;
}

QString DialogGSC::getProg() const
{
    return QDir(QApplication::applicationDirPath()).filePath(QLatin1String("gsc"));
}

QStringList DialogGSC::getArgs() const
{
    QStringList args;

    if ( ! ui->lineEditVcf->text().isEmpty() )
        args << QLatin1String("--vcf") << ui->lineEditVcf->text();

    if ( ! ui->lineEditTop->text().isEmpty() )
        args << QLatin1String("--top") << ui->lineEditTop->text();

    QString prefix = QLatin1String("gsc.out.");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    args << QLatin1String("--out") << QDir(Parameter::work).absoluteFilePath(prefix);

    return args;
}

void DialogGSC::apply()
{
    Parameter::geno = ui->lineEditVcf->text();
}

void DialogGSC::on_pushButtonVcf_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose VCF file"), Parameter::open);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditVcf->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}
