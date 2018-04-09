#include <QDir>
#include <QDateTime>
#include <QFileDialog>
#include "dialogsnpldb.h"
#include "ui_dialogsnpldb.h"
#include "parameter.h"

DialogSNPLDB::DialogSNPLDB(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogSNPLDB)
{
    ui->setupUi(this);

    ui->lineEditVcf->setText(Parameter::geno);
    ui->lineEditBlock->setText(Parameter::block);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));
}

DialogSNPLDB::~DialogSNPLDB()
{
    delete ui;
}

QString DialogSNPLDB::getProg() const
{
    return QDir(QApplication::applicationDirPath()).filePath(QLatin1String("snpldb"));
}

QStringList DialogSNPLDB::getArgs() const
{
    QStringList args;

    if ( ! ui->lineEditVcf->text().isEmpty() )
        args << QLatin1String("--vcf") << ui->lineEditVcf->text();

    if ( ! ui->lineEditBlock->text().isEmpty() )
        args << QLatin1String("--blk") << ui->lineEditBlock->text();

    if ( ! ui->lineEditMaf->text().isEmpty() )
        args << QLatin1String("--maf") << ui->lineEditMaf->text();

    if ( ! ui->lineEditMaxlen->text().isEmpty() )
        args << QLatin1String("--maxlen") << ui->lineEditMaxlen->text();

    if ( ! ui->lineEditLlim->text().isEmpty() )
        args << QLatin1String("--llim") << ui->lineEditLlim->text();

    if ( ! ui->lineEditUlim->text().isEmpty() )
        args << QLatin1String("--ulim") << ui->lineEditUlim->text();

    if ( ! ui->lineEditRecomb->text().isEmpty() )
        args << QLatin1String("--recomb") << ui->lineEditRecomb->text();

    if ( ! ui->lineEditInform->text().isEmpty() )
        args << QLatin1String("--inform") << ui->lineEditInform->text();

    QString prefix = QLatin1String("snpldb.out.");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    args << QLatin1String("--out") << QDir(Parameter::work).absoluteFilePath(prefix);

    return args;
}

void DialogSNPLDB::apply()
{
    Parameter::geno = ui->lineEditVcf->text();
    Parameter::block = ui->lineEditBlock->text();
}

void DialogSNPLDB::on_pushButtonVcf_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose VCF file"), Parameter::open);
    if (!fileName.isEmpty()) {
        ui->lineEditVcf->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}

void DialogSNPLDB::on_pushButtonBlock_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose block file"), Parameter::open);
    if (!fileName.isEmpty()) {
        ui->lineEditBlock->setText(fileName);
        Parameter::open = QFileInfo(fileName).absolutePath();
    }
}
