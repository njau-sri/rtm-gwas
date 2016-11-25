#include <QDir>
#include <QDateTime>
#include <QFileInfo>
#include <QFileDialog>

#include "dialogdata.h"
#include "ui_dialogdata.h"

#include "params.h"

DialogData::DialogData(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogData)
{
    ui->setupUi(this);

    ui->comboBoxGenotype->setCurrentIndex(Params::gen_type);
    ui->lineEditGenotype->setText(Params::geno);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));

    setWindowTitle(tr("Genotype Data"));
}

DialogData::~DialogData()
{
    delete ui;
}

QStringList DialogData::arguments() const
{
    QStringList argv;

    argv << QLatin1String("data");

    if ( ! ui->lineEditGenotype->text().isEmpty() ) {
        int indx = ui->comboBoxGenotype->currentIndex();
        if (indx == 1)
            argv << QLatin1String("--ped");
        else if (indx == 2)
            argv << QLatin1String("--hmp");
        else if (indx == 3)
            argv << QLatin1String("--geno");
        else
            argv << QLatin1String("--vcf");
        argv << ui->lineEditGenotype->text();
        argv << QLatin1String("--format");
        indx = ui->comboBoxOutputFormat->currentIndex();
        if (indx == 1)
            argv << QLatin1String("ped");
        else if (indx == 2)
            argv << QLatin1String("hmp");
        else if (indx == 3)
            argv << QLatin1String("geno");
        else
            argv << QLatin1String("vcf");
    }

    if ( ! ui->lineEditLocus->text().isEmpty() )
        argv << QLatin1String("--loc") << ui->lineEditLocus->text();

    if ( ! ui->lineEditIndividual->text().isEmpty() )
        argv << QLatin1String("--ind") << ui->lineEditIndividual->text();

    if ( ! ui->lineEditChromosome->text().isEmpty() )
        argv << QLatin1String("--chr") << ui->lineEditChromosome->text();

    if ( ! ui->lineEditMinMAF->text().isEmpty() )
        argv << QLatin1String("--maf") << ui->lineEditMinMAF->text();

    if ( ! ui->lineEditMaxLocusMissing->text().isEmpty() )
        argv << QLatin1String("--mloc") << ui->lineEditMaxLocusMissing->text();

    if ( ! ui->lineEditMaxIndivMissing->text().isEmpty() )
        argv << QLatin1String("--mind") << ui->lineEditMaxIndivMissing->text();

    if ( ! ui->lineEditMaxHeteroz->text().isEmpty() )
        argv << QLatin1String("--het") << ui->lineEditMaxHeteroz->text();

    if ( ! ui->lineEditMinAllele->text().isEmpty() )
        argv << QLatin1String("--min-alleles") << ui->lineEditMinAllele->text();

    if ( ! ui->lineEditMaxAllele->text().isEmpty() )
        argv << QLatin1String("--max-alleles") << ui->lineEditMaxAllele->text();

    if ( ui->checkBoxDiploid->isChecked() )
        argv << QLatin1String("--diploid");

    QString prefix = QLatin1String("appdata_");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    argv << QLatin1String("--out") << QDir(Params::work_dir).absoluteFilePath(prefix);

    return argv;
}

void DialogData::apply()
{
    Params::gen_type = ui->comboBoxGenotype->currentIndex();
    Params::geno = ui->lineEditGenotype->text();
}

void DialogData::on_pushButtonGenotype_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Genotype File"), Params::open_dir);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditGenotype->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}

void DialogData::on_pushButtonLocus_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose LocusList File"), Params::open_dir);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditLocus->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}

void DialogData::on_pushButtonIndividual_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose IndivList File"), Params::open_dir);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditIndividual->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}
