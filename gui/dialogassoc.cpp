#include <QDir>
#include <QDateTime>
#include <QFileInfo>
#include <QFileDialog>

#include "dialogassoc.h"
#include "ui_dialogassoc.h"

#include "params.h"


DialogAssoc::DialogAssoc(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogAssoc)
{
    ui->setupUi(this);

    ui->comboBoxGenotype->setCurrentIndex(Params::gen_type);
    ui->lineEditGenotype->setText(Params::geno);
    ui->lineEditPhenotype->setText(Params::pheno);
    ui->lineEditCovariate->setText(Params::covar);
    ui->lineEditKinship->setText(Params::kin);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));

    setWindowTitle(tr("Association Analysis"));
}

DialogAssoc::~DialogAssoc()
{
    delete ui;
}

QStringList DialogAssoc::arguments() const
{
    QStringList argv;

    argv << QLatin1String("assoc");

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
    }

    if ( ! ui->lineEditPhenotype->text().isEmpty() )
        argv << QLatin1String("--pheno") << ui->lineEditPhenotype->text();

    if ( ! ui->lineEditCovariate->text().isEmpty() )
        argv << QLatin1String("--covar") << ui->lineEditCovariate->text();

    if (ui->comboBoxMethod->currentIndex() == 2 && ! ui->lineEditKinship->text().isEmpty())
        argv << QLatin1String("--kin") << ui->lineEditKinship->text();

    argv << QLatin1String("--method") << ui->comboBoxMethod->currentText();

    if ( ! ui->checkBoxGxE->isChecked() )
        argv << QLatin1String("--no-gxe");

    if (ui->comboBoxMethod->currentIndex() == 0) {
        if ( ! ui->lineEditSignificanceLevel->text().isEmpty() )
            argv << QLatin1String("--alpha") << ui->lineEditSignificanceLevel->text();

        if ( ! ui->lineEditPreselection->text().isEmpty() )
            argv << QLatin1String("--preselect") << ui->lineEditPreselection->text();

        if ( ! ui->lineEditMaxNumQTLs->text().isEmpty() )
            argv << QLatin1String("--maxqtl") << ui->lineEditMaxNumQTLs->text();

        if ( ! ui->lineEditMaxModelRSquare->text().isEmpty() )
            argv << QLatin1String("--maxrsq") << ui->lineEditMaxModelRSquare->text();

        if (ui->comboBoxMultTest->currentIndex() != 0)
            argv << QLatin1String("--multtest") << ui->comboBoxMultTest->currentText();

        argv << QLatin1String("--sstype");
        argv << ( ui->comboBoxSumSqrType->currentIndex() == 0 ? QLatin1String("1") : QLatin1String("3") );
    }

    QString prefix = QLatin1String("appassoc_");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    QDir dir(Params::work_dir);
    argv << QLatin1String("--out") << dir.absoluteFilePath(prefix);

    return argv;
}

void DialogAssoc::apply()
{
    Params::gen_type = ui->comboBoxGenotype->currentIndex();
    Params::geno = ui->lineEditGenotype->text();
    Params::pheno = ui->lineEditPhenotype->text();
    Params::covar = ui->lineEditCovariate->text();
    Params::kin = ui->lineEditKinship->text();
}

void DialogAssoc::on_pushButtonGenotype_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Genotype File"), Params::open_dir);
    if (!fileName.isEmpty()) {
        ui->lineEditGenotype->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}

void DialogAssoc::on_pushButtonPhenotype_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Phenotype File"), Params::open_dir);
    if (!fileName.isEmpty()) {
        ui->lineEditPhenotype->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}

void DialogAssoc::on_pushButtonCovariate_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Covariate File"), Params::open_dir);
    if (!fileName.isEmpty()) {
        ui->lineEditCovariate->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}

void DialogAssoc::on_pushButtonKinship_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Kinship File"), Params::open_dir);
    if (!fileName.isEmpty()) {
        ui->lineEditKinship->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}
