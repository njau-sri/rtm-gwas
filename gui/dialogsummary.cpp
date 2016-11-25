#include <QDateTime>
#include <QFileDialog>
#include "dialogsummary.h"
#include "ui_dialogsummary.h"
#include "params.h"

DialogSummary::DialogSummary(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogSummary)
{
    ui->setupUi(this);

    ui->comboBoxGenotype->setCurrentIndex(Params::gen_type);
    ui->lineEditGenotype->setText(Params::geno);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));

    setWindowTitle(tr("Genotype Summary"));
}

DialogSummary::~DialogSummary()
{
    delete ui;
}


QStringList DialogSummary::arguments() const
{
    QStringList argv;

    argv << QLatin1String("sum");

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

    if ( ui->radioButtonAllele->isChecked() )
        argv << QLatin1String("--allele");
    else if ( ui->radioButtonIndiv->isChecked() )
        argv << QLatin1String("--indiv");

    QString prefix = QLatin1String("appsum_");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    argv << QLatin1String("--out") << QDir(Params::work_dir).absoluteFilePath(prefix);

    return argv;
}

void DialogSummary::apply()
{
    Params::gen_type = ui->comboBoxGenotype->currentIndex();
    Params::geno = ui->lineEditGenotype->text();
}

void DialogSummary::on_pushButtonGenotype_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Genotype File"), Params::open_dir);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditGenotype->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}
