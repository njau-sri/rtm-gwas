#include <QDir>
#include <QDateTime>
#include <QFileInfo>
#include <QFileDialog>

#include "dialogeigen.h"
#include "ui_dialogeigen.h"

#include "params.h"

DialogEigen::DialogEigen(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogEigen)
{
    ui->setupUi(this);

    ui->comboBoxGenotype->setCurrentIndex(Params::gen_type);
    ui->lineEditGenotype->setText(Params::geno);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));

    setWindowTitle(tr("Eigen GSC"));
}

DialogEigen::~DialogEigen()
{
    delete ui;
}

QStringList DialogEigen::arguments() const
{
    QStringList argv;

    argv << QLatin1String("gsc");

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

    if ( ! ui->lineEditTop->text().isEmpty() )
        argv << QLatin1String("--top") << ui->lineEditTop->text();

    QString prefix = QLatin1String("appgsc_");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    argv << QLatin1String("--out") << QDir(Params::work_dir).absoluteFilePath(prefix);

    return argv;
}

void DialogEigen::apply()
{
    Params::gen_type = ui->comboBoxGenotype->currentIndex();
    Params::geno = ui->lineEditGenotype->text();
}

void DialogEigen::on_pushButtonGenotype_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Genotype File"), Params::open_dir);
    if ( ! fileName.isEmpty() ) {
        ui->lineEditGenotype->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}
