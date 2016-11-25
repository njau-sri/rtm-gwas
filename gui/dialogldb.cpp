#include <QDir>
#include <QDateTime>
#include <QFileInfo>
#include <QFileDialog>

#include "dialogldb.h"
#include "ui_dialogldb.h"

#include "params.h"

DialogLDB::DialogLDB(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DialogLDB)
{
    ui->setupUi(this);

    ui->comboBoxGenotype->setCurrentIndex(Params::gen_type);
    ui->lineEditGenotype->setText(Params::geno);
    ui->lineEditBlock->setText(Params::block);

    connect(ui->buttonBox->button(QDialogButtonBox::Ok), SIGNAL(clicked()), this, SLOT(apply()));

    setWindowTitle(tr("SNPLDB"));
}

DialogLDB::~DialogLDB()
{
    delete ui;
}

QStringList DialogLDB::arguments() const
{
    QStringList argv;

    argv << QLatin1String("ldb");

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

    if ( ! ui->lineEditBlock->text().isEmpty() )
        argv << QLatin1String("--block") << ui->lineEditBlock->text();

    if ( ! ui->lineEditMAF->text().isEmpty() )
        argv << QLatin1String("--mhf") << ui->lineEditMAF->text();

    if ( ! ui->lineEditMaxDist->text().isEmpty() )
        argv << QLatin1String("--maxlen") << ui->lineEditMaxDist->text();

    if ( ! ui->lineEditLowCI->text().isEmpty() )
        argv << QLatin1String("--low") << ui->lineEditLowCI->text();

    if ( ! ui->lineEditHighCI->text().isEmpty() )
        argv << QLatin1String("--high") << ui->lineEditHighCI->text();

    if ( ! ui->lineEditRecHighCI->text().isEmpty() )
        argv << QLatin1String("--rec") << ui->lineEditRecHighCI->text();

    if ( ! ui->lineEditInfoFrac->text().isEmpty() )
        argv << QLatin1String("--frac") << ui->lineEditInfoFrac->text();

    QString prefix = QLatin1String("appldb_");
    prefix += QDateTime::currentDateTime().toString(QLatin1String("yyMMdd_hhmmsszzz"));

    argv << QLatin1String("--out") << QDir(Params::work_dir).absoluteFilePath(prefix);

    return argv;
}

void DialogLDB::apply()
{
    Params::gen_type = ui->comboBoxGenotype->currentIndex();
    Params::geno = ui->lineEditGenotype->text();
    Params::block = ui->lineEditBlock->text();
}

void DialogLDB::on_pushButtonGenotype_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Genotype File"), Params::open_dir);
    if (!fileName.isEmpty()) {
        ui->lineEditGenotype->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}

void DialogLDB::on_pushButtonBlock_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Block File"), Params::open_dir);
    if (!fileName.isEmpty()) {
        ui->lineEditBlock->setText(fileName);
        Params::open_dir = QFileInfo(fileName).absolutePath();
    }
}
