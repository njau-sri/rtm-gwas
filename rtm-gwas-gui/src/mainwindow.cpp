#include <QDir>
#include <QUrl>
#include <QTextCodec>
#include <QFileInfo>
#include <QFileDialog>
#include <QMessageBox>
#include <QDesktopServices>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "parameter.h"
#include "dialogsnpldb.h"
#include "dialoggsc.h"
#include "dialogassoc.h"
#include "version.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    proc_(0), bar_(0)
{
    ui->setupUi(this);

    ui->actionExit->setMenuRole(QAction::QuitRole);
    ui->actionAbout->setMenuRole(QAction::AboutRole);

    ui->plainTextEditLog->setMaximumBlockCount(Parameter::logsize);

    bar_ = new QProgressBar;
    bar_->setTextVisible(false);
    bar_->setMaximumHeight(statusBar()->height() / 2);
    bar_->setMaximumWidth(150);
    statusBar()->addPermanentWidget(bar_);

    proc_ = new QProcess(this);
    proc_->setProcessChannelMode(QProcess::MergedChannels);
    connect(proc_, SIGNAL(readyReadStandardOutput()), this, SLOT(slot_proc_readyReadStandardOutput()));
    connect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));

    statusBar()->showMessage(tr("Ready"));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::closeEvent(QCloseEvent *e)
{
    if ( isProcessRunning() ) {
        e->ignore();
        return;
    }

    if ( Parameter::delete_onexit ) {
        foreach (const QString &name, mapfile_)
            QFile::remove(name);
    }

    e->accept();
}

void MainWindow::on_actionChangeDir_triggered()
{
    QString dir = QFileDialog::getExistingDirectory(this, tr("Change working directory"), Parameter::work);
    if ( ! dir.isEmpty() ) {
        Parameter::work = dir;
        QDir::setCurrent(dir);
    }
}

void MainWindow::on_actionSNPLDB_triggered()
{
    DialogSNPLDB d(this);
    if (d.exec() != QDialog::Accepted)
        return;
    startProcess(d.getProg(), d.getArgs());
}

void MainWindow::on_actionGSC_triggered()
{
    DialogGSC d(this);
    if (d.exec() != QDialog::Accepted)
        return;
    startProcess(d.getProg(), d.getArgs());
}

void MainWindow::on_actionAssociation_triggered()
{
    DialogAssoc d(this);
    if (d.exec() != QDialog::Accepted)
        return;
    startProcess(d.getProg(), d.getArgs());
}

void MainWindow::on_actionContents_triggered()
{
    QDir dir = QApplication::applicationDirPath();
    QString fileName = dir.filePath(QLatin1String("RTM-GWAS_UserGuide.pdf"));
    QString url = QLatin1String("https://github.com/njau-sri/rtm-gwas/wiki");

    if ( ! QFile::exists(fileName) )
        QDesktopServices::openUrl(QUrl(url));
    else if ( ! QDesktopServices::openUrl(QUrl(QLatin1String("file:///") + fileName)) ) {
        QMessageBox::critical(this, tr("ERROR"),
            tr("Can't find suitable application to open: %1<br>"
               "Please visit the wiki: <a href=%2>%2</a>").arg(fileName, url));
        QDesktopServices::openUrl(QUrl(url));
    }
}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("About RTM-GWAS"),
        tr(
    "<h3>RTM-GWAS " RTM_GWAS_VERSION " (Built on %1 at %2)</h3>"
    "<p><a href=%3>%3</a></p>"
        ).arg(QLatin1String(__DATE__), QLatin1String(__TIME__),
              QLatin1String("https://github.com/njau-sri/rtm-gwas")));
}

void MainWindow::on_listWidgetFile_currentRowChanged(int currentRow)
{
    if (currentRow >= 0 && currentRow < ui->listWidgetFile->count()) {
        QString fileName = ui->listWidgetFile->item(currentRow)->text();
        showFile(mapfile_.value(fileName));
    }
    else {
        ui->plainTextEditData->clear();
    }
}

void MainWindow::slot_proc_readyReadStandardOutput()
{
    QByteArray v = proc_->readAllStandardOutput();
    ui->plainTextEditLog->moveCursor(QTextCursor::End);
    ui->plainTextEditLog->insertPlainText(QString::fromLocal8Bit(v.data(),v.size()));
    ui->plainTextEditLog->moveCursor(QTextCursor::End);
}

void MainWindow::slot_proc_finished(int code, QProcess::ExitStatus status)
{
    statusBar()->clearMessage();
    bar_->setRange(0,1);
    bar_->reset();

    if (status != QProcess::NormalExit || code != 0) {
        foreach (const QString &name, getOutputFiles())
            QFile::remove(name);
        QMessageBox::critical(this, tr("ERROR"), tr("Process exited unexpectedly: code %1, status %2.").arg(code).arg(status));
        return;
    }

    loadFiles(getOutputFiles());
}

void MainWindow::loadFiles(const QStringList &fileNames)
{
    foreach (const QString &name, fileNames) {
        QFileInfo fi(name);
        QFile file(name);
        if (!mapfile_.contains(fi.fileName()) && file.open(QIODevice::ReadOnly)) {
            mapfile_.insert(fi.fileName(),name);
            ui->listWidgetFile->addItem(fi.fileName());
        }
    }
    ui->listWidgetFile->setCurrentRow(ui->listWidgetFile->count() - 1);
}

void MainWindow::showFile(const QString &fileName)
{
    QString name = mapfile_.value(fileName, fileName);

    QFile file(name);
    if ( ! file.open(QIODevice::ReadOnly) ) {
        QMessageBox::critical(this, tr("ERROR"), tr("Can't open file: %1").arg(name));
        return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);

    qint64 maxlen = Parameter::txtsize * 1024 * 1024;

    QByteArray v = file.read(maxlen);

    QTextCodec::ConverterState state;
    QTextCodec *codec = QTextCodec::codecForName("UTF-8");
    QString text = codec->toUnicode(v.data(), v.size(), &state);
    if (state.invalidChars > 0)
        text = QString::fromLocal8Bit(v);

    ui->plainTextEditData->setPlainText(text);

    if (v.size() < file.size())
        ui->plainTextEditData->appendHtml(tr("<br>&lt;<i>Not entire contents of file being displayed</i>&gt;"));

    QApplication::restoreOverrideCursor();
}

bool MainWindow::isProcessRunning()
{
    if (proc_->state() != QProcess::NotRunning) {
        if (QMessageBox::question(this, tr("Terminate"), tr("Stop currently running computations?"), QMessageBox::Ok | QMessageBox::Cancel) != QMessageBox::Ok)
            return true;

        disconnect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
        proc_->close();

        statusBar()->clearMessage();
        bar_->setRange(0,1);
        bar_->reset();

        foreach (const QString &name, getOutputFiles())
            QFile::remove(name);

        connect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
    }

    return false;
}

QStringList MainWindow::getOutputFiles() const
{
    QStringList output;

    int pos = args_.indexOf(QLatin1String("--out"));
    if (pos == -1)
        return output;

    QFileInfo fi(args_.at(pos+1));
    QStringList filter(fi.fileName().append(QLatin1Char('*')));
    QDir dir = fi.absoluteDir();
    QStringList fileNames = dir.entryList(filter, QDir::Files | QDir::NoSymLinks | QDir::CaseSensitive);

    foreach (const QString &name, fileNames)
        output.append(dir.filePath(name));

    return output;
}

void MainWindow::startProcess(const QString &prog, const QStringList &args)
{
    if ( isProcessRunning() )
        return;

    args_ = args;
    proc_->setWorkingDirectory(Parameter::work);

    if (Parameter::openmp > 0) {
        QStringList env = QProcess::systemEnvironment();
        env << QString(QLatin1String("OMP_NUM_THREADS=%1")).arg(Parameter::openmp);
        proc_->setEnvironment(env);
    }

    proc_->start(prog, args);
    if ( ! proc_->waitForStarted() ) {
        QMessageBox::critical(this, tr("ERROR"), tr("Can't start process: %1").arg(prog));
        return;
    }

    statusBar()->showMessage(tr("Busy"));
    bar_->setRange(0,0);
    bar_->reset();
}
