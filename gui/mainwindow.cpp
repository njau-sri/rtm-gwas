#include <QDir>
#include <QUrl>
#include <QFile>
#include <QTextCodec>
#include <QFileDialog>
#include <QMessageBox>
#include <QStatusBar>
#include <QDesktopServices>

#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "params.h"

#include "dialogldb.h"
#include "dialogeigen.h"
#include "dialogassoc.h"
#include "dialogdata.h"
#include "dialogsummary.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    process(0),
    progressBar(0)
{
    ui->setupUi(this);

    ui->actionExit->setMenuRole(QAction::QuitRole);
    ui->actionAbout->setMenuRole(QAction::AboutRole);

    ui->plainTextEditLogViewer->setMaximumBlockCount(Params::log_size);
    ui->plainTextEditLogViewer->setStyleSheet(
        QLatin1String("QPlainTextEdit{background-color:rgb(0,0,0);color:rgb(0,255,0);}"));

    progressBar = new QProgressBar;
    progressBar->setTextVisible(false);
    progressBar->setMaximumHeight(statusBar()->height()/2);
    progressBar->setMaximumWidth(150);
    statusBar()->addPermanentWidget(progressBar);

    process = new QProcess(this);
    process->setProcessChannelMode(QProcess::MergedChannels);
    connect(process, SIGNAL(readyReadStandardOutput()), this, SLOT(slot_process_readyReadStandardOutput()));
    connect(process, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_process_finished(int,QProcess::ExitStatus)));

    setWindowTitle(QLatin1String("RTM-GWAS"));
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

    if ( false ) {
        foreach (const QString &name, fileMap)
            QFile::remove(name);
    }

    e->accept();
}

void MainWindow::on_actionChangeDir_triggered()
{
    QString dir = QFileDialog::getExistingDirectory(this, tr("Change Working Directory"), Params::work_dir);
    if ( ! dir.isEmpty() ) {
        Params::work_dir = dir;
        QDir::setCurrent(dir);
    }
}

void MainWindow::on_actionSNPLDB_triggered()
{
    DialogLDB dlg(this);
    if (dlg.exec() != QDialog::Accepted)
        return;
    startProcess(Params::exe, dlg.arguments());
}

void MainWindow::on_actionEigenGSC_triggered()
{
    DialogEigen dlg(this);
    if (dlg.exec() != QDialog::Accepted)
        return;
    startProcess(Params::exe, dlg.arguments());
}

void MainWindow::on_actionAssociation_triggered()
{
    DialogAssoc dlg(this);
    if (dlg.exec() != QDialog::Accepted)
        return;
    startProcess(Params::exe, dlg.arguments());
}

void MainWindow::on_actionGenotype_triggered()
{
    DialogData dlg(this);
    if (dlg.exec() != QDialog::Accepted)
        return;
    startProcess(Params::exe, dlg.arguments());
}

void MainWindow::on_actionSummary_triggered()
{
    DialogSummary dlg(this);
    if (dlg.exec() != QDialog::Accepted)
        return;
    startProcess(Params::exe, dlg.arguments());
}

void MainWindow::on_actionContents_triggered()
{
    QDir dir = QApplication::applicationDirPath();
    QString fileName = dir.filePath(QLatin1String("RTM-GWAS_UserGuide.pdf"));
    QString urlWiki = QLatin1String("https://github.com/njau-sri/rtm-gwas/wiki");

    if ( ! QFile::exists(fileName) )
        QDesktopServices::openUrl(QUrl(urlWiki));
    else if ( ! QDesktopServices::openUrl(QUrl(QLatin1String("file:///") + fileName)) ) {
        QMessageBox::critical(this, tr("ERROR"),
            tr("Can't find suitable application to open: %1<br>"
               "Please visit the wiki: <a href=%2>%2</a>").arg(fileName, urlWiki));
        QDesktopServices::openUrl(QUrl(urlWiki));
    }
}

void MainWindow::on_actionAbout_triggered()
{
    QMessageBox::about(this, tr("About RTM-GWAS"),
        tr(
    "<h3>RTM-GWAS v1.0 (Built on %1 at %2)</h3>"
    "<p><a href=%3>%3</a></p>"
        ).arg(QLatin1String(__DATE__), QLatin1String(__TIME__),
              QLatin1String("https://github.com/njau-sri/rtm-gwas")));
}

void MainWindow::on_listWidgetDataList_currentRowChanged(int currentRow)
{
    if (currentRow >= 0 && currentRow < ui->listWidgetDataList->count()) {
        QString fileName = ui->listWidgetDataList->item(currentRow)->text();
        showFile(fileMap.value(fileName));
    }
    else {
        ui->plainTextEditDataViewer->clear();
    }
}

void MainWindow::slot_process_readyReadStandardOutput()
{
    QByteArray ba = process->readAllStandardOutput();
    ui->plainTextEditLogViewer->moveCursor(QTextCursor::End);
    ui->plainTextEditLogViewer->insertPlainText(QString::fromLocal8Bit(ba.data(),ba.size()));
    ui->plainTextEditLogViewer->moveCursor(QTextCursor::End);
}

void MainWindow::slot_process_finished(int code, QProcess::ExitStatus status)
{
    statusBar()->clearMessage();
    progressBar->setRange(0,1);
    progressBar->reset();

    if (status != QProcess::NormalExit || code != 0) {
        foreach (const QString &name, getOutputFiles())
            QFile::remove(name);
        QMessageBox::critical(this, tr("ERROR"),
            tr("The process exited unexpectedly: code %1, status %2.").arg(code).arg(status));
        return;
    }

    loadFiles(getOutputFiles());
}

void MainWindow::loadFiles(const QStringList &fileNames)
{
    foreach (const QString &name, fileNames) {
        QFileInfo fi(name);
        QFile file(name);
        if (!fileMap.contains(fi.fileName()) && file.open(QIODevice::ReadOnly)) {
            fileMap.insert(fi.fileName(),name);
            ui->listWidgetDataList->addItem(fi.fileName());
        }
    }
    ui->listWidgetDataList->setCurrentRow(ui->listWidgetDataList->count() - 1);
}

void MainWindow::showFile(const QString &fileName)
{
    QString name = fileMap.value(fileName, fileName);

    QFile file(name);
    if ( ! file.open(QIODevice::ReadOnly) ) {
        QMessageBox::critical(this, tr("ERROR"), tr("Can't open file: %1").arg(name));
        return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);

    qint64 maxlen = Params::txt_size * 1024 * 1024;

    QByteArray ba = file.read(maxlen);

    QTextCodec::ConverterState state;
    QTextCodec *codec = QTextCodec::codecForName("UTF-8");
    QString text = codec->toUnicode(ba.data(), ba.size(), &state);
    if (state.invalidChars > 0)
        text = QString::fromLocal8Bit(ba);

    ui->plainTextEditDataViewer->setPlainText(text);

    if (ba.size() < file.size())
        ui->plainTextEditDataViewer->appendHtml(tr("<br>&lt;<i>Not entire contents of file being displayed</i>&gt;"));

    QApplication::restoreOverrideCursor();
}

bool MainWindow::isProcessRunning()
{
    if (process->state() != QProcess::NotRunning) {
        QMessageBox::StandardButton ret = QMessageBox::question(this,
            tr("Terminate"), tr("Stop currently running process?"),
            QMessageBox::Ok | QMessageBox::Cancel);

        if (ret != QMessageBox::Ok)
            return true;

        disconnect(process, SIGNAL(finished(int,QProcess::ExitStatus)),
                this, SLOT(slot_process_finished(int,QProcess::ExitStatus)));

        process->close();

        statusBar()->clearMessage();
        progressBar->setRange(0,1);
        progressBar->reset();

        foreach (const QString &name, getOutputFiles())
            QFile::remove(name);

        connect(process, SIGNAL(finished(int,QProcess::ExitStatus)),
                this, SLOT(slot_process_finished(int,QProcess::ExitStatus)));
    }

    return false;
}

QStringList MainWindow::getOutputFiles() const
{
    QStringList output;

    int pos = argv.indexOf(QLatin1String("--out"));
    if (pos == -1)
        return output;

    QFileInfo fi(argv.at(pos+1));
    QStringList filter(fi.fileName().append(QLatin1Char('*')));
    QDir dir = fi.absoluteDir();
    QStringList fileNames = dir.entryList(filter, QDir::Files | QDir::NoSymLinks | QDir::CaseSensitive);

    foreach (const QString &name, fileNames)
        output.append(dir.filePath(name));

    return output;
}

void MainWindow::startProcess(const QString &program, const QStringList &arguments)
{
    if ( isProcessRunning() )
        return;

    argv = arguments;

    ui->plainTextEditLogViewer->moveCursor(QTextCursor::End);
    ui->plainTextEditLogViewer->insertPlainText(
        QLatin1String("CMD: ") + arguments.join(QLatin1String(" ")) + QLatin1Char('\n'));
    ui->plainTextEditLogViewer->moveCursor(QTextCursor::End);

    process->setWorkingDirectory(Params::work_dir);

    process->start(program, arguments);
    if ( ! process->waitForStarted() ) {
        QMessageBox::critical(this, tr("ERROR"), tr("Can't start process: %1").arg(program));
        return;
    }

    statusBar()->showMessage(tr("Busy"));
    progressBar->setRange(0,0);
    progressBar->reset();
}
