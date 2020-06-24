#include "mainwindow.h"

#include <QDir>
#include <QUrl>
#include <QMenu>
#include <QLabel>
#include <QMenuBar>
#include <QStatusBar>
#include <QDockWidget>
#include <QFileDialog>
#include <QListWidget>
#include <QMessageBox>
#include <QTextStream>
#include <QApplication>
#include <QProgressBar>
#include <QPlainTextEdit>
#include <QDesktopServices>

#include "parameter.h"
#include "dialogsnpldb.h"
#include "dialoggsc.h"
#include "dialogassoc.h"

#ifndef RTM_GWAS_VERSION
#define RTM_GWAS_VERSION "unknown"
#endif // RTM_GWAS_VERSION

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
    setWindowTitle("RTM-GWAS");

    setup_menu_file();
    setup_menu_analysis();
    setup_menu_help();

    setup_central_widget();
    steup_dock_widget();
    setup_status_bar();

    create_process();

    status_->setText("Ready");
}

void MainWindow::closeEvent(QCloseEvent *e)
{
    if (is_process_running()) {
        e->ignore();
        return;
    }
    e->accept();
}

void MainWindow::setup_menu_file()
{
    QMenu *file = menuBar()->addMenu("&File");
    file->addAction("&Set Working Directory", this, SLOT(set_working_directory()));
    file->addAction("Show &Working Directory", this, SLOT(show_working_directory()));
    file->addSeparator();
    QAction *act = file->addAction("E&xit", this, SLOT(close()), QKeySequence("Ctrl+Q"));
    act->setMenuRole(QAction::QuitRole); // macOS
}

void MainWindow::setup_menu_analysis()
{
    QMenu *analysis = menuBar()->addMenu("A&nalysis");
    analysis->addAction("&SNPLDB", this, SLOT(show_dialog_snpldb()));
    analysis->addAction("&GSC", this, SLOT(show_dialog_gsc()));
    analysis->addAction("&Association", this, SLOT(show_dialog_assoc()));
}

void MainWindow::setup_menu_help()
{
    QMenu *help = menuBar()->addMenu("&Help");
    help->addAction("Con&tents", this, SLOT(show_help_content()), QKeySequence::HelpContents);
    help->addSeparator();
    QAction *act = help->addAction("A&bout", this, SLOT(show_dialog_about()));
    act->setMenuRole(QAction::AboutRole); // macOS
}

void MainWindow::setup_central_widget()
{
    wnd_file_ = new QPlainTextEdit(this);
    wnd_file_->setReadOnly(true);
    setCentralWidget(wnd_file_);
}

void MainWindow::steup_dock_widget()
{
    QDockWidget *dock_left = new QDockWidget("Result", this);
    dock_left->setFeatures(QDockWidget::NoDockWidgetFeatures);
    wnd_list_ = new QListWidget(dock_left);
    dock_left->setWidget(wnd_list_);

    QDockWidget *dock_bottom = new QDockWidget("Log", this);
    dock_bottom->setFeatures(QDockWidget::NoDockWidgetFeatures);
    wnd_log_ = new QPlainTextEdit(dock_bottom);
    wnd_log_->setReadOnly(true);
    wnd_log_->setMaximumBlockCount(par->log_size);
    wnd_log_->setStyleSheet("background-color: rgb(0,0,0); color: rgb(0,255,0);");
    dock_bottom->setWidget(wnd_log_);

    addDockWidget(Qt::LeftDockWidgetArea, dock_left);
    addDockWidget(Qt::BottomDockWidgetArea, dock_bottom);

    connect(wnd_list_, SIGNAL(currentRowChanged(int)), this, SLOT(on_wnd_list_currentRowChanged(int)));
}

void MainWindow::setup_status_bar()
{
    QStatusBar *sb = statusBar();

    status_ = new QLabel;
    sb->addWidget(status_);

    progress_ = new QProgressBar;
    progress_->setTextVisible(false);
    progress_->setMaximumHeight(sb->height() * 3 / 4);
    progress_->setMaximumWidth(200);
    sb->addPermanentWidget(progress_);
}

void MainWindow::create_process()
{
    process_ = new QProcess(this);

    process_->setProcessChannelMode(QProcess::MergedChannels);

    connect(process_, SIGNAL(readyReadStandardOutput()), this, SLOT(read_process_stdout()));

    connect(process_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(on_process_finished(int,QProcess::ExitStatus)));
}

void MainWindow::show_file_content(const QString &filename)
{
    qint64 maxlen = par->txt_size * 1024 * 1024;

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) {
        QMessageBox::critical(this, "Error - RTM-GWAS", "Can't open file: " + filename);
        return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);

    QTextStream in(&file);
    QString text = in.read(maxlen);

    wnd_file_->setPlainText(text);
    if (text.size() < file.size())
        wnd_file_->appendPlainText("\n*** Only part of the file is displayed");

    QApplication::restoreOverrideCursor();
}

QStringList MainWindow::get_process_output_files() const
{
    QStringList out;

    int pos = arguments_.indexOf("--out");

    if (pos != -1 && (pos + 1) < arguments_.size()) {
        QFileInfo fi(arguments_.at(pos + 1));
        QDir dir = fi.absoluteDir();
        dir.setNameFilters(QStringList(fi.fileName().append('*')));
        QStringList filenames = dir.entryList(QDir::Files | QDir::NoSymLinks | QDir::CaseSensitive);
        foreach (const QString &name, filenames)
            out << dir.absoluteFilePath(name);
    }

    return out;
}

bool MainWindow::is_process_running()
{
    if (process_->state() == QProcess::NotRunning)
        return false;

    QMessageBox::StandardButton ret =
            QMessageBox::question(this, "Terminate - RTM-GWAS",
                                  "Stop currently running computation?",
                                  QMessageBox::Ok | QMessageBox::Cancel);
    if (ret != QMessageBox::Ok)
        return true;

    disconnect(process_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(on_process_finished(int,QProcess::ExitStatus)));

    process_->close();

    status_->setText("Ready");
    progress_->setRange(0, 1);
    progress_->reset();

    foreach (const QString &name, get_process_output_files())
        QFile::remove(name);

    connect(process_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(on_process_finished(int,QProcess::ExitStatus)));

    return false;
}

void MainWindow::start_process()
{
    if (is_process_running())
        return;

    QString msg = "Start Process - ";
    msg.append(program_).append(' ').append(arguments_.join(" ")).append("\n\n");
    wnd_log_->moveCursor(QTextCursor::End);
    wnd_log_->appendPlainText(msg);
    wnd_log_->moveCursor(QTextCursor::End);

    process_->setWorkingDirectory(par->working_directory);
    process_->start(program_, arguments_);

    if (!process_->waitForStarted()) {
        QMessageBox::critical(this, "Error - RTM-GWAS", "Can't start process: " + program_);
        return;
    }

    status_->setText("Computing...");
    progress_->setRange(0, 0);
    progress_->reset();
}

void MainWindow::set_working_directory()
{
    QString dir = QFileDialog::getExistingDirectory(this, QString(), par->working_directory);
    if (!dir.isEmpty()) {
        par->working_directory = dir;
        par->file_dialog_directory = dir;
    }
}

void MainWindow::show_working_directory()
{
    QDesktopServices::openUrl(QUrl::fromLocalFile(par->working_directory));
}

void MainWindow::show_dialog_snpldb()
{
    DialogSNPLDB d(this);
    d.resize(640, 480);
    if (d.exec() == QDialog::Accepted) {
        program_ = d.program();
        arguments_ = d.arguments();
        start_process();
    }
}

void MainWindow::show_dialog_gsc()
{
    DialogGSC d(this);
    d.resize(640, 480);
    if (d.exec() == QDialog::Accepted) {
        program_ = d.program();
        arguments_ = d.arguments();
        start_process();
    }
}

void MainWindow::show_dialog_assoc()
{
    DialogAssoc d(this);
    d.resize(640, 480);
    if (d.exec() == QDialog::Accepted) {
        program_ = d.program();
        arguments_ = d.arguments();
        start_process();
    }
}

void MainWindow::show_help_content()
{
    QString pdf = "RTM-GWAS_UserGuide.pdf";
    QString www = "https://github.com/njau-sri/rtm-gwas/wiki";

    pdf = QDir(QApplication::applicationDirPath()).absoluteFilePath(pdf);
    if (QFile::exists(pdf)) {
        if (QDesktopServices::openUrl(QUrl::fromLocalFile(pdf)))
            return;
        QMessageBox::critical(this, "Error - RTM-GWAS", "Can't open " + pdf);
    }

    QDesktopServices::openUrl(QUrl(www));
}

void MainWindow::show_dialog_about()
{
    QString desc =
            "<h3>RTM-GWAS %1</h3>"
            "<p>Built on %2 %3</p>"
            "<p>Home: <a href=%4>%4</a></p>"
            "<p>Copyright (C) 2020 Nanjing Agricultural University</p>";

    QMessageBox::about(this, "About RTM-GWAS",
                       desc.arg(
                           RTM_GWAS_VERSION,
                           __DATE__,
                           __TIME__,
                           "https://github.com/njau-sri/rtm-gwas")
                       );
}

void MainWindow::on_wnd_list_currentRowChanged(int currentRow)
{
    if (currentRow >= 0 && currentRow < wnd_list_->count()) {
        QString filename = wnd_list_->item(currentRow)->toolTip();
        show_file_content(filename);
    }
    else
        wnd_file_->clear();
}

void MainWindow::read_process_stdout()
{
    QByteArray v = process_->readAllStandardOutput();
    wnd_log_->moveCursor(QTextCursor::End);
    wnd_log_->insertPlainText(QString::fromLocal8Bit(v.data(),v.size()));
    wnd_log_->moveCursor(QTextCursor::End);
}

void MainWindow::on_process_finished(int code, QProcess::ExitStatus status)
{
    status_->setText("Ready");
    progress_->setRange(0, 1);
    progress_->reset();

    QStringList filenames = get_process_output_files();

    if (status != QProcess::NormalExit || code != 0) {
        foreach (const QString &name, filenames)
            QFile::remove(name);
        QString msg = "Process exited unexpectedly: code %1, status %2";
        QMessageBox::critical(this, "Error - RTM-GWAS", msg.arg(code).arg(status));
        return;
    }

    foreach (const QString &name, filenames) {
        wnd_list_->addItem(QFileInfo(name).fileName());
        wnd_list_->item(wnd_list_->count() - 1)->setToolTip(name);
    }

    wnd_list_->setCurrentRow(wnd_list_->count() - 1);
}
