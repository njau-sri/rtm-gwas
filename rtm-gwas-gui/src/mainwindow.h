#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMap>
#include <QProcess>
#include <QMainWindow>
#include <QProgressBar>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);

    ~MainWindow();

protected:

    void closeEvent(QCloseEvent *e);

private slots:

    void on_actionChangeDir_triggered();

    void on_actionSNPLDB_triggered();

    void on_actionGSC_triggered();

    void on_actionAssociation_triggered();

    void on_actionContents_triggered();

    void on_actionAbout_triggered();

    void on_listWidgetFile_currentRowChanged(int currentRow);

    void slot_proc_readyReadStandardOutput();

    void slot_proc_finished(int code, QProcess::ExitStatus status);

private:

    void loadFiles(const QStringList &fileNames);

    void showFile(const QString &fileName);

    bool isProcessRunning();

    QStringList getOutputFiles() const;

    void startProcess(const QString &prog, const QStringList &args);

private:

    Ui::MainWindow *ui;
    QProcess *proc_;
    QProgressBar *bar_;
    QStringList args_;
    QMap<QString,QString> mapfile_;
};

#endif // MAINWINDOW_H
