#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QProcess>
#include <QMainWindow>
#include <QStringList>

class QLabel;
class QListWidget;
class QProgressBar;
class QPlainTextEdit;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);

protected:
    void closeEvent(QCloseEvent *e);

private:
    void setup_menu_file();
    void setup_menu_analysis();
    void setup_menu_tools();
    void setup_menu_help();
    void setup_central_widget();
    void steup_dock_widget();
    void setup_status_bar();
    void create_process();
    void show_file_content(const QString &filename);
    QStringList get_process_output_files() const;
    bool is_process_running();
    void start_process();

private slots:
    void show_dialog_gconv();
    void set_working_directory();
    void show_working_directory();
    void show_dialog_snpldb();
    void show_dialog_gsc();
    void show_dialog_assoc();
    void show_dialog_ld();
    void show_help_content();
    void show_dialog_about();
    void on_wnd_list_currentRowChanged(int currentRow);
    void read_process_stdout();
    void on_process_finished(int code, QProcess::ExitStatus status);

private:
    QListWidget *wnd_list_;
    QPlainTextEdit *wnd_file_;
    QPlainTextEdit *wnd_log_;
    QLabel *status_;
    QProgressBar *progress_;
    QProcess *process_;
    QString program_;
    QStringList arguments_;
};

#endif // MAINWINDOW_H
