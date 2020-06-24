#include <QDir>
#include <QApplication>

#include "parameter.h"
#include "mainwindow.h"

namespace {

void handle_high_dpi_scaling()
{
#ifndef Q_OS_DARWIN
#if QT_VERSION >= QT_VERSION_CHECK(5,6,0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif
#if QT_VERSION >= QT_VERSION_CHECK(5,14,0)
    QGuiApplication::setHighDpiScaleFactorRoundingPolicy(Qt::HighDpiScaleFactorRoundingPolicy::Round);
#endif
#endif // Q_OS_DARWIN
}

void tweak_font()
{
#ifdef Q_OS_WIN
    QFont font = QApplication::font();
    if (font.family() == QLatin1String("SimSun")) {
        font.setFamily(QLatin1String("Microsoft YaHei"));
        QApplication::setFont(font);
    }
#endif // Q_OS_WIN
}

void set_bin_path()
{
    par->bin_path = QApplication::applicationDirPath();

#ifdef Q_OS_DARWIN
    if (par->bin_path.endsWith(".app/Contents/MacOS")) {
        QDir dir(par->bin_path);
        dir.cdUp();
        dir.cdUp();
        dir.cdUp();
        par->bin_path = dir.absolutePath();
    }
#endif // Q_OS_DARWIN
}

void set_working_directory()
{
    QDir dir = QDir::home();
#ifdef Q_OS_WIN
    dir.cd("Documents");
#endif
    dir.mkdir("RTM-GWAS");
    dir.cd("RTM-GWAS");

    par->working_directory = dir.absolutePath();
    par->file_dialog_directory = dir.absolutePath();
}

} // namespace

Parameter *par = nullptr;

int main(int argc, char *argv[])
{
    handle_high_dpi_scaling();

    QApplication a(argc, argv);

    Parameter parameter;
    par = &parameter;
    par->txt_size = 1;
    par->log_size = 1000;

    tweak_font();

    set_bin_path();

    set_working_directory();

    MainWindow w;
    w.resize(800, 600);
    w.show();

    return a.exec();
}
