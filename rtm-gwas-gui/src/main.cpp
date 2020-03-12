#include <QDir>
#include <QApplication>
#include "mainwindow.h"
#include "parameter.h"

int main(int argc, char *argv[])
{
    // handle high DPI scaling
#ifndef Q_OS_DARWIN
#if QT_VERSION >= QT_VERSION_CHECK(5,6,0)
    QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif
#if QT_VERSION >= QT_VERSION_CHECK(5,14,0)
    QGuiApplication::setHighDpiScaleFactorRoundingPolicy(Qt::HighDpiScaleFactorRoundingPolicy::Round);
#endif
#endif // Q_OS_DARWIN

    QApplication a(argc, argv);

    // tweak font for zh-CN
#ifdef Q_OS_WIN
    QFont font = QApplication::font();
    if (font.family() == QLatin1String("SimSun")) {
        font.setFamily(QLatin1String("Microsoft YaHei"));
        QApplication::setFont(font);
    }
#endif // Q_OS_WIN

    // set kernel path
    Parameter::exe = QApplication::applicationDirPath();
#ifdef Q_OS_DARWIN
    if (Parameter::exe.endsWith(QLatin1String(".app/Contents/MacOS"))) {
        QDir exe(Parameter::exe);
        exe.cdUp();
        exe.cdUp();
        exe.cdUp();
        Parameter::exe = exe.absolutePath();
    }
#endif // Q_OS_DARWIN

    // set working directory
    QDir home = QDir::home();
#ifdef Q_OS_WIN
    home.cd(QLatin1String("Documents"));
#endif
    home.mkdir(QLatin1String("RTM-GWAS"));
    home.cd(QLatin1String("RTM-GWAS"));

    Parameter::work = home.absolutePath();
    Parameter::open = Parameter::work;
    QDir::setCurrent(Parameter::work);

    MainWindow w;
    w.show();

    return a.exec();
}
