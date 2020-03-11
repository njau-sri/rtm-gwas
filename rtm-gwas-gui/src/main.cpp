#include <QDir>
#include <QApplication>
#include "mainwindow.h"
#include "parameter.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QFont font = QApplication::font();
    if (font.defaultFamily() == QLatin1String("SimSun")) {
        font.setFamily(QLatin1String("Microsoft YaHei"));
        QApplication::setFont(font);
    }

    Parameter::exe = QApplication::applicationDirPath();
#ifdef Q_OS_DARWIN
    if (Parameter::exe.endsWith(QLatin1String(".app/Contents/MacOS"))) {
        QDir exe(Parameter::exe);
        exe.cdUp();
        exe.cdUp();
        exe.cdUp();
        Parameter::exe = exe.absolutePath();
    }
#endif

    QDir home = QDir::home();
#ifdef Q_OS_WIN32
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
