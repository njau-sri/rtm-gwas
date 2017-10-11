#include <QDir>
#include <QApplication>
#include <QDesktopWidget>
#include "mainwindow.h"
#include "parameter.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QDir::setCurrent(QDir::homePath());

    Parameter::root = QApplication::applicationDirPath();
    Parameter::work = QDir::homePath();
    Parameter::open = QDir::homePath();

    QRect s = QApplication::desktop()->screenGeometry();

    MainWindow w;

    w.resize(s.size() * 0.75);

    int x = (s.width() - w.width()) / 2;
    int y = (s.height() - w.height()) / 2;
    w.move(x, y);

    w.show();

    return a.exec();
}
