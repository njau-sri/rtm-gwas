#include <QDir>
#include <QApplication>
#include <QDesktopWidget>
#include "mainwindow.h"
#include "params.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QDir::setCurrent(QDir::homePath());

    Params::exe = QDir(QApplication::applicationDirPath()).filePath(QLatin1String("rtm-gwas"));

    Params::work_dir = QDir::homePath();
    Params::open_dir = QDir::homePath();

    QRect s = QApplication::desktop()->screenGeometry();

    MainWindow w;

    w.resize(s.size()*0.8);

    int x = (s.width() - w.width()) / 2;
    int y = (s.height() - w.height()) / 2;
    w.move(x, y);

    w.show();

    return a.exec();
}
