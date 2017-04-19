QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = rtm-gwas-gui
TEMPLATE = app
CONFIG += static

DEFINES += QT_NO_CAST_FROM_ASCII QT_NO_CAST_TO_ASCII

SOURCES += main.cpp\
        mainwindow.cpp \
    dialogassoc.cpp \
    dialogdata.cpp \
    dialogeigen.cpp \
    dialogldb.cpp \
    params.cpp \
    dialogsummary.cpp

HEADERS  += mainwindow.h \
    dialogassoc.h \
    dialogdata.h \
    dialogeigen.h \
    dialogldb.h \
    params.h \
    dialogsummary.h

FORMS    += mainwindow.ui \
    dialogassoc.ui \
    dialogdata.ui \
    dialogeigen.ui \
    dialogldb.ui \
    dialogsummary.ui
