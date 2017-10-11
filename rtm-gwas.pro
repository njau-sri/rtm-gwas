QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = rtm-gwas
TEMPLATE = app

CONFIG += static

DEFINES += QT_NO_CAST_FROM_ASCII QT_NO_CAST_TO_ASCII
DEFINES += QT_DEPRECATED_WARNINGS

SOURCES += main.cpp \
    mainwindow.cpp \
    parameter.cpp \
    dialogsnpldb.cpp \
    dialoggsc.cpp \
    dialogassoc.cpp

HEADERS += mainwindow.h \
    parameter.h \
    dialogsnpldb.h \
    dialoggsc.h \
    dialogassoc.h

FORMS += mainwindow.ui \
    dialogsnpldb.ui \
    dialoggsc.ui \
    dialogassoc.ui
