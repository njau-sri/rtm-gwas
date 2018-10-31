#-------------------------------------------------
#
# Project created by QtCreator 2018-10-31T11:59:38
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = rtm-gwas
TEMPLATE = app

DEFINES += QT_NO_CAST_FROM_ASCII
DEFINES += QT_NO_CAST_TO_ASCII

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    dialogassoc.cpp \
    dialoggsc.cpp \
    dialogsnpldb.cpp \
    parameter.cpp

HEADERS += \
        mainwindow.h \
    dialogassoc.h \
    dialoggsc.h \
    dialogsnpldb.h \
    parameter.h

FORMS += \
        mainwindow.ui \
    dialogassoc.ui \
    dialoggsc.ui \
    dialogsnpldb.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
