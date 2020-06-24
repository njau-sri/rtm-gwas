QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

DEFINES += QT_DEPRECATED_WARNINGS

RTM_GWAS_VERSION = $$cat(../VERSION)
DEFINES += "RTM_GWAS_VERSION=\\\"$$RTM_GWAS_VERSION\\\""

SOURCES += \
    src/dialogassoc.cpp \
    src/dialoggsc.cpp \
    src/dialogsnpldb.cpp \
    src/main.cpp \
    src/mainwindow.cpp

HEADERS += \
    src/dialogassoc.h \
    src/dialoggsc.h \
    src/dialogsnpldb.h \
    src/mainwindow.h \
    src/parameter.h
