QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = rtm-gwas-gui
TEMPLATE = app

DEFINES += QT_NO_CAST_FROM_ASCII
DEFINES += QT_NO_CAST_TO_ASCII

RTM_GWAS_VERSION = $$(RTM_GWAS_VERSION)
isEmpty(RTM_GWAS_VERSION) {
    RTM_GWAS_VERSION = unknown
}
DEFINES += RTM_GWAS_VERSION=$$RTM_GWAS_VERSION

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
    parameter.h \
    version.h

FORMS += \
    mainwindow.ui \
    dialogassoc.ui \
    dialoggsc.ui \
    dialogsnpldb.ui
