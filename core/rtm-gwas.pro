TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -static

LIBS += -Le:/devel/deps/OpenBLAS-win32 -lopenblas -lgfortran -lquadmath -lgomp
LIBS += -Le:/devel/deps/Rmath-win32 -lRmath

SOURCES += \
    src/appassoc.cpp \
    src/appdata.cpp \
    src/appgsc.cpp \
    src/appldb.cpp \
    src/appsum.cpp \
    src/assoclm.cpp \
    src/assoclmm.cpp \
    src/assocrtm.cpp \
    src/cmdline.cpp \
    src/emma.cpp \
    src/glm.cpp \
    src/hapmapio.cpp \
    src/lapack.cpp \
    src/lstsqr.cpp \
    src/main.cpp \
    src/plinkio.cpp \
    src/rtmio.cpp \
    src/stat.cpp \
    src/stepreg.cpp \
    src/util.cpp \
    src/vcfio.cpp \
    src/zeroin.c

HEADERS += \
    src/appassoc.h \
    src/appdata.h \
    src/appgsc.h \
    src/appldb.h \
    src/appsum.h \
    src/assoclm.h \
    src/assoclmm.h \
    src/assocrtm.h \
    src/cmdline.h \
    src/emma.h \
    src/glm.h \
    src/haploprob.h \
    src/hapmapio.h \
    src/lapack.h \
    src/lstsqr.h \
    src/main.h \
    src/plinkio.h \
    src/rtmio.h \
    src/stat.h \
    src/stepreg.h \
    src/strsplit.h \
    src/util.h \
    src/vcfio.h
