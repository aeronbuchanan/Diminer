TEMPLATE = app
CONFIG += console
CONFIG -= qt
CONFIG += x11
CONFIG += thread

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += "../CImg/"

SOURCES += main.cpp \
    inpainters.cpp \
    diminer.cpp \
    boundaryChains.cpp \

HEADERS += \
    inpainters.h \
    patch.h \
    patch.h \
    texSynth.h \
    texSynth.hpp \
    vecn.h \
    table.h \
    dinimer.h \
    circularSeam.h \
    boundaryChains.h \
    smoothGradient.hpp


