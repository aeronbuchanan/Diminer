TEMPLATE = app
CONFIG += console
CONFIG -= qt
CONFIG += x11
CONFIG += thread

QMAKE_CXXFLAGS += -std=c++11

INCLUDEPATH += "../../CImg-1.5.4/"

SOURCES += main.cpp \
    inpainters.cpp

HEADERS += \
    inpainters.h \
    patch.h \
    patch.h \
    texSynth.h \
    texSynth.hpp \
    vecn.h \
    table.h \
    common.h \
    circularSeam.h

