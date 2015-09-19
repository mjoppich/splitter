#-------------------------------------------------
#
# Project created by QtCreator 2014-12-01T14:44:06
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = gtxtools-gui
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui

unix|win32: LIBS += -L$$PWD/../../../lib/ -lsplitutils -lgtx

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp

INCLUDEPATH += $$PWD/../../../include/
