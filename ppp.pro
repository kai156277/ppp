#-------------------------------------------------
#
# Project created by QtCreator 2016-04-03T14:58:50
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = ppp
CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++11

TEMPLATE = app


SOURCES += main.cpp \
    o_date.cpp \
    file_read.cpp \
    sp3_date.cpp \
    clock_date.cpp \
    math_function.cpp \
    matching.cpp \
    readfilepath.cpp \
    gc_gpss.cpp \
    ppp_calculate.cpp \
    ppp_date.cpp \
    snx_date.cpp \
    coordinate.cpp

HEADERS += \
    o_date.h \
    file_read.h \
    sp3_date.h \
    clock_date.h \
    math_function.h \
    matching.h \
    readfilepath.h \
    gc_gpss.h \
    ppp_calculate.h \
    ppp_date.h \
    snx_date.h \
    coordinate.h
