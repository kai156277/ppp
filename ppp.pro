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
INCLUDEPATH += C:\QtLibrary\eigen

TEMPLATE = app


SOURCES += main.cpp \
    file_read.cpp \
    math_function.cpp \
    matching.cpp \
    readfilepath.cpp \
    ppp_calculate.cpp \
    antmod_data.cpp \
    clock_data.cpp \
    o_data.cpp \
    ppp_data.cpp \
    snx_data.cpp \
    sp3_data.cpp \
    time_system.cpp \
    coordinate_system.cpp \
    sun_moon_position.cpp \
    ocean_date.cpp \
    tide_correction.cpp \
    erp_data.cpp

HEADERS += \
    file_read.h \
    math_function.h \
    matching.h \
    readfilepath.h \
    ppp_calculate.h \
    antmod_data.h \
    clock_data.h \
    o_data.h \
    ppp_data.h \
    snx_data.h \
    sp3_data.h \
    time_system.h \
    coordinate_system.h \
    sun_moon_position.h \
    ocean_date.h \
    tide_correction.h \
    erp_data.h
