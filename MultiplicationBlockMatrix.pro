######################################################################
# Automatically generated by qmake (2.01a) ?? ??? 30 09:04:21 2017
######################################################################

CONFIG += release

TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += "headers/" "src/"

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

# Input
SOURCES += \
    src/main.cpp \
    src/utils.cpp

HEADERS += \
    headers/utils.h
