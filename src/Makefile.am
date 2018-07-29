CXXFLAGS = -std=c++11 -O2 -fopenmp

lib_LTLIBRARIES = libtylib.la libss.la

libtylib_la_SOURCES = bit.cpp observable_base.cpp ss.cpp observable.cpp translation.cpp free_fermion.cpp
libtylib_la_LDFLAGS = -version-info 0:0:0

libss_la_SOURCES = ss.cpp
libss_la_CXXFLAGS = -fopenmp
libss_la_LDFLAGS = -version-info 0:0:0 

include_HEADERS = log.h bit.h observable_base.h observable.h translation.h ss.h free_fermion.h
