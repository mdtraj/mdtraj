/* -*- C -*-  (not really, but good for syntax highlighting) */
%module xtc

%{
  #include <xtc.hh>
  #define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
  import_array();
%}

/*This translates the C++ Error StopIteration into a python error */
%typemap(throws) StopIteration %{
  PyErr_SetNone(PyExc_StopIteration);
  SWIG_fail;
%}

%typemap(throws) TypeError %{
  PyErr_SetNone(PyExc_TypeError);
  SWIG_fail;
%}

%include "xtc.hh"
