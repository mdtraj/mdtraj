.. _faq:

Frequently Asked Questions
==========================

1. *Why am I seeing weird "symbol not found" errors while trying to import
   mdtraj?*

   Let me guess -- you're using a Linux cluster running CentOS 5? Your
   system compiler is gcc 4.1?

   You need a more recent version of gcc (>= 4.4). You can set the ``CC``
   and ``CXX`` environment variables to point to a different compiler. 

   Unfortunately, installing from conda will not help you here. You must
   compile from source on CentOS 5 (or similar).

   For Pande lab members on ``vsp-compute``, set the following environment
   variables and compile from source ::

       $ export CC=gcc44
       $ export CXX=g++44
       $ python setup.py install
