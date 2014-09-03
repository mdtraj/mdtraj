.. _faq:

Frequently Asked Questions
==========================

1. Why am I seeing an error message about my CPU and SSE 4.1?::

     RuntimeError: This CPU does not support the required instruction set (SSE4.1)
    
   Let me guess -- you're using a Linux cluster running CentOS 5? Your system
   compiler is gcc 4.1? Some functions in MDTraj require specific CPU features
   (and compiler support for them) to enable faster performance.
   
   It's your machine actually DOES support SSE4.1 instructions, but your compiler
   is too old to emit them. Try installing MDTraj from the prebuilt
   :ref:`conda <install-with-conda>` packages, or getting a more recent
   compiler. (When building MDTraj from source, you can set the ``CC``
   environment variable to point to a different compiler.)
   
