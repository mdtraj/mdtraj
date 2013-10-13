#ifndef __MSVC_COMPAT_H__
#define __MSVC_COMPAT_H__

#ifdef _MSC_VER
# ifdef _M_IX86_FP
#  if _M_IX86_FP >= 1
#   ifndef __SSE__
#    define __SSE__ 1
#   endif
#  endif
#  if _M_IX86_FP >= 2
#   ifndef __SSE2__
#    define __SSE2__ 1
#   endif
#  endif
# elif defined(_M_AMD64)
// If the target is x86_64 then SSE2 is guaranteed
#  ifndef __SSE__
#   define __SSE__ 1
#  endif
#  ifndef __SSE2__
#   define __SSE2__ 1
#  endif
# endif
#endif

#ifndef INLINE
  #if defined(__GNUC__)
    #define INLINE __inline__
  #elif defined(_MSC_VER)
    #define INLINE __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define INLINE inline
  #else
    #define INLINE
  #endif
#endif

/**
 * Bitwise casts of SSE registers
 * http://stackoverflow.com/questions/13631951/bitwise-cast-from-m128-to-m128i-on-msvc/13632812#13632812
 **/
#ifdef _MSC_VER
 #define CAST__M128(x)  ( _mm_castsi128_ps(x) )
 #define CAST__M128I(x) ( _mm_castps_si128(x) )
#else
 #define CAST__M128(x)  ( (__m128)  x )
 #define CAST__M128I(x) ( (__m128i) x )
#endif

/**
 * Alignment of stack variables
 * http://stackoverflow.com/questions/7895869/cross-platform-alignx-macro
 **/
#if defined(_MSC_VER)
#define _ALIGNED(x) __declspec(align(x))
#else
#if defined(__GNUC__)
#define _ALIGNED(x) __attribute__ ((aligned(x)))
#endif
#endif



#endif
