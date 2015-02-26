/*
 * sse_swizzle.h: SSE "intrinsics" for SP float/32-bit int vector swizzling
 *                using SSE2 PSHUFD instruction
 * modified to add sse shuffle instructions too
 *
 * MIT License
 * Copyright 2011 Imran Haque (ihaque@cs.stanford.edu)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * This file defines SSE intrinsic-style functions for doing vector swizzling
 * according to the usual notation. For example, swizzles that in GLSL might
 * be written as:
 *
 *     float4 dst1 = src.xzwx;
 *
 * Can be written using this module as:
 *
 *     __m128 dst1 = _mm_swizzle_ps_xzwx(src);
 *
 * Similarly, a shuffle:
 *      
 *     float4 dst2 = src1.xy,src2.yz
 *
 * Becomes
 *
 *     __m128 dst2 = _mm_shuffle_ps_xyyz(src1,src2);
 *
 * Variants are defined for XYZW naming conventions only for floats.
 *
 */
#ifndef _SSE_SWIZZLE_H_
#define _SSE_SWIZZLE_H_
#include "msvccompat.h"
#include <emmintrin.h>
#ifdef __SSE2__
static INLINE __m128  _mm_swizzle_ps_xxxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x00)); }
static INLINE __m128  _mm_swizzle_ps_xxxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x40)); }
static INLINE __m128  _mm_swizzle_ps_xxxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x80)); }
static INLINE __m128  _mm_swizzle_ps_xxxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC0)); }
static INLINE __m128  _mm_swizzle_ps_xxyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x10)); }
static INLINE __m128  _mm_swizzle_ps_xxyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x50)); }
static INLINE __m128  _mm_swizzle_ps_xxyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x90)); }
static INLINE __m128  _mm_swizzle_ps_xxyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD0)); }
static INLINE __m128  _mm_swizzle_ps_xxzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x20)); }
static INLINE __m128  _mm_swizzle_ps_xxzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x60)); }
static INLINE __m128  _mm_swizzle_ps_xxzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA0)); }
static INLINE __m128  _mm_swizzle_ps_xxzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE0)); }
static INLINE __m128  _mm_swizzle_ps_xxwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x30)); }
static INLINE __m128  _mm_swizzle_ps_xxwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x70)); }
static INLINE __m128  _mm_swizzle_ps_xxwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB0)); }
static INLINE __m128  _mm_swizzle_ps_xxww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF0)); }
static INLINE __m128  _mm_swizzle_ps_xyxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x04)); }
static INLINE __m128  _mm_swizzle_ps_xyxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x44)); }
static INLINE __m128  _mm_swizzle_ps_xyxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x84)); }
static INLINE __m128  _mm_swizzle_ps_xyxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC4)); }
static INLINE __m128  _mm_swizzle_ps_xyyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x14)); }
static INLINE __m128  _mm_swizzle_ps_xyyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x54)); }
static INLINE __m128  _mm_swizzle_ps_xyyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x94)); }
static INLINE __m128  _mm_swizzle_ps_xyyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD4)); }
static INLINE __m128  _mm_swizzle_ps_xyzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x24)); }
static INLINE __m128  _mm_swizzle_ps_xyzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x64)); }
static INLINE __m128  _mm_swizzle_ps_xyzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA4)); }
static INLINE __m128  _mm_swizzle_ps_xyzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE4)); }
static INLINE __m128  _mm_swizzle_ps_xywx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x34)); }
static INLINE __m128  _mm_swizzle_ps_xywy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x74)); }
static INLINE __m128  _mm_swizzle_ps_xywz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB4)); }
static INLINE __m128  _mm_swizzle_ps_xyww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF4)); }
static INLINE __m128  _mm_swizzle_ps_xzxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x08)); }
static INLINE __m128  _mm_swizzle_ps_xzxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x48)); }
static INLINE __m128  _mm_swizzle_ps_xzxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x88)); }
static INLINE __m128  _mm_swizzle_ps_xzxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC8)); }
static INLINE __m128  _mm_swizzle_ps_xzyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x18)); }
static INLINE __m128  _mm_swizzle_ps_xzyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x58)); }
static INLINE __m128  _mm_swizzle_ps_xzyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x98)); }
static INLINE __m128  _mm_swizzle_ps_xzyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD8)); }
static INLINE __m128  _mm_swizzle_ps_xzzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x28)); }
static INLINE __m128  _mm_swizzle_ps_xzzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x68)); }
static INLINE __m128  _mm_swizzle_ps_xzzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA8)); }
static INLINE __m128  _mm_swizzle_ps_xzzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE8)); }
static INLINE __m128  _mm_swizzle_ps_xzwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x38)); }
static INLINE __m128  _mm_swizzle_ps_xzwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x78)); }
static INLINE __m128  _mm_swizzle_ps_xzwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB8)); }
static INLINE __m128  _mm_swizzle_ps_xzww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF8)); }
static INLINE __m128  _mm_swizzle_ps_xwxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x0C)); }
static INLINE __m128  _mm_swizzle_ps_xwxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x4C)); }
static INLINE __m128  _mm_swizzle_ps_xwxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x8C)); }
static INLINE __m128  _mm_swizzle_ps_xwxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xCC)); }
static INLINE __m128  _mm_swizzle_ps_xwyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x1C)); }
static INLINE __m128  _mm_swizzle_ps_xwyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x5C)); }
static INLINE __m128  _mm_swizzle_ps_xwyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x9C)); }
static INLINE __m128  _mm_swizzle_ps_xwyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xDC)); }
static INLINE __m128  _mm_swizzle_ps_xwzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x2C)); }
static INLINE __m128  _mm_swizzle_ps_xwzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x6C)); }
static INLINE __m128  _mm_swizzle_ps_xwzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xAC)); }
static INLINE __m128  _mm_swizzle_ps_xwzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xEC)); }
static INLINE __m128  _mm_swizzle_ps_xwwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x3C)); }
static INLINE __m128  _mm_swizzle_ps_xwwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x7C)); }
static INLINE __m128  _mm_swizzle_ps_xwwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xBC)); }
static INLINE __m128  _mm_swizzle_ps_xwww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xFC)); }
static INLINE __m128  _mm_swizzle_ps_yxxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x01)); }
static INLINE __m128  _mm_swizzle_ps_yxxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x41)); }
static INLINE __m128  _mm_swizzle_ps_yxxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x81)); }
static INLINE __m128  _mm_swizzle_ps_yxxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC1)); }
static INLINE __m128  _mm_swizzle_ps_yxyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x11)); }
static INLINE __m128  _mm_swizzle_ps_yxyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x51)); }
static INLINE __m128  _mm_swizzle_ps_yxyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x91)); }
static INLINE __m128  _mm_swizzle_ps_yxyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD1)); }
static INLINE __m128  _mm_swizzle_ps_yxzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x21)); }
static INLINE __m128  _mm_swizzle_ps_yxzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x61)); }
static INLINE __m128  _mm_swizzle_ps_yxzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA1)); }
static INLINE __m128  _mm_swizzle_ps_yxzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE1)); }
static INLINE __m128  _mm_swizzle_ps_yxwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x31)); }
static INLINE __m128  _mm_swizzle_ps_yxwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x71)); }
static INLINE __m128  _mm_swizzle_ps_yxwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB1)); }
static INLINE __m128  _mm_swizzle_ps_yxww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF1)); }
static INLINE __m128  _mm_swizzle_ps_yyxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x05)); }
static INLINE __m128  _mm_swizzle_ps_yyxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x45)); }
static INLINE __m128  _mm_swizzle_ps_yyxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x85)); }
static INLINE __m128  _mm_swizzle_ps_yyxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC5)); }
static INLINE __m128  _mm_swizzle_ps_yyyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x15)); }
static INLINE __m128  _mm_swizzle_ps_yyyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x55)); }
static INLINE __m128  _mm_swizzle_ps_yyyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x95)); }
static INLINE __m128  _mm_swizzle_ps_yyyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD5)); }
static INLINE __m128  _mm_swizzle_ps_yyzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x25)); }
static INLINE __m128  _mm_swizzle_ps_yyzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x65)); }
static INLINE __m128  _mm_swizzle_ps_yyzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA5)); }
static INLINE __m128  _mm_swizzle_ps_yyzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE5)); }
static INLINE __m128  _mm_swizzle_ps_yywx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x35)); }
static INLINE __m128  _mm_swizzle_ps_yywy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x75)); }
static INLINE __m128  _mm_swizzle_ps_yywz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB5)); }
static INLINE __m128  _mm_swizzle_ps_yyww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF5)); }
static INLINE __m128  _mm_swizzle_ps_yzxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x09)); }
static INLINE __m128  _mm_swizzle_ps_yzxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x49)); }
static INLINE __m128  _mm_swizzle_ps_yzxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x89)); }
static INLINE __m128  _mm_swizzle_ps_yzxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC9)); }
static INLINE __m128  _mm_swizzle_ps_yzyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x19)); }
static INLINE __m128  _mm_swizzle_ps_yzyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x59)); }
static INLINE __m128  _mm_swizzle_ps_yzyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x99)); }
static INLINE __m128  _mm_swizzle_ps_yzyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD9)); }
static INLINE __m128  _mm_swizzle_ps_yzzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x29)); }
static INLINE __m128  _mm_swizzle_ps_yzzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x69)); }
static INLINE __m128  _mm_swizzle_ps_yzzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA9)); }
static INLINE __m128  _mm_swizzle_ps_yzzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE9)); }
static INLINE __m128  _mm_swizzle_ps_yzwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x39)); }
static INLINE __m128  _mm_swizzle_ps_yzwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x79)); }
static INLINE __m128  _mm_swizzle_ps_yzwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB9)); }
static INLINE __m128  _mm_swizzle_ps_yzww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF9)); }
static INLINE __m128  _mm_swizzle_ps_ywxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x0D)); }
static INLINE __m128  _mm_swizzle_ps_ywxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x4D)); }
static INLINE __m128  _mm_swizzle_ps_ywxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x8D)); }
static INLINE __m128  _mm_swizzle_ps_ywxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xCD)); }
static INLINE __m128  _mm_swizzle_ps_ywyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x1D)); }
static INLINE __m128  _mm_swizzle_ps_ywyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x5D)); }
static INLINE __m128  _mm_swizzle_ps_ywyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x9D)); }
static INLINE __m128  _mm_swizzle_ps_ywyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xDD)); }
static INLINE __m128  _mm_swizzle_ps_ywzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x2D)); }
static INLINE __m128  _mm_swizzle_ps_ywzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x6D)); }
static INLINE __m128  _mm_swizzle_ps_ywzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xAD)); }
static INLINE __m128  _mm_swizzle_ps_ywzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xED)); }
static INLINE __m128  _mm_swizzle_ps_ywwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x3D)); }
static INLINE __m128  _mm_swizzle_ps_ywwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x7D)); }
static INLINE __m128  _mm_swizzle_ps_ywwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xBD)); }
static INLINE __m128  _mm_swizzle_ps_ywww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xFD)); }
static INLINE __m128  _mm_swizzle_ps_zxxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x02)); }
static INLINE __m128  _mm_swizzle_ps_zxxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x42)); }
static INLINE __m128  _mm_swizzle_ps_zxxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x82)); }
static INLINE __m128  _mm_swizzle_ps_zxxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC2)); }
static INLINE __m128  _mm_swizzle_ps_zxyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x12)); }
static INLINE __m128  _mm_swizzle_ps_zxyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x52)); }
static INLINE __m128  _mm_swizzle_ps_zxyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x92)); }
static INLINE __m128  _mm_swizzle_ps_zxyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD2)); }
static INLINE __m128  _mm_swizzle_ps_zxzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x22)); }
static INLINE __m128  _mm_swizzle_ps_zxzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x62)); }
static INLINE __m128  _mm_swizzle_ps_zxzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA2)); }
static INLINE __m128  _mm_swizzle_ps_zxzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE2)); }
static INLINE __m128  _mm_swizzle_ps_zxwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x32)); }
static INLINE __m128  _mm_swizzle_ps_zxwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x72)); }
static INLINE __m128  _mm_swizzle_ps_zxwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB2)); }
static INLINE __m128  _mm_swizzle_ps_zxww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF2)); }
static INLINE __m128  _mm_swizzle_ps_zyxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x06)); }
static INLINE __m128  _mm_swizzle_ps_zyxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x46)); }
static INLINE __m128  _mm_swizzle_ps_zyxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x86)); }
static INLINE __m128  _mm_swizzle_ps_zyxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC6)); }
static INLINE __m128  _mm_swizzle_ps_zyyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x16)); }
static INLINE __m128  _mm_swizzle_ps_zyyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x56)); }
static INLINE __m128  _mm_swizzle_ps_zyyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x96)); }
static INLINE __m128  _mm_swizzle_ps_zyyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD6)); }
static INLINE __m128  _mm_swizzle_ps_zyzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x26)); }
static INLINE __m128  _mm_swizzle_ps_zyzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x66)); }
static INLINE __m128  _mm_swizzle_ps_zyzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA6)); }
static INLINE __m128  _mm_swizzle_ps_zyzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE6)); }
static INLINE __m128  _mm_swizzle_ps_zywx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x36)); }
static INLINE __m128  _mm_swizzle_ps_zywy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x76)); }
static INLINE __m128  _mm_swizzle_ps_zywz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB6)); }
static INLINE __m128  _mm_swizzle_ps_zyww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF6)); }
static INLINE __m128  _mm_swizzle_ps_zzxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x0A)); }
static INLINE __m128  _mm_swizzle_ps_zzxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x4A)); }
static INLINE __m128  _mm_swizzle_ps_zzxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x8A)); }
static INLINE __m128  _mm_swizzle_ps_zzxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xCA)); }
static INLINE __m128  _mm_swizzle_ps_zzyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x1A)); }
static INLINE __m128  _mm_swizzle_ps_zzyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x5A)); }
static INLINE __m128  _mm_swizzle_ps_zzyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x9A)); }
static INLINE __m128  _mm_swizzle_ps_zzyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xDA)); }
static INLINE __m128  _mm_swizzle_ps_zzzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x2A)); }
static INLINE __m128  _mm_swizzle_ps_zzzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x6A)); }
static INLINE __m128  _mm_swizzle_ps_zzzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xAA)); }
static INLINE __m128  _mm_swizzle_ps_zzzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xEA)); }
static INLINE __m128  _mm_swizzle_ps_zzwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x3A)); }
static INLINE __m128  _mm_swizzle_ps_zzwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x7A)); }
static INLINE __m128  _mm_swizzle_ps_zzwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xBA)); }
static INLINE __m128  _mm_swizzle_ps_zzww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xFA)); }
static INLINE __m128  _mm_swizzle_ps_zwxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x0E)); }
static INLINE __m128  _mm_swizzle_ps_zwxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x4E)); }
static INLINE __m128  _mm_swizzle_ps_zwxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x8E)); }
static INLINE __m128  _mm_swizzle_ps_zwxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xCE)); }
static INLINE __m128  _mm_swizzle_ps_zwyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x1E)); }
static INLINE __m128  _mm_swizzle_ps_zwyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x5E)); }
static INLINE __m128  _mm_swizzle_ps_zwyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x9E)); }
static INLINE __m128  _mm_swizzle_ps_zwyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xDE)); }
static INLINE __m128  _mm_swizzle_ps_zwzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x2E)); }
static INLINE __m128  _mm_swizzle_ps_zwzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x6E)); }
static INLINE __m128  _mm_swizzle_ps_zwzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xAE)); }
static INLINE __m128  _mm_swizzle_ps_zwzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xEE)); }
static INLINE __m128  _mm_swizzle_ps_zwwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x3E)); }
static INLINE __m128  _mm_swizzle_ps_zwwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x7E)); }
static INLINE __m128  _mm_swizzle_ps_zwwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xBE)); }
static INLINE __m128  _mm_swizzle_ps_zwww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xFE)); }
static INLINE __m128  _mm_swizzle_ps_wxxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x03)); }
static INLINE __m128  _mm_swizzle_ps_wxxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x43)); }
static INLINE __m128  _mm_swizzle_ps_wxxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x83)); }
static INLINE __m128  _mm_swizzle_ps_wxxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC3)); }
static INLINE __m128  _mm_swizzle_ps_wxyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x13)); }
static INLINE __m128  _mm_swizzle_ps_wxyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x53)); }
static INLINE __m128  _mm_swizzle_ps_wxyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x93)); }
static INLINE __m128  _mm_swizzle_ps_wxyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD3)); }
static INLINE __m128  _mm_swizzle_ps_wxzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x23)); }
static INLINE __m128  _mm_swizzle_ps_wxzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x63)); }
static INLINE __m128  _mm_swizzle_ps_wxzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA3)); }
static INLINE __m128  _mm_swizzle_ps_wxzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE3)); }
static INLINE __m128  _mm_swizzle_ps_wxwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x33)); }
static INLINE __m128  _mm_swizzle_ps_wxwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x73)); }
static INLINE __m128  _mm_swizzle_ps_wxwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB3)); }
static INLINE __m128  _mm_swizzle_ps_wxww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF3)); }
static INLINE __m128  _mm_swizzle_ps_wyxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x07)); }
static INLINE __m128  _mm_swizzle_ps_wyxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x47)); }
static INLINE __m128  _mm_swizzle_ps_wyxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x87)); }
static INLINE __m128  _mm_swizzle_ps_wyxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xC7)); }
static INLINE __m128  _mm_swizzle_ps_wyyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x17)); }
static INLINE __m128  _mm_swizzle_ps_wyyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x57)); }
static INLINE __m128  _mm_swizzle_ps_wyyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x97)); }
static INLINE __m128  _mm_swizzle_ps_wyyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xD7)); }
static INLINE __m128  _mm_swizzle_ps_wyzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x27)); }
static INLINE __m128  _mm_swizzle_ps_wyzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x67)); }
static INLINE __m128  _mm_swizzle_ps_wyzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xA7)); }
static INLINE __m128  _mm_swizzle_ps_wyzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xE7)); }
static INLINE __m128  _mm_swizzle_ps_wywx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x37)); }
static INLINE __m128  _mm_swizzle_ps_wywy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x77)); }
static INLINE __m128  _mm_swizzle_ps_wywz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xB7)); }
static INLINE __m128  _mm_swizzle_ps_wyww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xF7)); }
static INLINE __m128  _mm_swizzle_ps_wzxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x0B)); }
static INLINE __m128  _mm_swizzle_ps_wzxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x4B)); }
static INLINE __m128  _mm_swizzle_ps_wzxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x8B)); }
static INLINE __m128  _mm_swizzle_ps_wzxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xCB)); }
static INLINE __m128  _mm_swizzle_ps_wzyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x1B)); }
static INLINE __m128  _mm_swizzle_ps_wzyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x5B)); }
static INLINE __m128  _mm_swizzle_ps_wzyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x9B)); }
static INLINE __m128  _mm_swizzle_ps_wzyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xDB)); }
static INLINE __m128  _mm_swizzle_ps_wzzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x2B)); }
static INLINE __m128  _mm_swizzle_ps_wzzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x6B)); }
static INLINE __m128  _mm_swizzle_ps_wzzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xAB)); }
static INLINE __m128  _mm_swizzle_ps_wzzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xEB)); }
static INLINE __m128  _mm_swizzle_ps_wzwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x3B)); }
static INLINE __m128  _mm_swizzle_ps_wzwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x7B)); }
static INLINE __m128  _mm_swizzle_ps_wzwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xBB)); }
static INLINE __m128  _mm_swizzle_ps_wzww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xFB)); }
static INLINE __m128  _mm_swizzle_ps_wwxx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x0F)); }
static INLINE __m128  _mm_swizzle_ps_wwxy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x4F)); }
static INLINE __m128  _mm_swizzle_ps_wwxz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x8F)); }
static INLINE __m128  _mm_swizzle_ps_wwxw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xCF)); }
static INLINE __m128  _mm_swizzle_ps_wwyx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x1F)); }
static INLINE __m128  _mm_swizzle_ps_wwyy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x5F)); }
static INLINE __m128  _mm_swizzle_ps_wwyz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x9F)); }
static INLINE __m128  _mm_swizzle_ps_wwyw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xDF)); }
static INLINE __m128  _mm_swizzle_ps_wwzx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x2F)); }
static INLINE __m128  _mm_swizzle_ps_wwzy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x6F)); }
static INLINE __m128  _mm_swizzle_ps_wwzz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xAF)); }
static INLINE __m128  _mm_swizzle_ps_wwzw(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xEF)); }
static INLINE __m128  _mm_swizzle_ps_wwwx(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x3F)); }
static INLINE __m128  _mm_swizzle_ps_wwwy(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0x7F)); }
static INLINE __m128  _mm_swizzle_ps_wwwz(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xBF)); }
static INLINE __m128  _mm_swizzle_ps_wwww(__m128 reg)      {  return CAST__M128(_mm_shuffle_epi32(CAST__M128I(reg), 0xFF)); }
static INLINE __m128  _mm_shuffle_ps_xxxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x00); }
static INLINE __m128  _mm_shuffle_ps_xxxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x40); }
static INLINE __m128  _mm_shuffle_ps_xxxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x80); }
static INLINE __m128  _mm_shuffle_ps_xxxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC0); }
static INLINE __m128  _mm_shuffle_ps_xxyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x10); }
static INLINE __m128  _mm_shuffle_ps_xxyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x50); }
static INLINE __m128  _mm_shuffle_ps_xxyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x90); }
static INLINE __m128  _mm_shuffle_ps_xxyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD0); }
static INLINE __m128  _mm_shuffle_ps_xxzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x20); }
static INLINE __m128  _mm_shuffle_ps_xxzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x60); }
static INLINE __m128  _mm_shuffle_ps_xxzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA0); }
static INLINE __m128  _mm_shuffle_ps_xxzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE0); }
static INLINE __m128  _mm_shuffle_ps_xxwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x30); }
static INLINE __m128  _mm_shuffle_ps_xxwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x70); }
static INLINE __m128  _mm_shuffle_ps_xxwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB0); }
static INLINE __m128  _mm_shuffle_ps_xxww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF0); }
static INLINE __m128  _mm_shuffle_ps_xyxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x04); }
static INLINE __m128  _mm_shuffle_ps_xyxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x44); }
static INLINE __m128  _mm_shuffle_ps_xyxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x84); }
static INLINE __m128  _mm_shuffle_ps_xyxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC4); }
static INLINE __m128  _mm_shuffle_ps_xyyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x14); }
static INLINE __m128  _mm_shuffle_ps_xyyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x54); }
static INLINE __m128  _mm_shuffle_ps_xyyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x94); }
static INLINE __m128  _mm_shuffle_ps_xyyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD4); }
static INLINE __m128  _mm_shuffle_ps_xyzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x24); }
static INLINE __m128  _mm_shuffle_ps_xyzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x64); }
static INLINE __m128  _mm_shuffle_ps_xyzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA4); }
static INLINE __m128  _mm_shuffle_ps_xyzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE4); }
static INLINE __m128  _mm_shuffle_ps_xywx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x34); }
static INLINE __m128  _mm_shuffle_ps_xywy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x74); }
static INLINE __m128  _mm_shuffle_ps_xywz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB4); }
static INLINE __m128  _mm_shuffle_ps_xyww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF4); }
static INLINE __m128  _mm_shuffle_ps_xzxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x08); }
static INLINE __m128  _mm_shuffle_ps_xzxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x48); }
static INLINE __m128  _mm_shuffle_ps_xzxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x88); }
static INLINE __m128  _mm_shuffle_ps_xzxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC8); }
static INLINE __m128  _mm_shuffle_ps_xzyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x18); }
static INLINE __m128  _mm_shuffle_ps_xzyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x58); }
static INLINE __m128  _mm_shuffle_ps_xzyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x98); }
static INLINE __m128  _mm_shuffle_ps_xzyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD8); }
static INLINE __m128  _mm_shuffle_ps_xzzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x28); }
static INLINE __m128  _mm_shuffle_ps_xzzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x68); }
static INLINE __m128  _mm_shuffle_ps_xzzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA8); }
static INLINE __m128  _mm_shuffle_ps_xzzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE8); }
static INLINE __m128  _mm_shuffle_ps_xzwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x38); }
static INLINE __m128  _mm_shuffle_ps_xzwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x78); }
static INLINE __m128  _mm_shuffle_ps_xzwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB8); }
static INLINE __m128  _mm_shuffle_ps_xzww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF8); }
static INLINE __m128  _mm_shuffle_ps_xwxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x0C); }
static INLINE __m128  _mm_shuffle_ps_xwxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x4C); }
static INLINE __m128  _mm_shuffle_ps_xwxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x8C); }
static INLINE __m128  _mm_shuffle_ps_xwxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xCC); }
static INLINE __m128  _mm_shuffle_ps_xwyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x1C); }
static INLINE __m128  _mm_shuffle_ps_xwyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x5C); }
static INLINE __m128  _mm_shuffle_ps_xwyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x9C); }
static INLINE __m128  _mm_shuffle_ps_xwyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xDC); }
static INLINE __m128  _mm_shuffle_ps_xwzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x2C); }
static INLINE __m128  _mm_shuffle_ps_xwzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x6C); }
static INLINE __m128  _mm_shuffle_ps_xwzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xAC); }
static INLINE __m128  _mm_shuffle_ps_xwzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xEC); }
static INLINE __m128  _mm_shuffle_ps_xwwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x3C); }
static INLINE __m128  _mm_shuffle_ps_xwwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x7C); }
static INLINE __m128  _mm_shuffle_ps_xwwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xBC); }
static INLINE __m128  _mm_shuffle_ps_xwww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xFC); }
static INLINE __m128  _mm_shuffle_ps_yxxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x01); }
static INLINE __m128  _mm_shuffle_ps_yxxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x41); }
static INLINE __m128  _mm_shuffle_ps_yxxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x81); }
static INLINE __m128  _mm_shuffle_ps_yxxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC1); }
static INLINE __m128  _mm_shuffle_ps_yxyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x11); }
static INLINE __m128  _mm_shuffle_ps_yxyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x51); }
static INLINE __m128  _mm_shuffle_ps_yxyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x91); }
static INLINE __m128  _mm_shuffle_ps_yxyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD1); }
static INLINE __m128  _mm_shuffle_ps_yxzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x21); }
static INLINE __m128  _mm_shuffle_ps_yxzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x61); }
static INLINE __m128  _mm_shuffle_ps_yxzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA1); }
static INLINE __m128  _mm_shuffle_ps_yxzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE1); }
static INLINE __m128  _mm_shuffle_ps_yxwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x31); }
static INLINE __m128  _mm_shuffle_ps_yxwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x71); }
static INLINE __m128  _mm_shuffle_ps_yxwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB1); }
static INLINE __m128  _mm_shuffle_ps_yxww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF1); }
static INLINE __m128  _mm_shuffle_ps_yyxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x05); }
static INLINE __m128  _mm_shuffle_ps_yyxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x45); }
static INLINE __m128  _mm_shuffle_ps_yyxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x85); }
static INLINE __m128  _mm_shuffle_ps_yyxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC5); }
static INLINE __m128  _mm_shuffle_ps_yyyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x15); }
static INLINE __m128  _mm_shuffle_ps_yyyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x55); }
static INLINE __m128  _mm_shuffle_ps_yyyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x95); }
static INLINE __m128  _mm_shuffle_ps_yyyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD5); }
static INLINE __m128  _mm_shuffle_ps_yyzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x25); }
static INLINE __m128  _mm_shuffle_ps_yyzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x65); }
static INLINE __m128  _mm_shuffle_ps_yyzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA5); }
static INLINE __m128  _mm_shuffle_ps_yyzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE5); }
static INLINE __m128  _mm_shuffle_ps_yywx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x35); }
static INLINE __m128  _mm_shuffle_ps_yywy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x75); }
static INLINE __m128  _mm_shuffle_ps_yywz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB5); }
static INLINE __m128  _mm_shuffle_ps_yyww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF5); }
static INLINE __m128  _mm_shuffle_ps_yzxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x09); }
static INLINE __m128  _mm_shuffle_ps_yzxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x49); }
static INLINE __m128  _mm_shuffle_ps_yzxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x89); }
static INLINE __m128  _mm_shuffle_ps_yzxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC9); }
static INLINE __m128  _mm_shuffle_ps_yzyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x19); }
static INLINE __m128  _mm_shuffle_ps_yzyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x59); }
static INLINE __m128  _mm_shuffle_ps_yzyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x99); }
static INLINE __m128  _mm_shuffle_ps_yzyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD9); }
static INLINE __m128  _mm_shuffle_ps_yzzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x29); }
static INLINE __m128  _mm_shuffle_ps_yzzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x69); }
static INLINE __m128  _mm_shuffle_ps_yzzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA9); }
static INLINE __m128  _mm_shuffle_ps_yzzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE9); }
static INLINE __m128  _mm_shuffle_ps_yzwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x39); }
static INLINE __m128  _mm_shuffle_ps_yzwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x79); }
static INLINE __m128  _mm_shuffle_ps_yzwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB9); }
static INLINE __m128  _mm_shuffle_ps_yzww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF9); }
static INLINE __m128  _mm_shuffle_ps_ywxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x0D); }
static INLINE __m128  _mm_shuffle_ps_ywxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x4D); }
static INLINE __m128  _mm_shuffle_ps_ywxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x8D); }
static INLINE __m128  _mm_shuffle_ps_ywxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xCD); }
static INLINE __m128  _mm_shuffle_ps_ywyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x1D); }
static INLINE __m128  _mm_shuffle_ps_ywyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x5D); }
static INLINE __m128  _mm_shuffle_ps_ywyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x9D); }
static INLINE __m128  _mm_shuffle_ps_ywyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xDD); }
static INLINE __m128  _mm_shuffle_ps_ywzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x2D); }
static INLINE __m128  _mm_shuffle_ps_ywzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x6D); }
static INLINE __m128  _mm_shuffle_ps_ywzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xAD); }
static INLINE __m128  _mm_shuffle_ps_ywzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xED); }
static INLINE __m128  _mm_shuffle_ps_ywwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x3D); }
static INLINE __m128  _mm_shuffle_ps_ywwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x7D); }
static INLINE __m128  _mm_shuffle_ps_ywwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xBD); }
static INLINE __m128  _mm_shuffle_ps_ywww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xFD); }
static INLINE __m128  _mm_shuffle_ps_zxxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x02); }
static INLINE __m128  _mm_shuffle_ps_zxxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x42); }
static INLINE __m128  _mm_shuffle_ps_zxxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x82); }
static INLINE __m128  _mm_shuffle_ps_zxxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC2); }
static INLINE __m128  _mm_shuffle_ps_zxyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x12); }
static INLINE __m128  _mm_shuffle_ps_zxyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x52); }
static INLINE __m128  _mm_shuffle_ps_zxyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x92); }
static INLINE __m128  _mm_shuffle_ps_zxyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD2); }
static INLINE __m128  _mm_shuffle_ps_zxzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x22); }
static INLINE __m128  _mm_shuffle_ps_zxzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x62); }
static INLINE __m128  _mm_shuffle_ps_zxzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA2); }
static INLINE __m128  _mm_shuffle_ps_zxzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE2); }
static INLINE __m128  _mm_shuffle_ps_zxwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x32); }
static INLINE __m128  _mm_shuffle_ps_zxwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x72); }
static INLINE __m128  _mm_shuffle_ps_zxwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB2); }
static INLINE __m128  _mm_shuffle_ps_zxww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF2); }
static INLINE __m128  _mm_shuffle_ps_zyxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x06); }
static INLINE __m128  _mm_shuffle_ps_zyxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x46); }
static INLINE __m128  _mm_shuffle_ps_zyxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x86); }
static INLINE __m128  _mm_shuffle_ps_zyxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC6); }
static INLINE __m128  _mm_shuffle_ps_zyyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x16); }
static INLINE __m128  _mm_shuffle_ps_zyyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x56); }
static INLINE __m128  _mm_shuffle_ps_zyyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x96); }
static INLINE __m128  _mm_shuffle_ps_zyyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD6); }
static INLINE __m128  _mm_shuffle_ps_zyzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x26); }
static INLINE __m128  _mm_shuffle_ps_zyzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x66); }
static INLINE __m128  _mm_shuffle_ps_zyzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA6); }
static INLINE __m128  _mm_shuffle_ps_zyzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE6); }
static INLINE __m128  _mm_shuffle_ps_zywx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x36); }
static INLINE __m128  _mm_shuffle_ps_zywy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x76); }
static INLINE __m128  _mm_shuffle_ps_zywz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB6); }
static INLINE __m128  _mm_shuffle_ps_zyww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF6); }
static INLINE __m128  _mm_shuffle_ps_zzxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x0A); }
static INLINE __m128  _mm_shuffle_ps_zzxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x4A); }
static INLINE __m128  _mm_shuffle_ps_zzxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x8A); }
static INLINE __m128  _mm_shuffle_ps_zzxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xCA); }
static INLINE __m128  _mm_shuffle_ps_zzyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x1A); }
static INLINE __m128  _mm_shuffle_ps_zzyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x5A); }
static INLINE __m128  _mm_shuffle_ps_zzyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x9A); }
static INLINE __m128  _mm_shuffle_ps_zzyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xDA); }
static INLINE __m128  _mm_shuffle_ps_zzzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x2A); }
static INLINE __m128  _mm_shuffle_ps_zzzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x6A); }
static INLINE __m128  _mm_shuffle_ps_zzzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xAA); }
static INLINE __m128  _mm_shuffle_ps_zzzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xEA); }
static INLINE __m128  _mm_shuffle_ps_zzwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x3A); }
static INLINE __m128  _mm_shuffle_ps_zzwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x7A); }
static INLINE __m128  _mm_shuffle_ps_zzwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xBA); }
static INLINE __m128  _mm_shuffle_ps_zzww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xFA); }
static INLINE __m128  _mm_shuffle_ps_zwxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x0E); }
static INLINE __m128  _mm_shuffle_ps_zwxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x4E); }
static INLINE __m128  _mm_shuffle_ps_zwxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x8E); }
static INLINE __m128  _mm_shuffle_ps_zwxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xCE); }
static INLINE __m128  _mm_shuffle_ps_zwyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x1E); }
static INLINE __m128  _mm_shuffle_ps_zwyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x5E); }
static INLINE __m128  _mm_shuffle_ps_zwyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x9E); }
static INLINE __m128  _mm_shuffle_ps_zwyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xDE); }
static INLINE __m128  _mm_shuffle_ps_zwzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x2E); }
static INLINE __m128  _mm_shuffle_ps_zwzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x6E); }
static INLINE __m128  _mm_shuffle_ps_zwzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xAE); }
static INLINE __m128  _mm_shuffle_ps_zwzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xEE); }
static INLINE __m128  _mm_shuffle_ps_zwwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x3E); }
static INLINE __m128  _mm_shuffle_ps_zwwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x7E); }
static INLINE __m128  _mm_shuffle_ps_zwwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xBE); }
static INLINE __m128  _mm_shuffle_ps_zwww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xFE); }
static INLINE __m128  _mm_shuffle_ps_wxxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x03); }
static INLINE __m128  _mm_shuffle_ps_wxxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x43); }
static INLINE __m128  _mm_shuffle_ps_wxxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x83); }
static INLINE __m128  _mm_shuffle_ps_wxxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC3); }
static INLINE __m128  _mm_shuffle_ps_wxyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x13); }
static INLINE __m128  _mm_shuffle_ps_wxyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x53); }
static INLINE __m128  _mm_shuffle_ps_wxyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x93); }
static INLINE __m128  _mm_shuffle_ps_wxyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD3); }
static INLINE __m128  _mm_shuffle_ps_wxzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x23); }
static INLINE __m128  _mm_shuffle_ps_wxzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x63); }
static INLINE __m128  _mm_shuffle_ps_wxzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA3); }
static INLINE __m128  _mm_shuffle_ps_wxzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE3); }
static INLINE __m128  _mm_shuffle_ps_wxwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x33); }
static INLINE __m128  _mm_shuffle_ps_wxwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x73); }
static INLINE __m128  _mm_shuffle_ps_wxwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB3); }
static INLINE __m128  _mm_shuffle_ps_wxww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF3); }
static INLINE __m128  _mm_shuffle_ps_wyxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x07); }
static INLINE __m128  _mm_shuffle_ps_wyxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x47); }
static INLINE __m128  _mm_shuffle_ps_wyxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x87); }
static INLINE __m128  _mm_shuffle_ps_wyxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xC7); }
static INLINE __m128  _mm_shuffle_ps_wyyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x17); }
static INLINE __m128  _mm_shuffle_ps_wyyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x57); }
static INLINE __m128  _mm_shuffle_ps_wyyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x97); }
static INLINE __m128  _mm_shuffle_ps_wyyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xD7); }
static INLINE __m128  _mm_shuffle_ps_wyzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x27); }
static INLINE __m128  _mm_shuffle_ps_wyzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x67); }
static INLINE __m128  _mm_shuffle_ps_wyzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xA7); }
static INLINE __m128  _mm_shuffle_ps_wyzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xE7); }
static INLINE __m128  _mm_shuffle_ps_wywx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x37); }
static INLINE __m128  _mm_shuffle_ps_wywy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x77); }
static INLINE __m128  _mm_shuffle_ps_wywz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xB7); }
static INLINE __m128  _mm_shuffle_ps_wyww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xF7); }
static INLINE __m128  _mm_shuffle_ps_wzxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x0B); }
static INLINE __m128  _mm_shuffle_ps_wzxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x4B); }
static INLINE __m128  _mm_shuffle_ps_wzxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x8B); }
static INLINE __m128  _mm_shuffle_ps_wzxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xCB); }
static INLINE __m128  _mm_shuffle_ps_wzyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x1B); }
static INLINE __m128  _mm_shuffle_ps_wzyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x5B); }
static INLINE __m128  _mm_shuffle_ps_wzyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x9B); }
static INLINE __m128  _mm_shuffle_ps_wzyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xDB); }
static INLINE __m128  _mm_shuffle_ps_wzzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x2B); }
static INLINE __m128  _mm_shuffle_ps_wzzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x6B); }
static INLINE __m128  _mm_shuffle_ps_wzzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xAB); }
static INLINE __m128  _mm_shuffle_ps_wzzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xEB); }
static INLINE __m128  _mm_shuffle_ps_wzwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x3B); }
static INLINE __m128  _mm_shuffle_ps_wzwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x7B); }
static INLINE __m128  _mm_shuffle_ps_wzwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xBB); }
static INLINE __m128  _mm_shuffle_ps_wzww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xFB); }
static INLINE __m128  _mm_shuffle_ps_wwxx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x0F); }
static INLINE __m128  _mm_shuffle_ps_wwxy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x4F); }
static INLINE __m128  _mm_shuffle_ps_wwxz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x8F); }
static INLINE __m128  _mm_shuffle_ps_wwxw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xCF); }
static INLINE __m128  _mm_shuffle_ps_wwyx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x1F); }
static INLINE __m128  _mm_shuffle_ps_wwyy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x5F); }
static INLINE __m128  _mm_shuffle_ps_wwyz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x9F); }
static INLINE __m128  _mm_shuffle_ps_wwyw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xDF); }
static INLINE __m128  _mm_shuffle_ps_wwzx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x2F); }
static INLINE __m128  _mm_shuffle_ps_wwzy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x6F); }
static INLINE __m128  _mm_shuffle_ps_wwzz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xAF); }
static INLINE __m128  _mm_shuffle_ps_wwzw(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xEF); }
static INLINE __m128  _mm_shuffle_ps_wwwx(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x3F); }
static INLINE __m128  _mm_shuffle_ps_wwwy(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0x7F); }
static INLINE __m128  _mm_shuffle_ps_wwwz(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xBF); }
static INLINE __m128  _mm_shuffle_ps_wwww(__m128 r1, __m128 r2)      {  return _mm_shuffle_ps(r1,r2,0xFF); }
#endif
#endif
