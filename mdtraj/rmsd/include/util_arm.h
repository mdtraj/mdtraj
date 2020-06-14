#include <arm_neon.h>


static INLINE void aos_interleaved_store(float* p, float32x4_t x, float32x4_t y, float32x4_t z) {    
    vst3q_f32(p, float32x4x3_t{x, y, z});
}

static INLINE float32x4x3_t aos_deinterleaved_load(const float* S) {
    return vld3q_f32(S);

}

#define REDUCTION_EPILOGUE(xx, xy, xz, yx, yy, yz, zx, zy, zz) \
  xx = vpaddq_f32(xx,xy); /* xx = xx01 xx23 xy01 xy23 */\
  xz = vpaddq_f32(xz,yx); /* xz = xz01 xz23 yx01 yx23 */\
  yy = vpaddq_f32(yy,yz); /* yy = yy01 yy23 yz01 yz23 */\
  zx = vpaddq_f32(zx,zy); /* zx = zx01 zx23 zy01 zy23 */\
  zz = vpaddq_f32(zz,zy); /* zz = zz01 zz23 zy01 zy23 */\
  xx = vpaddq_f32(xx,xz); /* xx = xx0123 xy0123 xz0123 yx0123 */\
  yy = vpaddq_f32(yy,zx); /* yy = yy0123 yz0123 zx0123 zy0123 */\
  zz = vpaddq_f32(zz,xz); /* zz = zz0123 zy0123 xz0123 yx0123 */  
