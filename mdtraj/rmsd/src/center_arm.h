#include "util_arm.h"
#include <cstddef>

void inplace_center_and_trace_atom_major(float* coords, float* traces, const int n_frames, const int n_atoms)
{
    /* Center a trajectory containing multiple conformations inplace.
       The coordinates are store in float, but the accumulation is done in
       double.

       Also compute the traces of the centered conformations, which are necessary
       for RMSD.
    */
    long long i, k;
    float* confp;
    float64x2_t sx_, sy_, sz_, trace_;
    float32x4_t mux_, muy_, muz_;
    float sxf, syf, szf;
    double sx, sy, sz, trace;
    float32x4x3_t xyz;
    float32x4x3_t xyz2;

    #ifdef _OPENMP
    #pragma omp parallel for shared(coords, traces) \
        private(sx_, sy_, sz_, trace_, mux_, muy_, muz_, sxf, syf, szf, \
        confp, i, xyz, xyz2, sx, sy, sz, trace)
    #endif
    for (k = 0; k < n_frames; k++) {
        confp = &coords[k * n_atoms * 3];
        sx_ = sy_ = sz_ = trace_ =  vdupq_n_f64(0);
        for (i = 0; i < n_atoms/4; i++) {
            xyz = aos_deinterleaved_load(confp);

            /* accumulate the sums of each coordinate in double */
            /* get the first two values from each float4 */
            sx_ = vaddq_f64(sx_, vcvt_f64_f32(vget_low_f32(xyz.val[0])));
            sy_ = vaddq_f64(sy_, vcvt_f64_f32(vget_low_f32(xyz.val[1])));
            sz_ = vaddq_f64(sz_, vcvt_f64_f32(vget_low_f32(xyz.val[2])));
            
            /* and the second two values from each float4 */
            sx_ = vaddq_f64(sx_, vcvt_f64_f32(vget_high_f32(xyz.val[0])));
            sy_ = vaddq_f64(sy_, vcvt_f64_f32(vget_high_f32(xyz.val[1])));
            sz_ = vaddq_f64(sz_, vcvt_f64_f32(vget_high_f32(xyz.val[2])));
            confp += 12;
        }
        /* copy the summed coordinates out of the SSE registers */
        sx = vgetq_lane_f64(sx_, 0) + vgetq_lane_f64(sx_, 1);
        sy = vgetq_lane_f64(sy_, 0) + vgetq_lane_f64(sy_, 1);
        sz = vgetq_lane_f64(sz_, 0) + vgetq_lane_f64(sz_, 1);
        
        /* Add the last couple entries that weren't a factor of four */
        for (i = 0; i < n_atoms % 4; i++) {
            sx += confp[i*3 + 0];
            sy += confp[i*3 + 1];
            sz += confp[i*3 + 2];
        }

        /* Now we want mean x, y, and z positions */
        sx /= n_atoms;
        sy /= n_atoms;
        sz /= n_atoms;

        /* Load these mean positions back into the SSE registers */
        sxf = (float) sx;
        syf = (float) sy;
        szf = (float) sz;
        mux_ = vld1q_dup_f32(&sxf);
        muy_ = vld1q_dup_f32(&syf);
        muz_ = vld1q_dup_f32(&szf);

        /* And subtract them out */
        confp = &coords[k * n_atoms * 3];
        for (i = 0; i < n_atoms/4; i++) {
            xyz = aos_deinterleaved_load(confp);
            xyz.val[0] = vsubq_f32(xyz.val[0], mux_);
            xyz.val[1] = vsubq_f32(xyz.val[1], muy_);
            xyz.val[2] = vsubq_f32(xyz.val[2], muz_);

            xyz2.val[0] = vmulq_f32(xyz.val[0], xyz.val[0]);        
            xyz2.val[1] = vmulq_f32(xyz.val[1], xyz.val[1]);
            xyz2.val[2] = vmulq_f32(xyz.val[2], xyz.val[2]);
            trace_ = vaddq_f64(trace_, vcvt_f64_f32(vget_low_f32(xyz2.val[0])));
            trace_ = vaddq_f64(trace_, vcvt_f64_f32(vget_low_f32(xyz2.val[1])));
            trace_ = vaddq_f64(trace_, vcvt_f64_f32(vget_low_f32(xyz2.val[2])));
            trace_ = vaddq_f64(trace_, vcvt_f64_f32(vget_high_f32(xyz2.val[0])));
            trace_ = vaddq_f64(trace_, vcvt_f64_f32(vget_high_f32(xyz2.val[1])));
            trace_ = vaddq_f64(trace_, vcvt_f64_f32(vget_high_f32(xyz2.val[2])));

            aos_interleaved_store(confp, xyz.val[0], xyz.val[1], xyz.val[2]);
            confp += 12;
        }
        trace = vgetq_lane_f64(trace_, 0) + vgetq_lane_f64(trace_, 1);
        
        for (i = 0; i < n_atoms % 4; i++) {
            confp[i*3 + 0] -= sxf;
            confp[i*3 + 1] -= syf;
            confp[i*3 + 2] -= szf;
            trace += confp[i*3 + 0]*confp[i*3 + 0];
            trace += confp[i*3 + 1]*confp[i*3 + 1];
            trace += confp[i*3 + 2]*confp[i*3 + 2];
        }        
        if (traces != NULL)
            traces[k] = (float) trace;
    }
}
