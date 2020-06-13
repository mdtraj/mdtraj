#include "util_sse.h"

void inplace_center_and_trace_atom_major(
    float* coords,
    float* traces,
    const int n_frames,
    const int n_atoms
) {
    /* Center a trajectory containing multiple conformations inplace.
       The coordinates are store in float, but the accumulation is done in
       double.

       Also compute the traces of the centered conformations, which are necessary
       for RMSD.
    */
    long long i, k;
    float* confp;
    __m128d sx_, sy_, sz_, trace_;
    __m128 mux_, muy_, muz_;
    float sxf, syf, szf;
    double sx[2], sy[2], sz[2], trace[2];
    __m128 x, y, z, x2, y2, z2;

    #ifdef _OPENMP
    #pragma omp parallel for shared(coords, traces) \
        private(sx_, sy_, sz_, trace_, mux_, muy_, muz_, sxf, syf, szf, \
        confp, i, x, y, z, x2, y2, z2, sx, sy, sz, trace)
    #endif
    for (k = 0; k < n_frames; k++) {
        confp = &coords[k * n_atoms * 3];
        sx_ = sy_ = sz_ = trace_ = _mm_setzero_pd();
        for (i = 0; i < n_atoms/4; i++) {
            aos_deinterleaved_loadu(confp, &x, &y, &z);

            /* accumulate the sums of each coordinate in double */
            /* get the first two values from each float4 */
            sx_ = _mm_add_pd(sx_, _mm_cvtps_pd(x));
            sy_ = _mm_add_pd(sy_, _mm_cvtps_pd(y));
            sz_ = _mm_add_pd(sz_, _mm_cvtps_pd(z));
            /* and shuffle in the second two values */
            sx_ = _mm_add_pd(sx_, _mm_cvtps_pd(_mm_movehl_ps(x, x)));
            sy_ = _mm_add_pd(sy_, _mm_cvtps_pd(_mm_movehl_ps(y, y)));
            sz_ = _mm_add_pd(sz_, _mm_cvtps_pd(_mm_movehl_ps(z, z)));
            confp += 12;
        }
        /* copy the summed coordinates out of the SSE registers */
        _mm_storeu_pd(sx, sx_);
        _mm_storeu_pd(sy, sy_);
        _mm_storeu_pd(sz, sz_);

        /* Add the last couple entries that weren't a factor of four */
        for (i = 0; i < n_atoms % 4; i++) {
            sx[0] += confp[i*3 + 0];
            sy[0] += confp[i*3 + 1];
            sz[0] += confp[i*3 + 2];
        }

        /* Put everything into the first value. We're doing this here, as */
        /* opposed to using a SSE horizontal add. */
        sx[0] += sx[1];
        sy[0] += sy[1];
        sz[0] += sz[1];

        /* Now we want mean x, y, and z positions */
        sx[0] /= n_atoms;
        sy[0] /= n_atoms;
        sz[0] /= n_atoms;

        /* Load these mean positions back into the SSE registers */
        sxf = (float) sx[0];
        syf = (float) sy[0];
        szf = (float) sz[0];
        mux_ = _mm_load1_ps(&sxf);
        muy_ = _mm_load1_ps(&syf);
        muz_ = _mm_load1_ps(&szf);

        /* And subtract them out */
        confp = &coords[k * n_atoms * 3];
        for (i = 0; i < n_atoms/4; i++) {
            aos_deinterleaved_loadu(confp, &x, &y, &z);
            x = _mm_sub_ps(x, mux_);
            y = _mm_sub_ps(y, muy_);
            z = _mm_sub_ps(z, muz_);

            x2 = _mm_mul_ps(x, x);
            y2 = _mm_mul_ps(y, y);
            z2 = _mm_mul_ps(z, z);
            trace_ = _mm_add_pd(trace_, _mm_cvtps_pd(x2));
            trace_ = _mm_add_pd(trace_, _mm_cvtps_pd(y2));
            trace_ = _mm_add_pd(trace_, _mm_cvtps_pd(z2));
            trace_ = _mm_add_pd(trace_, _mm_cvtps_pd(_mm_movehl_ps(x2, x2)));
            trace_ = _mm_add_pd(trace_, _mm_cvtps_pd(_mm_movehl_ps(y2, y2)));
            trace_ = _mm_add_pd(trace_, _mm_cvtps_pd(_mm_movehl_ps(z2, z2)));

            aos_interleaved_storeu(confp, x, y, z);
            confp += 12;
        }
        _mm_storeu_pd(trace, trace_);

        for (i = 0; i < n_atoms % 4; i++) {
            confp[i*3 + 0] -= sxf;
            confp[i*3 + 1] -= syf;
            confp[i*3 + 2] -= szf;
            trace[0] += confp[i*3 + 0]*confp[i*3 + 0];
            trace[0] += confp[i*3 + 1]*confp[i*3 + 1];
            trace[0] += confp[i*3 + 2]*confp[i*3 + 2];
        }
        trace[0] += trace[1];
        if (traces != NULL)
            traces[k] = (float) trace[0];
    }
}
