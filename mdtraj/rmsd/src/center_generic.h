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
    float mux_, muy_, muz_;

    float* conf;

    for (k = 0; k < n_frames; k++) {
        double sx = 0, sy = 0, sz = 0, trace = 0;
        conf = &coords[k*n_atoms*3];
        for (i = 0; i < n_atoms; i++) {
            sx += conf[i*3];
            sy += conf[i*3+1];
            sz += conf[i*3+2];
        }

        mux_ = (float) (sx / n_atoms);
        muy_ = (float) (sy / n_atoms);
        muz_ = (float) (sz / n_atoms);



        conf = &coords[k * n_atoms * 3];
        for (i = 0; i < n_atoms; i++) {
            conf[i*3] -= mux_;
            conf[i*3 + 1] -= muy_;
            conf[i*3 + 2] -= muz_;
            trace += conf[i*3]*conf[i*3] + conf[i*3 + 1]*conf[i*3 + 1] + conf[i*3 + 2]*conf[i*3 + 2];

        }
        if (traces != NULL)
            traces[k] = (float) trace;
    }
}
