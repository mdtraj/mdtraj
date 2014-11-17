#ifndef MIXTAPE_RMSD_CENTER_H
#define MIXTAPE_RMSD_CENTER_H

void inplace_center_and_trace_atom_major(float* coords, float* traces,
    const int n_frames, const int n_atoms);

#endif