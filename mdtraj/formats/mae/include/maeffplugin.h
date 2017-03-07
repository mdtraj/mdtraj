#include <molfile_plugin.h>

void *open_file_read(const char *fname, const char *ftype, int *vmdatoms);
int read_structure(void* v, int* optflags, molfile_atom_t *atoms);
int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts);
void close_file_read( void *v);
