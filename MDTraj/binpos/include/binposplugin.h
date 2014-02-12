#ifndef BINPOS_H
#define BINPOS_H
#include "molfile_plugin.h"

void *open_binpos_read(const char *path, const char *filetype, int *natoms);
int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts);
void close_file_read(void *v);

void *open_binpos_write(const char *path, const char *filetype, int natoms);
int write_timestep(void *v, const molfile_timestep_t *ts);
void close_file_write(void *v);

int seek_timestep(void* v, long int offset, int origin);
long int tell_timestep(void* v);
#endif
