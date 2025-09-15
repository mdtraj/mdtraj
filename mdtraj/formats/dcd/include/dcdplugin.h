#include "fastio.h"       /* must come before others, for O_DIRECT...   */
#include "molfile_plugin.h"

#ifndef DCDPLUGIN_H
#define DCDPLUGIN_H

typedef struct {
  fio_fd fd;
  int natoms;
  int nsets;
  int setsread;
  int istart;
  int nsavc;
  double delta;
  int nfixed;
  float *x, *y, *z;
  int *freeind;
  float *fixedcoords;
  int reverse;
  int charmm;
  int first;
  int with_unitcell;
} dcdhandle;

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define RECSCALE32BIT 1
#define RECSCALE64BIT 2
#define RECSCALEMAX   2

/* Define error codes that may be returned by the DCD routines */
#define DCD_SUCCESS      0  /* No problems                     */
#define DCD_EOF         -1  /* Normal EOF                      */
#define DCD_DNE         -2  /* DCD file does not exist         */
#define DCD_OPENFAILED  -3  /* Open of DCD file failed         */
#define DCD_BADREAD     -4  /* read call on DCD file failed    */
#define DCD_BADEOF      -5  /* premature EOF found in DCD file */
#define DCD_BADFORMAT   -6  /* format of DCD file is wrong     */
#define DCD_FILEEXISTS  -7  /* output file already exists      */
#define DCD_BADMALLOC   -8  /* malloc failed                   */
#define DCD_BADWRITE    -9  /* write call on DCD file failed   */

/* Define feature flags for this DCD file */
#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

/* defines used by write_dcdstep */
#define NFILE_POS 8L
#define NSTEP_POS 20L

/* READ Macro to make porting easier */
#define READ(fd, buf, size)  fio_fread(((void *) buf), (size), 1, (fd))

/* WRITE Macro to make porting easier */
#define WRITE(fd, buf, size) fio_fwrite(((void *) buf), (size), 1, (fd))

/* XXX This is broken - fread never returns -1 */
#define CHECK_FREAD(X, msg) if (X==-1) { return(DCD_BADREAD); }
#define CHECK_FEOF(X, msg)  if (X==0)  { return(DCD_BADEOF); }

dcdhandle* open_dcd_read(const char *path, const char *filetype, int *natoms, int* nsets);
void close_file_read(dcdhandle *v);
int read_next_timestep(dcdhandle *v, int natoms, molfile_timestep_t *ts);

dcdhandle* open_dcd_write(const char *path, const char *filetype, const int natoms,
                          const int with_unitcell);
int write_timestep(dcdhandle *v, const molfile_timestep_t *ts);
void close_file_write(dcdhandle *v);
int dcd_nsets(dcdhandle* v);
int dcd_rewind(dcdhandle* dcd);

#endif
