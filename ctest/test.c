#include <stdlib.h>
#include <stdio.h>

#ifdef LOCAL
int write_xtc(void* xd, int natoms, int step, float time, float* box, float* x, float prec)
/* Write a frame to xtc file */
{
  int result;
  
  //robert hack
  printf("write_xtc: step %d, time %f, prec %f \n", step, time, prec);
  
  /*if ((result = xtc_header(xd,&natoms,&step,&time,FALSE)) != exdrOK)
    return result;

  if ((result = xtc_coord(xd,&natoms,box,x,&prec,0)) != exdrOK)
    return result;
  
  return exdrOK;*/
  return 1;
}
#endif


int main(int argc, char* argv[]) {
  long fh;


  int natoms = 10;
  int step=11;
  float time=12.0;
  float prec=13.0;
  float* box = (float*) malloc(sizeof(float)*9);
  float* x = (float*) malloc(sizeof(float)*natoms*3);
  
  printf("main: step %d, time %f, prec %f \n", step, time, prec);
   
  write_xtc(&fh, natoms, step, time, box, x, prec);

  return 1;
}

