extern "C" {
  #include "xdrfile_xtc.h"
  #include "xdrfile.h"
}

class StopIteration {};
class TypeError {};

class XTCReader {
  char* filename;
  int firstframe;
  int lastframe;
  int stride;

  int step;

  XDRFILE* xdrfile;

  //the array that we're going to put the output into
  float* xyz;
  int n_frames;
  int n_atoms;
  int n_dims;
  
public:
  XTCReader(char* filename_, int firstframe_, int lastframe_, int stride_);
  void set_outarray(float* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) throw (TypeError);
  //get the next frame
  int next() throw (StopIteration);

  //required for the iterator interface
  XTCReader __iter__();
};

