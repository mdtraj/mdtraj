#include <xtc.hh>
#include <iostream>
#include <typeinfo>
#include <assert.h>
#include <malloc.h>

extern "C" {
  #include "xdrfile_xtc.h"
  #include "xdrfile.h"
}

using namespace std;

XTCReader::XTCReader(char* filename_, int firstframe_, int lastframe_, int stride_) {
  filename = filename_;
  firstframe = firstframe_;
  lastframe = lastframe_;
  stride = stride_;
  step = firstframe_;

  read_xtc_natoms(filename, & n_atoms);
  xdrfile = xdrfile_open(filename, "r");

  cout << n_atoms << endl;
  cout << filename << endl;
}
void XTCReader::set_outarray(float* xyz_, int n_frames_, int n_atoms_, int n_dims_) throw (TypeError) {
  if (n_atoms_ != n_atoms) {
    cout << "natoms wrong" << endl;
    cout << "you gave " << n_atoms_ << "should be " << n_atoms << endl;
    throw TypeError();
  }
  if (n_dims_ != 3) {
    cout << "dims wrong" << endl;
    throw TypeError();
  }

  xyz = xyz_;
  n_frames = n_frames_;
  n_atoms = n_atoms_;
  n_dims = n_dims_;
}


int XTCReader::next() throw (StopIteration) {
  float time = 0;
  float (*box)[3] = (float(*)[3]) malloc(9*sizeof(float));
  float (*coords)[3] = (float(*)[3]) xyz;
  float precision = 0;
  int result;
  
  result = read_xtc(xdrfile, 22, &step, &time, box, coords, &precision);

  cout << "time " << time << endl;
  cout << "prec " << precision << endl;
  cout << "step " << step << endl;
  //xdrfile_close(xdrfile);
  
  free(box);
}

XTCReader XTCReader::__iter__() {
  return *this;
}
