//
// Version info for VMD plugin tree:
//   $Id: dtrplugin.cxx,v 1.22 2011/12/23 22:40:52 johns Exp $
//
// Version info for last sync with D. E. Shaw Research:
//  //depot/desrad/main/sw/libs/molfile/plugins/dtrplugin.cxx#30
//

/*
Copyright 2009, D. E. Shaw Research, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research, LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "dtrplugin.hxx"

using namespace desres::molfile;

#include <sstream>
#include <ios>
#include <iomanip>
#include <math.h>
#include <errno.h>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>
#include <ios>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "vmddir.h"
#include "endianswap.h"


static const char SERIALIZED_VERSION[] = "0006";
const char * desres::molfile::dtr_serialized_version() {
  return SERIALIZED_VERSION;
}

static bool badversion(const std::string& version) {
    return version != SERIALIZED_VERSION;
}

#ifndef DESRES_WIN32
static const char s_sep = '/';
#include <sys/mman.h>

#include <netinet/in.h> /* for htonl */
#if defined(_AIX)
#include <fcntl.h>
#else
#include <sys/fcntl.h>
#endif

#if defined(__sun)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#define PathIsRelative(x) (x[0] != '/')

#else
/// windows version

#define M_PI (3.1415926535897932385)
#define M_PI_2 (1.5707963267948966192)

#ifndef S_ISREG
#define S_ISREG(x) (((x) & S_IFMT) == S_IFREG)
#endif

#ifndef S_ISDIR
#define S_ISDIR(x) (((x) & S_IFMT) == S_IFDIR)
#endif

static const char s_sep = '\\';

#include "Shlwapi.h"

#endif

static const uint32_t magic_timekey = 0x4445534b;
static const uint32_t magic_frame   = 0x4445534d;
static const uint32_t s_version     = 0x00000100;
static const uint32_t s_irosetta    = 0x12345678;
static const float    s_frosetta    = 1234.5;
static const double   s_drosetta    = 1234.5e6;
static const uint32_t s_lrosetta_lo = 0x89abcdef;
static const uint32_t s_lrosetta_hi = 0x01234567;
static const uint32_t s_blocksize   = 4096;
static const uint32_t s_alignsize   = 8;

namespace {

  const double PEAKmassInAmu = 418.4;

  double sfxp_ulp32flt(int32_t x) {
    return ldexp(((double) x),-31);
  }


  uint64_t assemble64( uint32_t lo, uint32_t hi) {
    uint64_t hi64 = hi; 
    return (hi64 << 32) | lo; 
  }
  
  double assembleDouble(uint32_t lo, uint32_t hi) {
    union {
      uint64_t ival;
      double   dval;
    } u;
    u.ival = assemble64(lo,hi);
    return u.dval;
  }

  /// definitions of binary representation of frameset files

  typedef struct {
    uint32_t magic;           //!< Magic number
    uint32_t version;         //!< Version of creator
    uint32_t framesize_lo;    //!< bytes in frame (low)
    uint32_t framesize_hi;    //!< bytes in frame (high)
  } required_header_t;

  //! Header structure within file.
  typedef struct {
    required_header_t required; //!< 4 word mini-header

    uint32_t size_header_block; //!< Size of this header
    uint32_t unused0;         //!< not used in current implementation
    uint32_t irosetta;        //!< 32-bit integer rosetta value.
    float    frosetta;        //!< 32-bit float rosetta

    uint32_t drosetta_lo;     //!< 64-bit float rosetta (low)
    uint32_t drosetta_hi;     //!< 64-bit float rosetta (high)
    uint32_t lrosetta_lo;     //!< 64-bit integer rosetta (low)
    uint32_t lrosetta_hi;     //!< 64-bit integer rosetta (high)

    uint32_t endianism;       //!< Endianism of writer machine.
    uint32_t nlabels;         //!< Number of labeled fields.
    uint32_t size_meta_block; //!< Number of bytes of meta information (padded)
    uint32_t size_typename_block; //!< Number of bytes of typenames (padded)

    uint32_t size_label_block; //!< Number of bytes to store label strings (padded)
    uint32_t size_scalar_block; //!< Number of bytes of scalar storage (padded)
    uint32_t size_field_block_lo; //!< Number of bytes of field storage (padded)
    uint32_t size_field_block_hi; //!< Number of bytes of field storage (padded)

    uint32_t size_crc_block;  //!< Size of the trailing CRC field (unused!)
    uint32_t size_padding_block; //!< Number of ignored bytes to pad to pagesize boundary.
    uint32_t unused1;         //!< Not used in current implementation.
    uint32_t unused2;         //!< Not used in current implementation.

  } header_t;

  typedef struct {
    uint32_t type;            //!< \brief Typecode for this type.
    uint32_t elementsize;     //!< \brief Number of bytes in an element
    uint32_t count_lo;        //!< \brief Number of elements (low)
    uint32_t count_hi;        //!< \brief Number of elements (high)
  } metadisk_t;

  typedef struct key_prologue {
    uint32_t magic;           /* Magic number for frames */
    uint32_t frames_per_file; /* Number of frames in each file */
    uint32_t key_record_size; /* The size of each key record */
  } key_prologue_t;


  //// utility routines

  /*!
   * Extracts the low 32 bits of a 64 bit integer by masking.
   */
  uint32_t lobytes(const uint64_t& x) {
    uint32_t mask = 0xffffffff;
    return x & mask;
  }

  /*!
   * Extract the high 32 bits of a 64 bit integer by shifting.
   */
  uint32_t hibytes(const uint64_t& x) {
    return x >> 32;
  }

  /*!
   * Extract the low 32 bits of a 64 bit float as an integer.
   */
  uint32_t lobytes(const double& x) {
    union {
      uint64_t ival;
      double   dval;
    } u;
    u.dval = x;
    return lobytes(u.ival);
  }

  /*!
   * Extract the high 32 bits of a 64 bit float as an integer.
   */
  uint32_t hibytes(const double& x) {
    union {
      uint64_t ival;
      double   dval;
    } u;
    u.dval = x;
    return hibytes(u.ival);
  }

  /*!
   * The byte order associated with this machine.  We use
   * 1234 for little endian, 4321 for big endian, and
   * 3412 for the unlikely PDB endianism.
   */
  uint32_t machineEndianism() {
#if __BYTE_ORDER == __LITTLE_ENDIAN
    uint32_t byteorder = 1234;
#else
#if __BYTE_ORDER == __BIG_ENDIAN
    uint32_t byteorder = 4321;
#else
#ifdef PDB_ENDIAN
#if __BYTE_ORDER == __PDB_ENDIAN
    uint32_t byteorder = 3412;
#endif
#endif
#endif
#endif
    // If we get a compile error here, then __BYTE_ORDER
    // has an unexpected value.
    return byteorder;
  }

  uint64_t alignInteger( const uint64_t &x, unsigned border) {
    return x + (border - x%border)%border;
  }


  /*!
   * See RFC 1146 for Fletcher's Checksum (http://tools.ietf.org/html/rfc1146)
   */
  uint32_t fletcher( const uint16_t *data, size_t len ) {
    uint32_t sum1 = 0xffff, sum2 = 0xffff;
 
    while (len) {
      unsigned tlen = len > 360 ? 360 : len;
      len -= tlen;
      do {
        sum1 += *data++;
        sum2 += sum1;
      } while (--tlen);
      sum1 = (sum1 & 0xffff) + (sum1 >> 16);
      sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    }
    /* Second reduction step to reduce sums to 16 bits */
    sum1 = (sum1 & 0xffff) + (sum1 >> 16);
    sum2 = (sum2 & 0xffff) + (sum2 >> 16);
    return sum2 << 16 | sum1;
  }

  bool isfile(const std::string &name) {
    struct stat statbuf;
    return (stat(name.c_str(),&statbuf) == 0 && S_ISREG(statbuf.st_mode));
  }

  /*!
   * Remove a file or directory.  For directories,
   * we recurse through subfiles and remove those
   * before attempting the ::rmdir();
   */
  void recursivelyRemove(std::string path) {
    struct stat statbuf;

    // -----------------------------------------------
    // Only try to unlink if the file exists
    // We recurse through directories and unlink
    // other files.
    // -----------------------------------------------

#ifdef DESRES_WIN32
    // Use ::stat instead of ::lstat on windows since there are no symlinks
    if (stat(path.c_str(),&statbuf) == 0) {
#else
    if (::lstat(path.c_str(),&statbuf) == 0) {
#endif
      if (!S_ISDIR(statbuf.st_mode)) {
        if (::unlink(path.c_str()) != 0) {
            throw std::runtime_error(strerror(errno));
        }
      } else {
        VMDDIR* directory = NULL;
        try {
          directory = vmd_opendir(path.c_str());
          if (directory) {
            // Remove subfiles
            char * entry;
            while( (entry=vmd_readdir(directory)) != NULL ) {
              // Don't unlink . or ..
              if (entry[0] == '.') {
                if (entry[1] == 0) continue;
                if (entry[1] == '.' && entry[2] == 0) continue;
              }
              recursivelyRemove(path + s_sep + entry);
            }
            vmd_closedir(directory);
            directory = NULL;

            // Remove the actual directory
            if (::rmdir(path.c_str()) != 0) {
              throw std::runtime_error(strerror(errno));
            }
          }
        } catch(...) {
          if (directory) vmd_closedir(directory);
          throw;
        }
      }
    }
  }


  ////////
  // CRC
  ////////
  
  typedef uint32_t crc;
  
  #define POLYNOMIAL 0x04C11DB7
  #define WIDTH  (8 * sizeof(crc))
  #define TOPBIT (1 << (WIDTH - 1))
  #define FINAL_XOR_VALUE 0xFFFFFFFF

  crc processByte( crc remainder, char msg ) {
        remainder ^= (msg << (WIDTH - 8));
        for (uint8_t bit = 8; bit > 0; --bit)
        {
            if (remainder & TOPBIT) {
                remainder = (remainder << 1) ^ POLYNOMIAL;
            } else {
                remainder = (remainder << 1);
            }
        }
        return remainder;
  }

  crc processBytes(const char *message, int nBytes) {
    crc  remainder = 0;	
    for (int byte = 0; byte < nBytes; ++byte) {
        remainder = processByte( remainder, message[byte] );
    }
    return remainder;
  } 

  int32_t cksum(const std::string &s) {
    size_t len = s.size();
    int32_t result = processBytes( s.c_str(), len );
  
    for ( ; len; len >>= 8) {
      result = processByte( result, len & 0xff );
    }
    return result ^ FINAL_XOR_VALUE;
  }

}

bool Timekeys::init(const std::string& path ) {
    std::string timekeys_path = path;
    timekeys_path += s_sep;
    timekeys_path += "timekeys";
    FILE * fd = fopen( timekeys_path.c_str(), "rb" );
    if (!fd) {
      fprintf(stderr, "Could not find timekeys file at %s\n", 
          timekeys_path.c_str());
      return false;
    }
  
    /* check the magic number */
    key_prologue_t prologue[1];
    if (fread( prologue, sizeof(key_prologue_t), 1, fd )!=1) {
      fprintf(stderr, "Failed to read key prologue from %s\n",
          timekeys_path.c_str());
      fclose(fd);
      return false;
    }
    prologue->magic = htonl(prologue->magic);
    if (prologue->magic != magic_timekey) {
      fprintf(stderr, "timekeys magic number %x doesn't match %x\n",
          prologue->magic, magic_timekey);
      fclose(fd);
      return false;
    }
  
    /* get frames per file and key record size */
    prologue->frames_per_file = ntohl( prologue->frames_per_file );
    prologue->key_record_size = ntohl( prologue->key_record_size );
    m_fpf = prologue->frames_per_file;
  
    /* read all key records */
    fseek(fd, 0, SEEK_END);
    off_t keyfile_size = ftello(fd);
    size_t nframes = (keyfile_size-sizeof(key_prologue_t))/sizeof(key_record_t);
  
    keys.resize(nframes);
    fseek(fd, sizeof(key_prologue_t), SEEK_SET);
    if (fread(&keys[0], sizeof(key_record_t), nframes, fd)!=nframes) {
      fprintf(stderr, "Failed to read all timekeys records: %s\n",
          strerror(errno));
      fclose(fd);
      return false;
    }
    fclose(fd);
  
    /* Check that we didn't get zero-length frames; this would be a strong
     * indicator of file corruption! */
    int warning_count=0;
    size_t i;
    for (i=0; i<nframes; i++) {
        if (keys[i].size()==0) {
            if (++warning_count<10) {
              fprintf(stderr, "dtrplugin -- WARNING: timekey %d of dtr %s reports 0-length frame; file corruption likely.\n", (int)i, path.c_str());
            }
            if (warning_count==10) {
                fprintf(stderr, "dtrplugin -- WARNING: skipping remaining warnings in dtr %s\n", path.c_str());
            }
        }
    }
    if (warning_count) {
        fprintf(stderr, "dtrplugin -- WARNING: found %d likely corrupt timekeys in %s\n",
                warning_count, path.c_str());
    }

    m_size = m_fullsize = keys.size();
    if (!keys.size()) return true;

    m_first = keys[0].time();
    m_framesize = keys[0].size();
    if (keys.size()==1) {
        m_interval=0;
        keys.clear();
        return true;
    }
    m_interval=keys[1].time()-keys[0].time();
    for (i=1; i<keys.size(); i++) {
        if (keys[i].size() == 0) {
            /* ignore obviously corrupt frames */
            continue;
        }
        /* constant frame size */
        if (keys[i].size() != m_framesize) {
            fprintf(stderr, "non-constant framesize at frame %ld\n", i);
            printf("size %d framesize %d\n\n",
                    keys[i].size(), m_framesize);
            return true;
        }
        /* constant time interval */
        if (fabs((keys[i].time()-keys[i-1].time())-m_interval) > 1e-3) {
            if (getenv("DTRPLUGIN_VERBOSE")) fprintf(stderr, 
                    "non-constant time interval at frame %ld\n", i);
            return true;
        }
        /* constant offset */
        if (keys[i].offset() != m_framesize*( i % m_fpf)) {
            fprintf(stderr, "unexpected offset for frame %ld\n", i);
            return true;
        }
    }
    /* looks good!  Don't need the explicit key records anymore */
    keys.clear();
    return true;
}

key_record_t Timekeys::operator[](uint64_t i) const {
    if (i>m_fullsize) throw std::runtime_error("frame index out of range");
    if (keys.size()) return keys.at(i);

    key_record_t timekey;
#if defined(_MSC_VER)
    double time = m_first + ((__int64) i)*m_interval;
#else
    double time = m_first + i*m_interval;
#endif
    uint64_t offset = (i % m_fpf) * m_framesize;

    timekey.time_lo = htonl(lobytes(time));
    timekey.time_hi = htonl(hibytes(time));
    timekey.offset_lo = htonl(lobytes(offset));
    timekey.offset_hi = htonl(hibytes(offset));
    timekey.framesize_lo = htonl(lobytes(m_framesize));
    timekey.framesize_hi = htonl(hibytes(m_framesize));
    return timekey;
}

namespace {
    template <typename T> 
    void rawdump(std::ostream& out, const T& v) {
        out.write((char *)&v, sizeof(v));
    }

    template <typename T> 
    void rawload(std::istream& in, T& v) {
        in.read((char *)&v, sizeof(v));
    }
}

void Timekeys::dump(std::ostream& out) const {
    rawdump(out, m_first);
    rawdump(out, m_interval);
    rawdump(out, m_framesize);
    rawdump(out, m_size);
    rawdump(out, m_fullsize);
    rawdump(out, m_fpf);
    rawdump(out, keys.size());
    if (keys.size()) {
        out.write((const char *)&keys[0], keys.size()*sizeof(keys[0]));
    }
}

void Timekeys::load(std::istream& in) {
    size_t sz;
    rawload(in, m_first);
    rawload(in, m_interval);
    rawload(in, m_framesize);
    rawload(in, m_size);
    rawload(in, m_fullsize);
    rawload(in, m_fpf);
    rawload(in, sz);
    if (sz) {
        keys.resize(sz);
        in.read((char *)&keys[0], keys.size()*sizeof(keys[0]));
    }
}

namespace {
  struct Blob {
    std::string type;
    uint64_t count;
    const void *data;
    bool byteswap;

    Blob() : count(0), data(0) {}
    Blob( const std::string &type_, uint64_t count_, const void *data_,
          uint32_t frame_endianism )
    : type(type_), count(count_), data(data_), byteswap(false) {
      uint32_t my_endianism = machineEndianism();
      if (frame_endianism != my_endianism) {
        if ( (frame_endianism==1234 && my_endianism==4321) || 
             (frame_endianism==4321 && my_endianism==1234) ) {
          byteswap=true;
        } else {
          throw std::runtime_error("Unable to handle frame endianness");
        }
      }
    }

    std::string str() const {
      if (type=="char" && count>0) {
        const char *s=(const char *)data;
        return std::string(s, s+count);
      }
      return "";
    }
    void get_float(float *buf) const {
      if (type=="float") {
        memcpy(buf, data, count*sizeof(float));
      } else if (type=="double") {
        const double *p = reinterpret_cast<const double *>(data);
        for (uint64_t i=0; i<count; i++) buf[i] = p[i];
      } else {
        memset(buf, 0, count*sizeof(float));
      }
      if (byteswap) swap4_unaligned(buf, count);
    }
    void get_double(double *buf) const {
      if (type=="double") {
        memcpy(buf, data, count*sizeof(double));
      } else if (type=="float") {
        const float *p = reinterpret_cast<const float *>(data);
        for (uint64_t i=0; i<count; i++) buf[i] = p[i];
      } else {
        memset(buf, 0, count*sizeof(double));
      }
      if (byteswap) swap8_unaligned(buf, count);
    }
    void get_int32(int32_t *buf) const {
      if (type=="int32_t") {
        memcpy(buf, data, count*sizeof(int32_t));
      } else {
        memset(buf, 0, count*sizeof(int32_t));
      }
      if (byteswap) swap4_unaligned(buf, count);
    }
    void get_uint32(uint32_t *buf) const {
      if (type=="uint32_t") {
        memcpy(buf, data, count*sizeof(uint32_t));
      } else {
        memset(buf, 0, count*sizeof(uint32_t));
      }
      if (byteswap) swap4_unaligned(buf, count);
    }
  };

  typedef std::map<std::string, Blob> BlobMap;
}

static inline std::string addslash(const std::string& s){
    return (s.rbegin()[0] == '/') ? s : s + "/";
}

#define DD_RELPATH_MAXLEN (9) 
static std::string 
DDreldir(const std::string& fname, int ndir1, int ndir2){

    if( fname.find('/', 0) != std::string::npos ) {
      fprintf(stderr, "DDreldir: filename '%s' must not contain '/'\n",
          fname.c_str());
      return "";
    }

    uint32_t hash = cksum(fname);

    // uint32_t u1 = ndir1;
    // uint32_t u2 = ndir2;
    uint32_t d1, d2;
    char answer[DD_RELPATH_MAXLEN];
    if(ndir1 > 0){
	d1 = hash%ndir1;
	if(ndir2 > 0){
	    d2 = (hash/ndir1)%ndir2;
	    sprintf(answer, "%03x/%03x/", d1, d2);
	}else{
	    sprintf(answer, "%03x/", d1);
	}
    }else{
	sprintf(answer, "./");
    }
    return std::string(answer);
}

namespace {
  class DDException : public std::runtime_error{
  public:
      int eno;
      DDException(const std::string &text, int _eno=0) 
      : std::runtime_error(text + strerror(eno)), eno(_eno){}
  };
}

void DDmkdir(const std::string &dirpath, mode_t mode, int ndir1, int ndir2){
    std::string dpslash(addslash(dirpath));

    mode_t openmode = mode | 0300; // make sure we can write into the directory
    if( mkdir(dpslash.c_str(), openmode) < 0 )
	throw DDException("mkdir", errno);
	
    if( mkdir((dpslash + "not_hashed").c_str(), openmode) < 0 )
        throw DDException("mkdir not_hashed subdirectory", errno);

    FILE *fp = fopen((dpslash + "not_hashed/.ddparams").c_str(), "w");
    if(fp == NULL)
	throw DDException("fopen( .ddparams, \"w\" )", errno);
    if( fprintf(fp, "%d %d\n", ndir1, ndir2) < 0 ){
        fclose(fp);
	throw DDException("fprintf(.ddparams ...)", errno);
    }
    if( fclose(fp) )
	throw DDException("fclose(.ddparams)", errno);

    for(int i=0; i<ndir1; ++i){
	char sub[6];
	sprintf(sub, "%03x/", i);
	std::string dirsub = dpslash + sub;
        {
	    if( mkdir(dirsub.c_str(), openmode) < 0 )
		throw DDException("mkdir " + dirsub, errno);
	}
	for(int j=0; j<ndir2; ++j){
	    char subsub[6];
	    sprintf(subsub, "%03x", j);
	    std::string dirsubsub = dirsub + subsub;
	    if( mkdir(dirsubsub.c_str(), mode) < 0 ) // NOT openmode!
		throw DDException("mkdir " + dirsubsub, errno);
	}
        if( mode != openmode ){
            // change the mode back to what the user requested now
            // that we're done creating stuff...
            if( chmod(dirsub.c_str(), mode) < 0 )
                throw DDException("chmod " + dirsub, errno);
        }
    }
    if( mode != openmode ){
        // change the mode back to what the user requested now
        // that we're done creating stuff...
        if( chmod(dpslash.c_str(), mode) < 0 )
            throw DDException("chmod " + dpslash, errno);
        if( chmod((dpslash + "not_hashed").c_str(), mode) < 0 )
          throw DDException("chmod " + dpslash + "not_hashed", errno);
    }
}


static void 
DDgetparams(const std::string& dirpath, int *ndir1, int *ndir2) {
  // get ddparams, or assume (0,0) and let the frame file opens fail.
  *ndir1 = *ndir2 = 0;
  std::string dirslash(addslash(dirpath));
  // New convention - .ddparams is in not_hashed/.
  FILE *fp = fopen((dirslash + "not_hashed/.ddparams").c_str(), "r");
  // Allow the old convention of placing .ddparams in the top-level.
  if( fp == NULL && errno == ENOENT ) {
      fp = fopen((dirslash + ".ddparams").c_str(), "r");
  }
  if(fp != NULL) {
    if( fscanf(fp, "%d%d", ndir1, ndir2) != 2 ) 
      fprintf(stderr, "Failed to parse .ddparams; assuming flat structure\n");
    if( fclose(fp) ) {
      fprintf(stderr, "Warning: Failed to close .ddparams file: %s\n",
          strerror(errno));
    }
  }
}

static std::string framefile( const std::string &dtr,
                              size_t frameno, 
                              size_t frames_per_file,
                              int ndir1,
                              int ndir2) {
  unsigned frame_file = frameno / frames_per_file;
  std::ostringstream filename;
  filename << "frame" << std::setfill('0') << std::setw(9)
           << frame_file;
  std::string fname = filename.str();

  std::string fullpath(dtr);
  fullpath += "/";
  fullpath += DDreldir(fname, ndir1, ndir2);
  fullpath += fname;
  return fullpath;
}

static BlobMap read_frame( const void *mapping, uint64_t len ) {

    const char *base = reinterpret_cast<const char *>(mapping);
    const header_t *header = reinterpret_cast<const header_t*>(base);
    if (len<sizeof(header_t)) 
        throw std::runtime_error("Frame size is smaller than header_t");
    if (ntohl(header->required.magic) != magic_frame) {
        char buf[256];
        sprintf(buf, "invalid magic number: expected %d, got %d\n",
                magic_frame, ntohl(header->required.magic));
        throw std::runtime_error(buf);
    }

    uint32_t size_header_block = ntohl(header->size_header_block);
    uint32_t frames_endianism = ntohl(header->endianism);
    uint32_t frames_nlabels = ntohl(header->nlabels);
    uint32_t size_meta_block = ntohl(header->size_meta_block);
    uint32_t size_typename_block = ntohl(header->size_typename_block);
    uint32_t size_label_block = ntohl(header->size_label_block);
    uint32_t size_scalar_block = ntohl(header->size_scalar_block);
    uint32_t size_field_block_lo = ntohl(header->size_field_block_lo);
    uint32_t size_field_block_hi = ntohl(header->size_field_block_hi);
    uint64_t size_field_block = assemble64(size_field_block_lo,
                                           size_field_block_hi);

    uint64_t offset_header_block = 0;
    uint64_t offset_meta_block = offset_header_block + size_header_block;
    uint64_t offset_typename_block = offset_meta_block + size_meta_block;
    uint64_t offset_label_block = offset_typename_block + size_typename_block;
    uint64_t offset_scalar_block = offset_label_block + size_label_block;
    uint64_t offset_field_block = offset_scalar_block + size_scalar_block;
    uint64_t offset_crc_block = offset_field_block + size_field_block;

    const metadisk_t* diskmeta  = reinterpret_cast<const metadisk_t*>(base+offset_meta_block);
    const char* typenames = reinterpret_cast<const char*>(base+offset_typename_block);
    const char* labels    = reinterpret_cast<const char*>(base+offset_label_block); 
    const char* scalars   = reinterpret_cast<const char*>(base+offset_scalar_block);
    const char* fields    = reinterpret_cast<const char*>(base+offset_field_block);
    const uint32_t  * crc = reinterpret_cast<const uint32_t*>(base+offset_crc_block);
    if (*crc != 0) {
        uint32_t frame_crc = fletcher(reinterpret_cast<const uint16_t*>(base),offset_crc_block/2);
        if (frame_crc != *crc) {
            throw std::runtime_error("Checksum did not match");
        }
    }
    /* More sanity checks */
    if (len<offset_meta_block+size_meta_block)
        throw std::runtime_error("Frame size cannot contain meta block");
    if (len<offset_typename_block+size_typename_block)
        throw std::runtime_error("F size cannot contain meta block");
    if (len<offset_label_block+size_label_block)
        throw std::runtime_error("F size cannot contain meta block");
    if (len<offset_scalar_block+size_scalar_block)
        throw std::runtime_error("F size cannot contain meta block");
    if (len<offset_field_block+size_field_block)
        throw std::runtime_error("Frame size cannot contain meta block");

    std::vector<std::string> types;
    while(*typenames) {
      if (typenames >= labels) {
        fprintf(stderr, "More typenames than labels!\n");
        break;
      }
      std::string type(typenames);
      types.push_back(type);
      typenames += type.size()+1;
    }

    BlobMap blobs;

    for (size_t ii=0; ii<frames_nlabels; ++ii) {
      std::string label(labels);
      labels += label.size()+1;
      // Pull out the typecode, elementsize, and count
      uint32_t code = ntohl(diskmeta[ii].type);
      uint32_t elementsize = ntohl(diskmeta[ii].elementsize);
      uint32_t count_lo = ntohl(diskmeta[ii].count_lo);
      uint32_t count_hi = ntohl(diskmeta[ii].count_hi);
      uint64_t count = assemble64(count_lo,count_hi);
      uint64_t nbytes = elementsize*count;

      const char *addr=0;
      if (count <= 1) {
        addr=scalars;
        scalars += alignInteger(nbytes, s_alignsize);
      } else {
        addr=fields;
        fields += alignInteger(nbytes, s_alignsize);
      }
      try {
        blobs[label] = Blob( types.at(code), count, addr, frames_endianism );
      }
      catch (std::exception &e) {
        fprintf(stderr, "Failed fetching '%s' data from frame\n", 
            label.c_str());
      }
    }
    return blobs;
}

static void *read_file( int fd, off_t offset, size_t &framesize ) {
  if (fd<=0) {
    fprintf(stderr, "read_file: bad file descriptor\n");
    return NULL;
  }
  if (framesize==0) {
    struct stat statbuf;
    if (fstat(fd,&statbuf)!=0) {
      fprintf(stderr, "Could not stat file: %s\n", strerror(errno));
      return NULL;
    }
    framesize=statbuf.st_size-offset;
  }

  void *mapping = malloc(framesize);
  if (lseek(fd, offset, SEEK_SET)!=offset) {
      fprintf(stderr, "seek to specified offset failed: %s\n", strerror(errno));
      free(mapping);
      return NULL;
  }

  size_t rc = read(fd, mapping, framesize);
  if (rc==0) {
      free(mapping);
      return NULL;
  }
  if (rc==-1) {
      fprintf(stderr, "reading bytes from frame failed: %s\n", strerror(errno));
      free(mapping);
      return NULL;
  }
  if (rc != framesize) {
      fprintf(stderr, "unexpected short read\n");
      free(mapping);
      return NULL;
  }
  return mapping;
}

uint64_t key_record_t::size() const {
  return assemble64(ntohl(framesize_lo), ntohl(framesize_hi));
}
uint64_t key_record_t::offset() const {
  return assemble64(ntohl(offset_lo), ntohl(offset_hi));
}
double key_record_t::time() const {
  return assembleDouble(ntohl(time_lo), ntohl(time_hi));
}

static metadata_t * read_meta( const std::string& metafile, unsigned natoms,
                               bool with_invmass ) {

  metadata_t * meta = NULL;
  int meta_fd = open(metafile.c_str(), O_RDONLY|O_BINARY);
  size_t framesize=0;
  void *meta_mapping = read_file( meta_fd, 0, framesize );
  if (meta_mapping==NULL) {
    close(meta_fd);
    return meta;
  }
  BlobMap meta_blobs;
  try {
      meta_blobs = read_frame( meta_mapping, framesize );
  }
  catch (std::exception &e) {
      fprintf(stderr, "Reading metadata failed: %s\n", e.what());
      free(meta_mapping);
      close(meta_fd);
      return meta;
  }
  meta = new metadata_t;

  if (with_invmass && meta_blobs.find("INVMASS")!=meta_blobs.end()) {
    Blob blob=meta_blobs["INVMASS"];
    if (blob.count != natoms) {
      fprintf(stderr, "bad rmass count %d != %d\n", (int)blob.count, (int)natoms);
    } else {
      meta->invmass.resize(natoms);
      blob.get_float(&meta->invmass[0]);
    }
  }
  free(meta_mapping);
  close(meta_fd);
  return meta;
}

bool StkReader::recognizes(const std::string &path) {
      return path.size()>4 && 
             path.substr(path.size()-4)==".stk" &&
             isfile(path);
}

StkReader::StkReader(DtrReader *reader) {
  dtr=reader->path();
  framesets.push_back(reader);
  curframeset=0;
}

bool StkReader::init(const std::string &path, int * changed) {
  curframeset=0;
  dtr=path;

  if (changed) *changed = 0;
  /* process all the lines in the stk file */
  std::vector<std::string> fnames;
  std::ifstream input(path.c_str());
  if (!input) {
      fprintf(stderr, "Cannot open '%s' for reading\n", path.c_str());
      return false;
  }
  std::string fname;
  /* instantiate all the dtr readers */
  while (std::getline(input, fname)) {
      fnames.push_back(fname);
  }
  if (!fnames.size()) {
      fprintf(stderr, "Empty stk file\n");
      return false;
  }
  if (framesets.size()) {
      /* reloading an stk.  Find the dtrs that match the ones we already have,
       * and discard the rest */
      unsigned i=0; /* i will become the index of the last dtr that we've
                       already loaded */
      for (; i<fnames.size(); i++) {
          if (i==framesets.size() || fnames[i]!=framesets[i]->path()) break;
          if (getenv("DTRPLUGIN_VERBOSE")) 
              fprintf(stderr, "StkReader: Reusing dtr at %s\n", 
                      fnames[i].c_str());
      }
      /* delete any remaining framesets */
      for (unsigned j=i; j<framesets.size(); j++) delete framesets[j];
      framesets.erase(framesets.begin()+i, framesets.end());

      /* delete the filenames we've already loaded */
      fnames.erase(fnames.begin(), fnames.begin()+i);

      /* The set of overlapping frames may have changed!  Restore the keys
       * to their full, non-overlapping glory.  */
      for (i=0; i<framesets.size(); i++) {
          DtrReader * r = framesets[i];
          r->keys.restore_full_size();
      }
  }
  
  /* instantiate dtr readers */
  for (unsigned i=0; i<fnames.size(); i++) {
      DtrReader *reader = new DtrReader;
      if (getenv("DTRPLUGIN_VERBOSE"))
        fprintf(stderr, "StkReader: Loading timekeys from dtr at %s\n", 
                fnames[i].c_str());
      if (i>0) {
        const DtrReader * first = framesets[0];
        /* reuse information from earlier readers */
        reader->natoms = first->natoms;
        reader->with_velocity = first->with_velocity;
        reader->set_meta(first->get_meta());
      }
      if (!reader->init(fnames[i], NULL)) {
          delete reader;
          fprintf(stderr, "Failed opening frameset at %s\n", fnames[i].c_str());
          return false;
      }
      if (changed) *changed += 1;
      framesets.push_back(reader);
      if (i==0) this->with_velocity = reader->with_velocity;
  }

  natoms=framesets[0]->natoms;

  // now remove overlaps
  while (!framesets.back()->size()) {
      delete framesets.back();
      framesets.pop_back();
  }
  if (framesets.size()) {
    double first=framesets.back()->keys[0].time();
    size_t i=framesets.size()-1;
    while (i--) {
        /* find out how many frames to keep in frameset[i] */
        Timekeys& cur = framesets[i]->keys;
        size_t n = cur.size();
        while (n && cur[n-1].time() >= first) --n;
        cur.truncate( n );
        if (cur.size()) {
          double c0t = cur[0].time();
          first = (first < c0t) ? first : c0t;
        }
    }
  }
  return true;
}

size_t StkReader::size() const {
  size_t result=0;
  for (size_t i=0; i<framesets.size(); i++) 
    result += framesets[i]->keys.size();
  return result;
}

int StkReader::next(molfile_timestep_t *ts) {
  int rc=MOLFILE_EOF;
  while (curframeset < framesets.size() && 
         (rc=framesets[curframeset]->next(ts))==MOLFILE_EOF) {
    ++curframeset;
  }
  return rc;
}

const DtrReader * StkReader::component(size_t &n) const {
  for (size_t i=0; i<framesets.size(); i++) {
    size_t size = framesets[i]->size();
    if (n < size) return framesets[i];
    n -= size;
  }
  return NULL;
}

int StkReader::frame(size_t n, molfile_timestep_t *ts) const {
  const DtrReader *comp = component(n);
  if (!comp) return MOLFILE_EOF;
  return comp->frame(n, ts);
}

StkReader::~StkReader() {
  for (size_t i=0; i<framesets.size(); i++) 
    delete framesets[i];
}

std::string DtrReader::framefile(size_t n) const {
  return ::framefile( dtr, n, framesperfile(), ndir1(), ndir2() );
}

bool DtrReader::init(const std::string &path, int * changed) {
  dtr = path;
  /* Read the timekeys file */
  if (!keys.init(path)) return false;

  bool with_momentum = false;
  // read the first frame to see how many atoms there are, and whether 
  // there are any velocities.
  // Do this only if n_atoms isn't already set
  if (keys.size()>0 && natoms==0) {
    if (getenv("DTRPLUGIN_VERBOSE")) {
      fprintf(stderr, "reading first frame to get atom count\n");
    }
    std::string fname=::framefile(dtr, 0, keys.framesperfile(), 
            ndir1(), ndir2());
    int fd = open(fname.c_str(), O_RDONLY|O_BINARY);
    size_t framesize=0;
    unsigned i;
    void *mapping = read_file( fd, 0, framesize );
    if (mapping==NULL)  {
      fprintf(stderr, "Failed to find frame at %s\n", fname.c_str());
      close(fd);
      return false;
    }
    BlobMap blobs;
    try {
        blobs = read_frame(mapping, framesize);
    }
    catch (std::exception &e) {
        fprintf(stderr, "Warning: reading first frame failed, %s", e.what());
    }
    with_momentum = blobs.find("MOMENTUM")!=blobs.end();

    // I'm aware of three possible sources of positions: 
    //  "POSN" (the original frameset format)
    //  "POSITION" (the wrapped frameset formats)
    //  "POS" (anton trajectories)
    const char *posnames[] = { "POSN", "POSITION", "POS" };
    for (i=0; i<3; i++) {
      if (blobs.find(posnames[i])!=blobs.end()) {
        natoms = blobs[posnames[i]].count / 3;
        break;
      }
    }
    // similar for velocities
    const char *velnames[] = { "MOMENTUM", "VELOCITY" };
    for (i=0; i<2; i++) {
      if (blobs.find(velnames[i])!=blobs.end()) {
        with_velocity=true;
        break;
      }
    }
    free(mapping);
    close(fd);
  }

  if (natoms>0 && meta==NULL && owns_meta==false) {
    meta=read_meta( dtr + s_sep + "metadata",natoms, with_momentum );
    owns_meta = true;
  }

  /* we always reread the timekeys */
  if (changed) *changed = 1;
  return true;

}

size_t DtrReader::times(size_t start, size_t count, double *t) const {
    size_t remaining = keys.size()-start;
    count = (count < remaining) ? count : remaining;
    for (size_t j=0; j<count; j++) {
        t[j]=keys[start++].time();
    }
    return count;
}

size_t StkReader::times(size_t start, size_t count, double *t) const {
    size_t nread=0;
    size_t i=0,n=framesets.size();
    if (start<0) return 0;
    if (count<=0) return 0;
    /* Find the first frameset containing frames in the desired range */
    /* FIXME: could do this using a binary search... */
    for (; i<n; i++) {
        size_t sz = framesets[i]->size();
        if (start<sz) break;
        start -= sz;
    }
    /* Read times from framesets until count times are read. */
    for (; i<n; i++) {
        size_t sz = framesets[i]->times(start, count, t+nread);
        nread += sz;
        count -= sz;
        start=0;
        if (!count) break;
    }
    return nread;
}

static double dotprod(const double *x, const double *y) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

static void read_homebox( const double *box,
                          molfile_timestep_t *ts ) {

  ts->A = ts->B = ts->C = 0;

  double A[3] = { box[0], box[3], box[6] };
  double B[3] = { box[1], box[4], box[7] };
  double C[3] = { box[2], box[5], box[8] };

  // store lengths
  ts->A = sqrt(dotprod(A,A));
  ts->B = sqrt(dotprod(B,B));
  ts->C = sqrt(dotprod(C,C));

  if (ts->A == 0 || ts->B == 0 || ts->C == 0) {

    ts->alpha = ts->beta = ts->gamma = 90;

  } else {

    // compute angles
    double cosAB = dotprod(A,B)/(ts->A * ts->B);
    double cosAC = dotprod(A,C)/(ts->A * ts->C);
    double cosBC = dotprod(B,C)/(ts->B * ts->C);

    // clamp
    if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
    if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
    if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;

    // convert to angles using asin to avoid nasty rounding when we are
    // close to 90 degree angles.
    ts->alpha = 90.0 - asin(cosBC) * 90.0 / M_PI_2; /* cosBC */
    ts->beta  = 90.0 - asin(cosAC) * 90.0 / M_PI_2; /* cosAC */
    ts->gamma = 90.0 - asin(cosAB) * 90.0 / M_PI_2; /* cosAB */
  }
}

void write_homebox( const molfile_timestep_t * ts,
                    float * box ) {

  double A[3], B[3], C[3];

  // Convert VMD's unit cell information
  double cosBC = sin( ((90 - ts->alpha ) / 180) * M_PI );
  double cosAC = sin( ((90 - ts->beta  ) / 180) * M_PI );
  double cosAB = sin( ((90 - ts->gamma ) / 180) * M_PI );
  double sinAB = cos( ((90 - ts->gamma ) / 180) * M_PI );

  double Ax = ts->A;
  double Ay = 0;
  double Az = 0;
  double Bx = ts->B * cosAB;
  double By = ts->B * sinAB;
  double Bz = 0;
  double Cx,Cy,Cz;
  if (sinAB != 0) {
    Cx = cosAC;
    Cy = (cosBC - cosAC*cosAB) / sinAB;
    Cz = sqrt(1-Cx*Cx-Cy*Cy);
    Cx *= ts->C;
    Cy *= ts->C;
    Cz *= ts->C;
  } else {
    Cx=Cy=Cz=0;
  }
  A[0] = Ax; A[1] = Ay; A[2] = Az;
  B[0] = Bx; B[1] = By; B[2] = Bz;
  C[0] = Cx; C[1] = Cy; C[2] = Cz;

  // put vectors in column of homebox
  box[0] = A[0]; box[3] = A[1]; box[6] = A[2];
  box[1] = B[0]; box[4] = B[1]; box[7] = B[2];
  box[2] = C[0]; box[5] = C[1]; box[8] = C[2];
}

static int handle_wrapped_v2(
    BlobMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    molfile_timestep_t *ts ) {

  // just read POSITION in either single or double precision
  if (blobs.find("POSITION")==blobs.end()) {
    fprintf(stderr, "ERROR, Missing POSITION field in frame\n");
    return MOLFILE_ERROR;
  }
  Blob pos=blobs["POSITION"];
  if (pos.count != 3*natoms) {
    fprintf(stderr, "ERROR, Expected %d elements in POSITION; got %ld\n",
        3*natoms, pos.count);
    return MOLFILE_ERROR;
  }
  pos.get_float(ts->coords);

  if (with_velocity && ts->velocities && blobs.find("VELOCITY")!=blobs.end()) {
    Blob vel=blobs["VELOCITY"];
    if (vel.count != 3*natoms) {
      fprintf(stderr, "ERROR, Expected %d elements in VELOCITY; got %ld\n",
          3*natoms, vel.count);
      return MOLFILE_ERROR;
    }
    vel.get_float(ts->velocities);
  }

  if (blobs.find("UNITCELL")!=blobs.end()) {
    double box[9];
    blobs["UNITCELL"].get_double(box);
    read_homebox( box, ts );
  }
/*
#if defined(DESRES_READ_TIMESTEP2)
  if (blobs.find("ENERGY")!=blobs.end()) {
      blobs["ENERGY"].get_double(&ts->total_energy);
  }

  if (blobs.find("POT_ENERGY")!=blobs.end()) {
      blobs["POT_ENERGY"].get_double(&ts->potential_energy);
  }

  if (blobs.find("KIN_ENERGY")!=blobs.end()) {
      blobs["KIN_ENERGY"].get_double(&ts->kinetic_energy);
  }

  if (blobs.find("EX_ENERGY")!=blobs.end()) {
      blobs["EX_ENERGY"].get_double(&ts->extended_energy);
  }

  if (blobs.find("PRESSURE")!=blobs.end()) {
      blobs["PRESSURE"].get_double(&ts->pressure);
  }

  if (blobs.find("TEMPERATURE")!=blobs.end()) {
      blobs["TEMPERATURE"].get_double(&ts->temperature);
  }
#endif
*/
  return MOLFILE_SUCCESS;
}

namespace {

  inline void
  compute_center(int partition,
                 int nx, int ny, int nz,
                 float b0, float b1, float b2,
                 float b3, float b4, float b5,
                 float b6, float b7, float b8,
                 float* cx, float* cy, float* cz) {
    double nu, nv, nw, mu, mv, mw;
    double xc, yc, zc;

    // -----------------------------------------------
    // Map the partition number to its "mesh" position
    // (see define_mesh_collective in topology.c)
    // -----------------------------------------------
    int hmx = partition;
    int hmy  = hmx / nx;     /* y = y + ny*( z + nz*r ) */
    int hmz  = hmy / ny;     /* z = z + nz*r */
    hmx -= hmy * nx;         /* x = x */
    hmy -= hmz * ny;         /* y = y */

    nu = (double)nx;
    nv = (double)ny;
    nw = (double)nz;

    // -----------------------------------------------
    // Code adapted from configure_global_cell in
    // topology.c
    // -----------------------------------------------
    mu = -0.5*(nu-1) + (double)hmx;
    mv = -0.5*(nv-1) + (double)hmy;
    mw = -0.5*(nw-1) + (double)hmz;

    // We used to do FORCE_PRECISION(xc,float) here, but that
    // seems unnecessary in the context of trajectory writing.
    xc = b0*mu + b1*mv + b2*mw; 
    yc = b3*mu + b4*mv + b5*mw; 
    zc = b6*mu + b7*mv + b8*mw; 

    *cx = xc;
    *cy = yc;
    *cz = zc;
  }

  inline int 
  posn_momentum_v_1(int32_t nx, int32_t ny, int32_t nz,
                    uint64_t nparticles,
                    const double  * home_box,
                    const uint32_t* gid,
                    const uint32_t* npp,
                    const float   * rm, // reciprocal mass
                    const float* posn, const float* momentum,
                    /* returns */
                    float *position, float *velocity, double *box) {

    // bounding box is a straight multiple of the home box
    if (box) {
      box[0] = home_box[0]*nx;
      box[1] = home_box[1]*ny;
      box[2] = home_box[2]*nz;
        
      box[3] = home_box[3]*nx;
      box[4] = home_box[4]*ny;
      box[5] = home_box[5]*nz;

      box[6] = home_box[6]*nx;
      box[7] = home_box[7]*ny;
      box[8] = home_box[8]*nz;
    }


    int partition = 0;
    int remaining = 0;
    float cx = 0;
    float cy = 0;
    float cz = 0;
    float ux = home_box[0];
    float vx = home_box[1];
    float wx = home_box[2];
    float uy = home_box[3];
    float vy = home_box[4];
    float wy = home_box[5];
    float uz = home_box[6];
    float vz = home_box[7];
    float wz = home_box[8];

    for(uint64_t i=0; i<nparticles; ++i) {
      if (remaining == 0) {
        do {
          remaining = npp[partition];
          ++partition;
        } while (!remaining); // skip empty partitions
          compute_center(partition-1, nx,ny,nz, ux,vx,wx, uy,vy,wy, uz,vz,wz,
                        &cx,&cy,&cz);
      }
      uint32_t id = gid[i];
      if (id >= nparticles) {
        fprintf(stderr, "non-contiguous particles\n");
        return MOLFILE_ERROR;
      }

      if (posn) {
        float x = posn[3*i+0];
        float y = posn[3*i+1];
        float z = posn[3*i+2];

        position[3*id+0] = ux*x + vx*y + wx*z + cx;
        position[3*id+1] = uy*x + vy*y + wy*z + cy;
        position[3*id+2] = uz*x + vz*y + wz*z + cz;
      }

      if (velocity && momentum && rm) {
        velocity[3*id+0] = momentum[3*i+0]*rm[id];
        velocity[3*id+1] = momentum[3*i+1]*rm[id];
        velocity[3*id+2] = momentum[3*i+2]*rm[id];
      } else if (velocity) {
        velocity[3*id+0] = 0.0;
        velocity[3*id+1] = 0.0;
        velocity[3*id+2] = 0.0;
      }
      --remaining;
    }
    return MOLFILE_SUCCESS;
  }
}

static int handle_posn_momentum_v1(
    BlobMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    const float * rmass,
    molfile_timestep_t *ts ) {

  int32_t nx, ny, nz;
  double home_box[9], box[9];
  blobs["HOME_BOX"].get_double(home_box);
  blobs["NX"].get_int32(&nx);
  blobs["NY"].get_int32(&ny);
  blobs["NZ"].get_int32(&nz);
  
  std::vector<uint32_t> gid, npp;
  std::vector<float> pos, mtm;
  Blob gidblob=blobs["GID"];
  Blob nppblob=blobs["NPP"];
  Blob posblob=blobs["POSN"];
  Blob mtmblob=blobs["MOMENTUM"];

  if (gidblob.count != natoms) {
    fprintf(stderr, "Missing GID field\n");
    return MOLFILE_ERROR;
  }
  if (posblob.count != 3*natoms) {
    fprintf(stderr, "Missing POSN field\n");
    return MOLFILE_ERROR;
  }
  gid.resize(gidblob.count);
  npp.resize(nppblob.count);
  pos.resize(posblob.count);
  mtm.resize(mtmblob.count);

  gidblob.get_uint32(&gid[0]);
  nppblob.get_uint32(&npp[0]);
  posblob.get_float(&pos[0]);

  if (rmass && with_velocity) mtmblob.get_float(&mtm[0]);

  posn_momentum_v_1( nx, ny, nz, natoms, home_box, 
                     &gid[0], &npp[0], rmass,
                     &pos[0],
                     &mtm[0],
                     ts->coords,
                     ts->velocities,
                     box );

  read_homebox( box, ts );
  return MOLFILE_SUCCESS;
}

static int handle_wrapped_v1(
    BlobMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    molfile_timestep_t *ts ) {

  {
    // homebox
    double home_box[9], box[9];
    int32_t nx, ny, nz;
    blobs["HOME_BOX"].get_double(home_box);
    blobs["NX"].get_int32(&nx);
    blobs["NY"].get_int32(&ny);
    blobs["NZ"].get_int32(&nz);
    box[0] = home_box[0]*nx;
    box[1] = home_box[1]*ny;
    box[2] = home_box[2]*nz;
      
    box[3] = home_box[3]*nx;
    box[4] = home_box[4]*ny;
    box[5] = home_box[5]*nz;
 
    box[6] = home_box[6]*nx;
    box[7] = home_box[7]*ny;
    box[8] = home_box[8]*nz;
    read_homebox( box, ts );
  }

  Blob posblob=blobs["POSN"];
  Blob velblob=blobs["VELOCITY"];

  // get positions
  if (posblob.count != 3*natoms) {
    fprintf(stderr, "Missing POSN field\n");
    return MOLFILE_ERROR;
  }
  posblob.get_float(ts->coords);
  
  // if required, get velocities
  if (ts->velocities && velblob.count > 0) {
    if (velblob.count != 3*natoms) {
      fprintf(stderr, "VELOCITY field has %ld values; expected %d\n",
          velblob.count, 3*natoms);
      return MOLFILE_ERROR;
    }
    velblob.get_float(ts->velocities);
  }
  return MOLFILE_SUCCESS;
}

static int handle_anton_sfxp_v3(
    BlobMap &blobs,
    uint32_t natoms,
    bool with_velocity, 
    const float * rmass,
    molfile_timestep_t *ts ) {

  if (!rmass) {
    fprintf(stderr, "Cannot read anton_sfxp_v3 frame without rmass\n");
    return MOLFILE_ERROR;
  }

  double positionScale=0, momentumScale=0;
  // position scale...
  {
    Blob blob = blobs["POSITIONSCALE"];
    if (blob.count != 1) {
      fprintf(stderr, "Missing POSITIONSCALE field\n");
      return MOLFILE_ERROR;
    }
    blob.get_double(&positionScale);
  }
  // momentum scale
  if (ts->velocities) {
    Blob blob = blobs["MOMENTUMSCALE"];
    if (blob.count != 1) {
      fprintf(stderr, "Missing MOMENTUMSCALE field\n");
      return MOLFILE_ERROR;
    }
    blob.get_double(&momentumScale);
    momentumScale *= PEAKmassInAmu;
  }

  // box
  {
    double box[9] = { 0,0,0, 0,0,0, 0,0,0 };
    uint32_t anton_box[3];
    Blob boxblob = blobs["BOX"];
    if (boxblob.count != 3) {
      fprintf(stderr, "Missing BOX field\n");
      return MOLFILE_ERROR;
    }
    boxblob.get_uint32(anton_box);
    box[0] = sfxp_ulp32flt(anton_box[0])*positionScale;
    box[4] = sfxp_ulp32flt(anton_box[1])*positionScale;
    box[8] = sfxp_ulp32flt(anton_box[2])*positionScale;
    read_homebox( box, ts );
  }

  // velocities
  std::vector<int32_t> vel;
  if (ts->velocities) {
    Blob velblob = blobs["MOMENTUM"];
    if (velblob.count != 3*natoms) {
      fprintf(stderr, "Missing MOMENTUM field\n");
      return MOLFILE_ERROR;
    }
    vel.resize(3*natoms);
    velblob.get_int32(&vel[0]);
  }

  // positions
  std::vector<int32_t> pos(3*natoms);
  {
    Blob posblob = blobs["POS"];
    if (posblob.count != 3*natoms) {
      fprintf(stderr, "Missing POS field\n");
      return MOLFILE_ERROR;
    }
    posblob.get_int32(&pos[0]);
  }
  // convert and read into supplied storage
  for (unsigned i=0; i<natoms; i++) {
    ts->coords[3*i  ] = sfxp_ulp32flt(pos[3*i+0])*positionScale;
    ts->coords[3*i+1] = sfxp_ulp32flt(pos[3*i+1])*positionScale;
    ts->coords[3*i+2] = sfxp_ulp32flt(pos[3*i+2])*positionScale;
    if (ts->velocities) {
      const double rm = rmass[i] * momentumScale; // includes PEAKmassInAmu
      ts->velocities[3*i  ] = (float)(rm * sfxp_ulp32flt(vel[3*i  ]));
      ts->velocities[3*i+1] = (float)(rm * sfxp_ulp32flt(vel[3*i+1]));
      ts->velocities[3*i+2] = (float)(rm * sfxp_ulp32flt(vel[3*i+2]));
    }
  }
  return MOLFILE_SUCCESS;
}

int DtrReader::next(molfile_timestep_t *ts) {

  if (eof()) return MOLFILE_EOF;
  if (!ts) {
    ++m_curframe;
    return MOLFILE_SUCCESS;
  }
  size_t iframe = m_curframe;
  ++m_curframe;
  return frame(iframe, ts);
}

int DtrReader::ndir1() const {
  if (m_ndir1<0) DDgetparams(dtr, &m_ndir1, &m_ndir2);
  return m_ndir1;
}

int DtrReader::ndir2() const {
  if (m_ndir2<0) DDgetparams(dtr, &m_ndir1, &m_ndir2);
  return m_ndir2;
}

int DtrReader::frame(size_t iframe, molfile_timestep_t *ts) const {
  int rc = MOLFILE_SUCCESS;
  {
    off_t offset=0;
    size_t framesize=0;
    if (framesperfile() != 1) {
      offset = assemble64( ntohl(keys[iframe].offset_lo), 
                           ntohl(keys[iframe].offset_hi) );
      framesize = assemble64( ntohl(keys[iframe].framesize_lo), 
                              ntohl(keys[iframe].framesize_hi) );

    }
    ts->physical_time = keys[iframe].time();
    std::string fname=::framefile(dtr, iframe, framesperfile(), ndir1(), ndir2());
    int fd = open(fname.c_str(), O_RDONLY|O_BINARY);
    if (fd<0) return MOLFILE_EOF;
    void *mapping = read_file( fd, offset, framesize );
    if (mapping==NULL) {
      close(fd);
      return MOLFILE_EOF;
    }

    rc = frame_from_bytes( mapping, framesize, ts );

    free(mapping);
    close(fd);
  }
  return rc;
}

int DtrReader::frame_from_bytes(const void *buf, uint64_t len, 
                                molfile_timestep_t *ts) const {

  BlobMap blobs;
  try {
      blobs = read_frame(buf, len);
  }
  catch (std::exception &e) {
      fprintf(stderr, "Reading frame failed: %s\n", e.what());
      return MOLFILE_ERROR;
  }

  const float * rmass = meta && meta->invmass.size() ? 
      &meta->invmass[0] : NULL;

  // Now, dispatch to routines based on format
  std::string format = blobs["FORMAT"].str();
  if (format=="WRAPPED_V_2" || format == "DBL_WRAPPED_V_2") {
    return handle_wrapped_v2(blobs, natoms, with_velocity, ts);

  } else if (format=="POSN_MOMENTUM_V_1" || format=="DBL_POSN_MOMENTUM_V_1") {
    return handle_posn_momentum_v1(blobs, natoms, with_velocity, rmass, ts);

  } else if (format=="WRAPPED_V_1" || format == "DBL_WRAPPED_V_1") {
    return handle_wrapped_v1(blobs, natoms, with_velocity, ts);

  } else if (format=="ANTON_SFXP_V3") {
    return handle_anton_sfxp_v3(blobs, natoms, with_velocity, rmass, ts);
  }
  fprintf(stderr, "ERROR, can't handle format %s\n", format.c_str());
  return MOLFILE_ERROR;
}


namespace {
  struct meta_t {
    std::string label;
    std::string typecode;
    uint32_t elementsize;
    uint64_t count;
    const char *bytes;
    meta_t() {}
    meta_t(const std::string &l, const std::string &t, uint32_t e, uint32_t c,
           const void *b)
    : label(l), typecode(t), elementsize(e), count(c), 
    bytes(reinterpret_cast<const char *>(b)) {}
  };
  typedef std::vector<meta_t> MetaList;

  uint64_t typename_size(const MetaList &meta) {
    // just the set of distinct types
    uint64_t sz=0;
    typedef std::set<std::string> Typemap;
    Typemap types;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      types.insert(m->typecode);
    for (Typemap::const_iterator s=types.begin(); s!=types.end();++s)
      sz += s->size() + 1;
    sz += 1;
    return alignInteger(sz, s_alignsize);
  }

  uint64_t label_size(const MetaList &meta) {
    uint64_t sz=0;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      sz += m->label.size() + 1;
    sz += 1;
    return alignInteger(sz, s_alignsize);
  }

  uint64_t scalar_size(const MetaList &meta) {
    uint64_t sz=0;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      if (m->count <= 1) 
        sz += alignInteger( m->elementsize * m->count, s_alignsize );
    return sz;
  }
  uint64_t field_size(const MetaList &meta) {
    uint64_t sz=0;
    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m)
      if (m->count > 1) 
        sz += alignInteger( m->elementsize * m->count, s_alignsize );
    return sz;
  }

  void construct_frame( const std::vector<meta_t>& meta, 
                        std::vector<char>& bytes ) {
    uint64_t offset_header_block = 0;
    uint64_t size_header_block =
      alignInteger( sizeof(header_t), s_alignsize );

    uint64_t offset_meta_block = offset_header_block + size_header_block;
    uint64_t size_meta_block = 
      alignInteger( meta.size()*sizeof(metadisk_t), s_alignsize );

    uint64_t offset_typename_block = offset_meta_block + size_meta_block;
    uint64_t size_typename_block = typename_size(meta);

    uint64_t offset_label_block = offset_typename_block + size_typename_block;
    uint64_t size_label_block = label_size(meta);

    uint64_t offset_scalar_block = offset_label_block + size_label_block;
    uint64_t size_scalar_block = scalar_size(meta);

    uint64_t offset_field_block = offset_scalar_block + size_scalar_block;
    uint64_t size_field_block = field_size(meta);

    uint64_t offset_crc_block = offset_field_block + size_field_block;
    uint64_t size_crc_block = sizeof(uint32_t);

    uint64_t offset_padding_block = offset_crc_block + size_crc_block;
    uint64_t size_padding_block = 
      alignInteger(offset_padding_block,s_blocksize) - offset_padding_block;

    uint64_t framesize = offset_padding_block + size_padding_block;

    // construct the frame
    bytes.resize(framesize);
    char * base = &bytes[0];
    memset( base, 0, framesize );

    header_t *header = reinterpret_cast<header_t*>(base+offset_header_block);
    metadisk_t* diskmeta  = reinterpret_cast<metadisk_t*>(base+offset_meta_block);
    char*       typenames = reinterpret_cast<char*>(base+offset_typename_block);
    char*       labels    = reinterpret_cast<char*>(base+offset_label_block);
    char*       scalars   = reinterpret_cast<char*>(base+offset_scalar_block);
    char*       fields    = reinterpret_cast<char*>(base+offset_field_block);
    uint32_t*   crc       = reinterpret_cast<uint32_t*>(base+offset_crc_block);
    //char*       padding   = reinterpret_cast<char*>(base+offset_padding_block);

    /*** header ***/
    memset(header,0,sizeof(header_t));
    header->required.magic = htonl(magic_frame);
    header->required.version = htonl(s_version);

    header->required.framesize_lo = htonl(lobytes(framesize));
    header->required.framesize_hi = htonl(hibytes(framesize));

    header->size_header_block = htonl(size_header_block);
    header->unused0 = 0;
    uint64_t lrosetta = assemble64(s_lrosetta_lo,s_lrosetta_hi);
    header->irosetta = s_irosetta;
    header->frosetta = s_frosetta;

    header->drosetta_lo = lobytes(s_drosetta);
    header->drosetta_hi = hibytes(s_drosetta);
    header->lrosetta_lo = lobytes(lrosetta);
    header->lrosetta_hi = hibytes(lrosetta);

    header->endianism = htonl(machineEndianism());
    header->nlabels = htonl(meta.size());
    header->size_meta_block = htonl(size_meta_block);
    header->size_typename_block = htonl(size_typename_block);

    header->size_label_block = htonl(size_label_block);
    header->size_scalar_block = htonl(size_scalar_block);
    header->size_field_block_lo = htonl(lobytes(size_field_block));
    header->size_field_block_hi = htonl(hibytes(size_field_block));

    header->size_crc_block = htonl(size_crc_block);
    header->size_padding_block = htonl(size_padding_block);
    header->unused1 = 0;
    header->unused2 = 0;

    std::map<std::string,unsigned> typemap;

    for (MetaList::const_iterator m=meta.begin(); m!=meta.end(); ++m) {

      if (typemap.find(m->typecode)==typemap.end()) {
        unsigned code=typemap.size();
        typemap[m->typecode]=code;
        typenames=std::copy(m->typecode.begin(), m->typecode.end(), typenames);
        *typenames++ = 0;
      }

      diskmeta->type = htonl( typemap[m->typecode] );
      diskmeta->elementsize = htonl( m->elementsize );
      diskmeta->count_lo = htonl( lobytes( m->count ));
      diskmeta->count_hi = htonl( hibytes( m->count ));
      diskmeta++;

      labels=std::copy(m->label.begin(), m->label.end(), labels);
      *labels++ = 0;

      uint64_t nbytes = m->count*m->elementsize;
      if (m->count <= 1) {
        memcpy( scalars, m->bytes, nbytes );
        scalars += alignInteger( nbytes, s_alignsize );
      } else {
        memcpy( fields, m->bytes, nbytes );
        fields += alignInteger( nbytes, s_alignsize );
      }
    }
    *crc = fletcher(reinterpret_cast<uint16_t*>(base),offset_crc_block/2);
  }
}

void write_all( int fd, const char * buf, size_t count ) {
    while (count) {
        size_t n = ::write(fd, buf, count);
        if (n<0) {
            if (errno==EINTR) continue;
            throw std::runtime_error(strerror(errno));
        }
        buf += n;
        count -= n;
    }
}

bool DtrWriter::init(const std::string &path) {

  dtr=path;
  try {
    m_directory=path;
    char cwd[4096];

    while(m_directory.size() > 0 && m_directory[m_directory.size()-1] == s_sep) {
      m_directory.erase(m_directory.size()-1);
    }

    if (PathIsRelative(m_directory.c_str())) {
      if (! ::getcwd(cwd,sizeof(cwd))) {
        throw std::runtime_error(strerror(errno));
      }
      m_directory = std::string(cwd) + s_sep + m_directory;
    }

    recursivelyRemove(m_directory);
    ::DDmkdir(m_directory,0777,0, 0);

    // craft an empty clickme.dtr
    {
      std::string clickme_file = m_directory + s_sep + "clickme.dtr";
      FILE *fd = fopen(clickme_file.c_str(), "wb");
      fclose(fd);
    }

    // craft an empty metadata frame
    std::vector<meta_t> meta;
    std::vector<char> bytes;
    construct_frame( meta, bytes );

    {
      std::string metadata_file = m_directory + s_sep + "metadata";
      FILE *fd = fopen(metadata_file.c_str(), "wb");
      fwrite( &bytes[0], bytes.size(), 1, fd );
      fclose(fd);
    }

    // start writing timekeys file */
    std::string timekeys_path = dtr + s_sep + "timekeys";
    timekeys_file = fopen( timekeys_path.c_str(), "wb" );
    if (!timekeys_file) {
      fprintf(stderr, "Opening timekeys failed: %s\n", strerror(errno));
      return false;
    } else {
      key_prologue_t prologue[1];
      prologue->magic = htonl(magic_timekey);
      prologue->frames_per_file = htonl(frames_per_file);
      prologue->key_record_size = htonl(sizeof(key_record_t));
      fwrite( prologue, sizeof(key_prologue_t), 1, timekeys_file );
    }
  }
  catch (std::exception &e) {
    fprintf(stderr, "%s\n", e.what());
    return false;
  }
  return true;
}

int DtrWriter::next(const molfile_timestep_t *ts) {

  try {
    static const char *format = "WRAPPED_V_2";
    static const char *title = "written by VMD";
    float box[9];
    write_homebox( ts, box );

    double time = ts->physical_time;

    /* require increasing times (DESRESCode#1053) */
    if (last_time != HUGE_VAL && time <= last_time) {
      fprintf(stderr, 
          "dtrplugin: framesets require increasing times. previous %e, current %e\n", 
          last_time, time);
      return MOLFILE_ERROR;
    }

    std::vector<meta_t> meta;
    meta.push_back( 
        meta_t( "FORMAT", "char", sizeof(char), strlen(format), format ));
    meta.push_back( 
        meta_t( "TITLE", "char", sizeof(char), strlen(title), title ));
    meta.push_back( 
        meta_t( "CHEMICAL_TIME", "double", sizeof(double), 1, &time));
    meta.push_back( 
        meta_t( "UNITCELL", "float", sizeof(float), 9, box ));
    meta.push_back( 
        meta_t( "POSITION", "float", sizeof(float), 3*natoms, ts->coords ));
    if (ts->velocities) meta.push_back(
        meta_t( "VELOCITY", "float", sizeof(float), 3*natoms, ts->velocities ));
/*
#if defined(DESRES_READ_TIMESTEP2)
    meta.push_back(
        meta_t( "ENERGY", "double", sizeof(double), 1, &ts->total_energy ));
    meta.push_back(
        meta_t( "POT_ENERGY", "double", sizeof(double), 1, &ts->potential_energy ));
    meta.push_back(
        meta_t( "KIN_ENERGY", "double", sizeof(double), 1, &ts->kinetic_energy ));
    meta.push_back(
        meta_t( "EX_ENERGY", "double", sizeof(double), 1, &ts->extended_energy ));
    meta.push_back(
        meta_t( "TEMPERATURE", "double", sizeof(double), 1, &ts->temperature));
    meta.push_back(
        meta_t( "PRESSURE", "double", sizeof(double), 1, &ts->pressure));
#endif
*/
    std::vector<char> base;
    construct_frame(meta, base);
    uint64_t framesize = base.size();
    uint64_t keys_in_file = nwritten % frames_per_file;

    if (!keys_in_file) {
      if (frame_fd>0) ::close(frame_fd);

      framefile_offset = 0;
      std::string filepath=framefile(dtr, nwritten, frames_per_file, 0, 0);
      frame_fd = open(filepath.c_str(),O_WRONLY|O_CREAT|O_BINARY,0666);
      if (frame_fd<0) throw std::runtime_error(strerror(errno));
    }

    // write the data to disk
    write_all( frame_fd, &base[0], framesize );

    // add an entry to the keyfile list
    key_record_t timekey;
    timekey.time_lo = htonl(lobytes(time));
    timekey.time_hi = htonl(hibytes(time));
    timekey.offset_lo = htonl(lobytes(framefile_offset));
    timekey.offset_hi = htonl(hibytes(framefile_offset));
    timekey.framesize_lo = htonl(lobytes(framesize));
    timekey.framesize_hi = htonl(hibytes(framesize));

    if (fwrite(&timekey, sizeof(timekey), 1, timekeys_file)!=1) {
      fprintf(stderr, "Writing timekey failed\n");
      return MOLFILE_ERROR;
    }

#if defined(_MSC_VER)
    _commit(frame_fd);
#else
    fsync(frame_fd);
#endif
    fflush(timekeys_file);
#if defined(_MSC_VER)
    _commit(fileno(timekeys_file));
#else
    fsync(fileno(timekeys_file));
#endif

    ++nwritten;
    framefile_offset += framesize;
  } 
  catch (std::exception &e) {
    fprintf(stderr, "Write failed: %s\n",e.what());
    return MOLFILE_ERROR;
  }
  return MOLFILE_SUCCESS;
}

DtrWriter::~DtrWriter() {
  if (frame_fd>0) ::close(frame_fd);
  if (timekeys_file) fclose(timekeys_file);
}

/* compressed form: first write a -1 to indicate the new format.  Then
 * write number of ranges n, followed by (start, count) pairs. */
std::ostream& operator<<(std::ostream& out, const metadata_t& meta) {
    out << meta.invmass.size() << ' ';
    if (meta.invmass.size()) {
        out.write( (const char *)&meta.invmass[0], meta.invmass.size()*sizeof(meta.invmass[0]));
    }
    return out;
}

std::istream& operator>>(std::istream& in, metadata_t& meta) {
    uint32_t sz;
    char c;
    in >> sz;
    in.get(c);
    meta.invmass.resize(sz);
    if (sz) {
        in.read((char *)&meta.invmass[0], sz*sizeof(meta.invmass[0]));
    }
    return in;
}

std::ostream& DtrReader::dump(std::ostream &out) const {
  bool has_meta = meta!=NULL;
  out << SERIALIZED_VERSION << ' '
      << dtr << ' '
      << natoms << ' '
      << with_velocity << ' '
      << owns_meta << ' '
      << has_meta << ' ';
  if (owns_meta && has_meta) {
      out << *meta;
  }
  /* write raw m_ndir values so that we don't read them from .ddparams
   * if they haven't been read yet */
  out << m_ndir1 << ' '
      << m_ndir2 << ' ';
  keys.dump(out);
  return out;
}
std::istream& DtrReader::load(std::istream &in) {
  char c;
  bool has_meta;
  std::string version;
  in >> version;
  if (badversion(version)) {
    fprintf(stderr, "Bad version string\n");
    in.setstate( std::ios::failbit );
    return in;
  }
  in >> dtr
     >> natoms
     >> with_velocity
     >> owns_meta
     >> has_meta;
  if (owns_meta && has_meta) {
    delete meta;
    meta = new metadata_t;
    in.get(c);
    in >> *meta;
  }
  in >> m_ndir1
     >> m_ndir2;
  in.get(c);
  keys.load(in);
  return in;
}

std::ostream& StkReader::dump(std::ostream &out) const {
  out << dtr << ' '
      << framesets.size() << ' ';
  for (size_t i=0; i<framesets.size(); i++) framesets[i]->dump(out);
  return out;
}

std::istream& StkReader::load(std::istream &in) {
  in >> dtr;
  size_t size; in >> size; framesets.resize(size);
  char c; in.get(c);
  with_velocity=false;
  for (size_t i=0; i<framesets.size(); i++) {
    delete framesets[i];
    framesets[i] = new DtrReader;
    framesets[i]->load(in);
    if (i>0) framesets[i]->set_meta(framesets[0]->get_meta());
    if (i==0) with_velocity=framesets[i]->with_velocity;
  }
  if (framesets.size()) natoms=framesets[0]->natoms;
  return in;
}

///////////////////////////////////////////////////////////////////
//
// Plugin Interface
//
// ////////////////////////////////////////////////////////////////

void *
open_file_read( const char *filename, const char *filetype, int *natoms ) {

  FrameSetReader *h = NULL;
  std::string fname;

  // check for .stk file
  if (StkReader::recognizes(filename)) {
    h = new StkReader;

  } else {
    h = new DtrReader;
    // check for "clickme.dtr"
    fname=filename;
    std::string::size_type pos = fname.rfind( "clickme.dtr" );
    if (pos != std::string::npos) {
      fname.resize(pos);
      filename = fname.c_str();
    }
  }

  if (!h->init(filename)) {
    delete h;
    return NULL;
  }
  *natoms = h->natoms;
  return h;
}

int read_timestep_metadata(void *v, molfile_timestep_metadata *m) {
  FrameSetReader* handle = reinterpret_cast<FrameSetReader*>(v);
  m->has_velocities = handle->has_velocities();
  m->count = handle->size();
  return MOLFILE_SUCCESS;
}

int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  return h->next(ts);
}

#if defined(DESRES_READ_TIMESTEP2)
int read_timestep2(void *v, molfile_ssize_t n, molfile_timestep_t *ts) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  return h->frame(n, ts);
}

molfile_ssize_t read_times(void *v,
                                  molfile_ssize_t start, 
                                  molfile_ssize_t count,
                                  double * times) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  return h->times(start, count, times);
}
#endif
void close_file_read( void *v ) {
  FrameSetReader *h = reinterpret_cast<FrameSetReader *>(v);
  delete h;
}

void *open_file_write(const char *path, const char *type, int natoms) {
  DtrWriter *h = new DtrWriter(natoms);
  if (!h->init(path)) {
    delete h;
    h=NULL;
  }
  return h;
}

int write_timestep(void *v, const molfile_timestep_t *ts) {
  DtrWriter *h = reinterpret_cast<DtrWriter *>(v);
  return h->next(ts);
}

void close_file_write( void * v ) {
  DtrWriter *h = reinterpret_cast<DtrWriter *>(v);
  delete h;
}


#if defined(TEST_DTRPLUGIN)

int main(int argc, char *argv[]) {

  /* check input arguments */
  if (argc != 3) {
    fprintf(stderr, "Usage: %s input.dtr output.dtr\n", argv[0]);
    return 1;
  }
  int natoms;
  
  /* read in all frames */
  void *handle = open_file_read( argv[1], "dtr", &natoms);
  printf("got %d atoms\n", natoms);
  molfile_timestep_t ts[1];
  std::vector<float> timesteps;
  do {
    timesteps.resize( timesteps.size()+3*natoms );
    ts->coords = &timesteps[ timesteps.size() - 3*natoms ];
    ts->velocities = NULL;
  } while (read_next_timestep(handle, natoms, ts)==MOLFILE_SUCCESS);
  int nframes = timesteps.size()/(3*natoms) - 1;

  printf("read %d steps\n", nframes );
  close_file_read(handle);

  /* write out all frames */
  handle = open_file_write( argv[2], "dtr", natoms);
  if (!handle) return 1;
  ts->coords = &timesteps[0];
  for (int i=0; i<nframes; i++) {
    if (write_timestep( handle, ts )!=MOLFILE_SUCCESS) {
      fprintf(stderr, "failed to write timestep %d/%d\n", i+1, nframes);
      return 1;
    }
    ts->coords += 3*natoms;
  }
  printf("wrote %d steps\n", nframes);
  open_file_write(handle);

  /* now try to read it back in */
  int new_natoms;
  handle = open_file_read( argv[2], "dtr", &new_natoms);
  if (handle) return 1;
  if (new_natoms != natoms) {
    fprintf(stderr, "number of atoms changed: %d -> %d\n", natoms, new_natoms);
    return 1;
  }
  close_file_read(handle);
  return 0;
}

#endif

#if defined(TEST_DTR_DUMP)
int main(int argc, char *argv[]) {

  StkReader src, dst;
  src.init(argv[1]);

  std::ostringstream out;
  src.dump(out);
  assert(out);

  std::istringstream in(out.str());
  dst.load(in);
  assert(in);

  return 0;
}
#endif

