#ifndef READER2_HH
#define READER2_HH

#include <cstring>
#include "ivanp/io/mem_file.hh"
#include "ivanp/math/vec4.hh"

class reader {
  ivanp::mem_file f;
  const char *pos, *end;

public:
  reader(const char* filename)
  : f(ivanp::mem_file::mmap(filename)), pos(f.mem()), end(pos+f.size()) { }

  operator bool() const noexcept { return pos != end; }

  void skip(size_t len) { pos += len; }

  template <typename T>
  T& operator()(T& x) {
    memcpy(&x,pos,sizeof(T));
    pos += sizeof(T);
    return x;
  }
  ivanp::vec4<>& operator()(ivanp::vec4<>& x) {
    float mom[4];
    return x = { operator()(mom), ivanp::vec4<>::PtEtaPhiM_t{} };
  }
};

#endif
