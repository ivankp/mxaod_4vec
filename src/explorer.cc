#include <cstring>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <iostream>
#include <vector>

#include "ivanp/math/vec4.hh"
#include "ivanp/error.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;

using namespace ivanp;

class file {
  char *m, *pos, *end;
public:
  file(const char* name) {
    struct stat sb;
    int fd = open(name, O_RDONLY);
    if (fd == -1) throw error("open");
    if (fstat(fd, &sb) == -1) throw error("fstat");
    if (!S_ISREG(sb.st_mode)) throw error("not a file");
    size_t len = sb.st_size;
    m = reinterpret_cast<char*>(mmap(0,len,PROT_READ,MAP_SHARED,fd,0));
    if (m == MAP_FAILED) throw error("mmap");
    if (close(fd) == -1) throw error("close");
    pos = m;
    end = m + len;
  }
  ~file() { munmap(m,end-m); }
  char get() { return *(pos++); }
  template <typename T>
  friend inline file& operator>>(file& f, T& x) {
    memcpy(reinterpret_cast<char*>(&x),f.pos,sizeof(x));
    f.pos += sizeof(x);
    return f;
  }
  operator bool() const { return pos!=end; }
  void skip(size_t off) { pos += off; }
};

struct event {
  vec4<> y[2];
  std::vector<vec4<>> jets;
};

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cout << "usage: " << argv[0] << " file.dat\n";
    return 1;
  }

  file dat(argv[1]);
  for (int nbraces=0;;) { // skip header
    if (!dat) {
      cerr << "file ended in header\n";
      return 1;
    }
    const char c = dat.get();
    if (c=='{') ++nbraces;
    else if (c=='}') --nbraces;
    if (!nbraces) break;
  }

  cout << "[["
    "\"runNumber\",\"eventNumber\","
    "\"γγ pT [GeV]\",\"γγ m [GeV]\",\"pT γ1/pT γ2\""
  "]";

  uint32_t runNumber;
  uint64_t eventNumber;
  uint32_t njets;
  float mom[4];
  event e;

  bool need_jets = false;

  for (;dat;) {
    dat >> runNumber >> eventNumber;
    dat >> mom;
    e.y[0] = { mom, vec4<>::PtEtaPhiM };
    dat >> mom;
    e.y[1] = { mom, vec4<>::PtEtaPhiM };
    dat >> njets;
    if (need_jets) {
      e.jets.resize(njets);
      for (decltype(njets) i=0; i<njets; ++i) {
        dat >> mom;
        e.jets[i] = { mom, vec4<>::PtEtaPhiM };
      }
    } else {
      dat.skip(sizeof(mom)*njets);
    }

    const auto yy = e.y[0] + e.y[1];
    const auto pT_yy = yy.pt();
    const auto m_yy = yy.m();

    if (pT_yy < 250) continue;
    if (m_yy < 121 || 129 < m_yy) continue;

    const auto ptrat = e.y[0].pt() / e.y[1].pt();

    cout << ",\n["
      << runNumber << ','
      << eventNumber << ','
      << pT_yy << ','
      << m_yy << ','
      << ptrat << ']';
  }
  cout << "]";
}
