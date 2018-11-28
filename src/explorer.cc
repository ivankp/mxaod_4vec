#include <cstring>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <iostream>
#include <vector>
#include <map>

#include <nlohmann/json.hpp>

#include "ivanp/math/vec4.hh"
#include "ivanp/error.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;
using nlohmann::json;
using namespace ivanp;

struct less_str {
  bool operator()(const char* a, const char* b) const noexcept {
    return strcmp(a,b) < 0;
  }
};

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
} e;

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cerr << "usage: " << argv[0] << " file.dat\n";
    return 1;
  }

  struct fcn_t {
    double(*f)();
    bool need_jets;
  };
  std::map<const char*,fcn_t,less_str> fcns {
    { "pT_y1", { []{ return e.y[0].pt(); }, false } },
    { "pT_y2", { []{ return e.y[1].pt(); }, false } },
    { "pT_yy", { []{ return (e.y[0]+e.y[1]).pt(); }, false } },
    { "m_yy",  { []{ return (e.y[0]+e.y[1]).m(); }, false } },
    { "rat_pT_y1_y2", { []{ return e.y[0].pt()/e.y[1].pt(); }, false } }
  };

  json req;
  std::cin >> req;

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

  bool need_jets = false;
  const auto& req_cuts = req.at("cuts");
  std::vector<std::function<bool()>> cuts;
  cuts.reserve(req_cuts.size());
  for (const auto& cut : req_cuts) {
    const auto& fcn = fcns.at(cut.at(0).get_ref<const std::string&>().c_str());
    const auto& op = cut.at(1).get_ref<const std::string&>();
    bool cmp = false;
    if (op=="l") cmp = true; else
    if (op!="g") throw error("unexpected cut operator \"",op,"\"");
    const double x = cut.at(2).get<double>();
    cuts.emplace_back([=,f=fcn.f]{ return (f() < x) == cmp; });
    if (fcn.need_jets) need_jets = true;
  }

  std::vector<double(*)()> vars;
  const auto& req_vars = req.at("vars");
  vars.reserve(req_vars.size());
  for (const auto& var : req_vars) {
    const auto& fcn = fcns.at(var.get_ref<const std::string&>().c_str());
    vars.emplace_back(fcn.f);
    if (fcn.need_jets) need_jets = true;
  }

  uint32_t runNumber;
  uint64_t eventNumber;
  uint32_t njets;
  float mom[4];

  bool first = true;
  cout << '[';
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

    for (const auto& cut : cuts)
      if (!cut()) goto next_event;

    if (first) first = false;
    else cout << ',';
    cout << '[' << runNumber << ',' << eventNumber;
    for (const auto& var : vars)
      cout << ',' << var();
    cout << "]";

next_event: ;
  }
  cout << "]";
}
