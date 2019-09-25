#include <cstring>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <iostream>
#include <vector>

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

const unsigned nmax = 1000;
const double jetR = 0.4;

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

// ==================================================================
ivanp::vec4<> y[2], yy;
std::vector<ivanp::vec4<>> j;

auto nj() noexcept { return j.size(); }

#include "varfcns.hh"
// ==================================================================

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cerr << "usage: " << argv[0] << " file.dat\n";
    return 1;
  }

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
  std::vector<std::function<bool()>> cuts;
  const auto& req_cuts = req["cuts"];
  cuts.reserve(req_cuts.size());
  for (const auto& cut : req_cuts) {
    const auto& fcn = fcns.at(cut.at(0).get_ref<const std::string&>().c_str());
    const auto& op = cut.at(1).get_ref<const std::string&>();
    bool lt = false;
    if (op=="l") lt = true; else
    if (op!="g") throw error("unexpected cut operator \"",op,"\"");
    const double x = cut.at(2).get<double>();
    // cuts.emplace_back([=,f=fcn.f]{ return (f() < x) == cmp; });
    cuts.emplace_back([=,f=fcn.f]{ return lt ? (f() < x) : (f() > x); });
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

  unsigned nselected = 0;
  uint32_t runNumber;
  uint64_t eventNumber;
  uint32_t njets;
  float mom[4];

  bool first = true;
  cout << '[';
  for (;dat;) {
    dat >> runNumber >> eventNumber;
    dat >> mom;
    y[0] = { mom, vec4<>::PtEtaPhiM };
    dat >> mom;
    y[1] = { mom, vec4<>::PtEtaPhiM };
    yy = y[0] + y[1];
    dat >> njets;
    if (need_jets) {
      j.resize(njets);
      for (decltype(njets) i=0; i<njets; ++i) {
        dat >> mom;
        j[i] = { mom, vec4<>::PtEtaPhiM };
      }
    } else {
      dat.skip(sizeof(mom)*njets);
    }

    for (const auto& cut : cuts)
      if (!cut()) goto next_event;

    if (++nselected <= nmax) {
      if (first) first = false;
      else cout << ',';
      cout << '[' << runNumber << ',' << eventNumber;
      for (const auto& var : vars)
        cout << ',' << var();
      cout << ']';
    }

next_event: ;
  }
  if (nselected > nmax) cout << ',' << (nselected-nmax);
  cout << "]";
}
