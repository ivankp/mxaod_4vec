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

const unsigned nmax = 1000;
const double jetR = 0.4;

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

// ==================================================================
vec4<> y[2], yy;
std::vector<vec4<>> j;

auto nj() noexcept { return j.size(); }
bool nj(unsigned n) noexcept { return j.size() >= n; }

double f_HT_jets() noexcept {
  double HT = 0;
  for (const auto& j : j) HT += j.pt();
  return HT;
}
double f_HT_jets_yy() noexcept { return f_HT_jets() + yy.pt(); }
// ==================================================================

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cerr << "usage: " << argv[0] << " file.dat\n";
    return 1;
  }

  struct fcn_t {
    double(*f)();
    bool need_jets;
  };
  const std::map<const char*,fcn_t,less_str> fcns {
    { "pT_yy", { []{ return yy.pt(); }, false } },
    { "m_yy",  { []{ return yy.m(); }, false } },
    { "pT_y1", { []{ return y[0].pt(); }, false } },
    { "pT_y2", { []{ return y[1].pt(); }, false } },
    { "rat_pT_y1_y2", { []{ return y[0].pt()/y[1].pt(); }, false } },
    { "eta_y1", { []{ return y[0].eta(); }, false } },
    { "eta_y2", { []{ return y[1].eta(); }, false } },
    { "y_y1", { []{ return y[0].rap(); }, false } },
    { "y_y2", { []{ return y[1].rap(); }, false } },
    { "dy_y1_y2", { []{ return std::abs(y[0].rap()-y[1].rap()); }, false } },

    { "Njets", { []{ return (double)nj(); }, true } },
    { "pT_j1", { []{ return nj(1) ? j[0].pt() : 0; }, true } },
    { "pT_j2", { []{ return nj(2) ? j[1].pt() : 0; }, true } },
    { "pT_j3", { []{ return nj(3) ? j[2].pt() : 0; }, true } },
    { "eta_j1", { []{ return nj(1) ? j[0].eta() : 0; }, true } },
    { "eta_j2", { []{ return nj(2) ? j[1].eta() : 0; }, true } },
    { "eta_j3", { []{ return nj(3) ? j[2].eta() : 0; }, true } },
    { "y_j1", { []{ return nj(1) ? j[0].rap() : 0; }, true } },
    { "y_j2", { []{ return nj(2) ? j[1].rap() : 0; }, true } },
    { "y_j3", { []{ return nj(3) ? j[2].rap() : 0; }, true } },
    { "dy_j1_j2", { []{ return nj(2) ? std::abs(j[0].rap()-j[1].rap()) : 0; }, true } },
    { "m_jj", { []{ return nj(2) ? (j[0]+j[1]).m() : 0; }, true } },

    { "HT_jets", { f_HT_jets, true } },
    { "HT_jets_yy", { f_HT_jets_yy, true } },

    { "x_yy", { []{ return yy.pt()/f_HT_jets_yy(); }, true } },
    { "x_j1", { []{ return nj(1) ? j[0].pt()/f_HT_jets_yy() : 0; }, true } },
    { "x_j2", { []{ return nj(2) ? j[1].pt()/f_HT_jets_yy() : 0; }, true } },
    { "x_j3", { []{ return nj(3) ? j[2].pt()/f_HT_jets_yy() : 0; }, true } },

    { "pT_miss", { []{
        auto p = yy;
        for (const auto& j : j) p += j;
        return p.pt();
      }, true } },
    { "s2", { []{
        auto p = yy;
        for (const auto& j : j) p += j;
        return p.m();
      }, true } },
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
