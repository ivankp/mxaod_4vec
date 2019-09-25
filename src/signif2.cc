#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "ivanp/io/mem_file.hh"
#include "ivanp/math/vec4.hh"
#include "ivanp/timed_counter.hh"
#include "reader2.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;

using vec4 = ivanp::vec4<double>;

float lumi=0, weight=1;
uint32_t nevents_total = 0;
bool is_mc, need_jets=false;
uint8_t njets = 0;

// ==================================================================
vec4 y[2], j[4], yy;

auto nj() noexcept { return njets; }

#include "varfcns.hh"
// ==================================================================

struct bin {
  double w = 0, w2 = 0;
  void operator++() noexcept {
    ++w; ++w2;
  }
  void operator+=(double weight) noexcept {
    w  += weight;
    w2 += weight*weight;
  }
  bin& operator+=(const bin& b) noexcept {
    w  += b.w;
    w2 += b.w2;
    return *this;
  }
};

int main(int argc, char* argv[]) {
  if (argc!=4) {
    cout << "usage: " << argv[0] << " data.dat mc.dat bins.txt\n";
    return 1;
  }

  std::vector<bin> data, mc;
  struct vardef {
    double(*f)();
    std::vector<double> edges;
    std::string name;
  };
  std::vector<vardef> vars;
  { std::ifstream f(argv[3]);
    for (std::string line; getline(f,line); ) {
      if (line.empty()||line[0]=='#') continue;
      const size_t col = line.find(':');
      if (col==std::string::npos) continue;
      vars.emplace_back();
      auto& var = vars.back();
      var.name = line.substr(0,col);
      size_t a = col + 1, b = a;
      for (;;) {
        char c = line[b];
        if (std::isspace(c)||c==','||c=='\0') {
          if (b > a) var.edges.push_back(stod(line.substr(a,b+1-a)));
          if (!c) break;
          a = ++b;
        } else ++b;
      }
    }
  }

  size_t nbins = 1;
  for (auto& var : vars) {
    cout << var.name << ":";
    for (const auto& x : var.edges)
      cout << " " << x;
    cout << endl;
    try {
      const auto& fcn = fcns.at(var.name.c_str());
      var.f = fcn.f;
      need_jets |= fcn.need_jets;
    } catch (...) {
      cerr << "\033[31mvariable \""<< var.name <<"\" is not defined\033[0m\n";
      return 1;
    }
    nbins *= (var.edges.size()-1);
  }
  TEST(nbins)
  data.resize(nbins);
  mc  .resize(nbins);

  return 0;

  for (const char* fname : {argv[1],argv[2]}) {
    reader read(fname);

    char dm;
    switch (read(dm)) {
      case 'm': is_mc = true;
      case 'd': is_mc = false;
      default: {
        cerr << "\033[31mfile \"" << fname << "\" starts with \'" << dm
          << "\' instead of \'d\' or \'m\'\033[0m\n";
        return 1;
      }
    }
    if (!is_mc) {
      read(lumi);
      TEST(lumi);
      weight = 1;
    }
    read(nevents_total);
    TEST(nevents_total);

    { ivanp::timed_counter<> ent;
      for (; read; ++ent) {
        if (is_mc) read(weight);
        read(y[0]);
        read(y[1]);
        yy = y[0] + y[1];
        read(njets);
        if (njets>4) njets = 4;
        for (decltype(njets) i=0; i<njets; ++i) {
          read(j[i]);
        }
      }
      if (ent!=nevents_total) {
        cerr << "\033[31m" << nevents_total << " expected, "
          << ent << " events read\033[0m" << endl;
      }
    }
  }
}
