#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>

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

const double inf = std::numeric_limits<double>::infinity();

float lumi=0, weight=1;
uint32_t nevents_total = 0;
bool is_mc, need_jets=false;
uint8_t njets=0, njets_stored=0;

// ==================================================================
vec4 y[2], jets[4], yy;

auto nj() noexcept { return njets; }

#include "varfcns.hh"
// ==================================================================

struct bin_t {
  double w = 0, w2 = 0;
  void operator++() noexcept {
    ++w; ++w2;
  }
  void operator+=(double weight) noexcept {
    w  += weight;
    w2 += weight*weight;
  }
  bin_t& operator+=(const bin_t& b) noexcept {
    w  += b.w;
    w2 += b.w2;
    return *this;
  }
};

template <typename E>
size_t find_bin(double x, const E& edges) {
  return std::distance(
    edges.begin(), std::upper_bound(edges.begin(), edges.end(), x)
  );
}

int main(int argc, char* argv[]) {
  if (argc!=5) {
    cout << "usage: " << argv[0] << " data.dat mc.dat bins.txt out.json\n";
    return 1;
  }

  struct vardef {
    double(*f)();
    std::vector<double> edges;
    std::string name;
  };
  std::vector<vardef> vars;
  { std::ifstream f(argv[3]);
    size_t g = 0;
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
        if (std::isspace(c)||c==','||c=='('||c==')'||c=='\0') {
          if (b > a) {
            auto x = line.substr(a,b+1-a);
            if (x=="inf" || x=="infty" || x=="∞")
              var.edges.push_back(inf);
            else if (x=="-inf" || x=="-infty" || x=="-∞")
              var.edges.push_back(-inf);
            else
              var.edges.push_back(stod(x));
          }
          if (!c) break;
          else if (c=='(') g = var.edges.size(); // (n min max)
          else if (c==')') {
            if (var.edges.size()-g!=3) {
              cerr << "\033[31minvalid edges definition\033[0m\n";
              cerr << line << endl;
              return 1;
            }
            unsigned n = var.edges[g];
            const double a = var.edges[g+1], b = var.edges[g+2], d = (b-a)/n;
            var.edges.resize(g);
            var.edges.push_back(a);
            for (unsigned i=1; i<n; ++i)
              var.edges.push_back(a+d*i);
            var.edges.push_back(b);
            g = 0;
          }
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
  std::vector<bin_t> data(nbins), mc(nbins);

  for (const char* fname : {argv[1],argv[2]}) {
    reader read(fname);

    char dm;
    switch (read(dm)) {
      case 'm': is_mc = true; break;
      case 'd': is_mc = false; break;
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

    auto& bins = is_mc ? mc : data;
    { ivanp::timed_counter<> ent;
      for (; read; ++ent) {
        if (is_mc) read(weight);
        read(y[0]);
        read(y[1]);
        yy = y[0] + y[1];
        read(njets);
        njets_stored = njets>4 ? 4 : njets;
        if (need_jets) {
          for (decltype(njets) i=0; i<njets_stored; ++i)
            read(jets[i]);
        } else {
          read.skip(sizeof(float[4])*njets_stored);
        }
        // ----------------------------------------------------------
        size_t bin = 0;
        for (size_t i=0; i<vars.size(); ++i) {
          size_t b = find_bin(vars[i].f(),vars[i].edges);
          if (b==0 || b==vars[i].edges.size()) goto next_event;
          --b;
          if (b==0) continue;
          for (size_t j=0; j<i; ++j) b *= (vars[j].edges.size()-1);
          bin += b;
        }
        bins[bin] += weight;
        // ----------------------------------------------------------
next_event: ;
      }
      if (ent!=nevents_total) {
        cerr << "\033[31m" << nevents_total << " expected, "
          << ent << " events read\033[0m" << endl;
      }
    }
  }

  // write output ---------------------------------------------------
  std::ofstream out(argv[4]);
  out << "{\"lumi\": " << lumi << ",\n";
  out << "\"bins\":[\n";
  bool first = true;
  for (const auto& var : vars) {
    if (!first) out << ",\n";
    out << "[\"" << var.name << "\",[\n";
    first = true;
    for (const auto& x : var.edges) {
      if (first) first = false;
      else out << ',';
      out << x;
    }
    out << "\n]]";
  }
  out << '\n';
  for (bool is_mc : {false,true}) {
    out << "],\n\"" << (is_mc ? "mc" : "data") << "\":[\n";
    first = true;
    for (const auto& bin : (is_mc ? mc : data)) {
      if (first) first = false;
      else out << ",\n";
      out << '[' << bin.w << ',' << std::sqrt(bin.w2) << ']';
    }
    out << '\n';
  }
  out << "]\n}" << std::flush;
}
