#include <iostream>
#include <type_traits>
#include <vector>

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
bool is_mc;
uint8_t njets = 0;

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "usage: " << argv[0] << " data.dat mc.dat\n";
    return 1;
  }
  for (const char* fname : {argv[1],argv[2]}) {
    reader read(fname);

    char dm;
    switch (read(dm)) {
      case 'm': is_mc = true;
      case 'd': is_mc = false;
      default: {
        cerr << "file \"" << fname << "\" starts with \'" << dm
          << "\' instead of \'d\' or \'m\'\n";
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

    vec4 ph[2], jets[4], diph;
    { ivanp::timed_counter<> ent;
      for (; read; ++ent) {
        if (is_mc) read(weight);
        read(ph[0]);
        read(ph[1]);
        diph = ph[0] + ph[1];
        read(njets);
        if (njets>4) njets = 4;
        for (decltype(njets) i=0; i<njets; ++i) {
          read(jets[i]);
        }
      }
      if (ent!=nevents_total) {
        cerr << "\033[31m " << nevents_total << " expected, "
          << ent << " events read" << endl;
      }
    }
  }
}
