#include <iostream>

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
uint8_t njets = 0;

int main(int argc, char* argv[]) {
  for (bool is_mc : {false,true}) {
    reader read(is_mc ? "hgam_mc.dat" : "hgam_data.dat");

    char dm;
    read(dm);
    if ((is_mc && dm!='m') || (!is_mc && dm!='d')) {
      cerr << (is_mc ? "mc" : "data")
           << " file begins with \'" << dm << "\'" << endl;
      return 1;
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
        TEST(ent)
        if (is_mc) read(weight);
        read(ph[0]);
        read(ph[1]);
        diph = ph[0] + ph[1];
        TEST(diph.m())
        read(njets);
        TEST((unsigned)njets)
        if (njets>4) njets = 4;
        // jets.resize(njets);
        for (decltype(njets) i=0; i<njets; ++i) {
          read(jets[i]);
          TEST(jets[i].pt());
        }

        if (ent>=4u) break;
      }
    }
  }
}
