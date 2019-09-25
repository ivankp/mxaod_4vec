#include <iostream>
#include <fstream>
#include <type_traits>
#include <vector>

#include "ivanp/math/vec4.hh"
#include "ivanp/timed_counter.hh"
#include "reader2.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;

using vec4 = ivanp::vec4<double>;
using mom_t = float[4];

float lumi=0, weight=1;
uint32_t nevents_total=0, nevents=0;
mom_t mom;
uint8_t njets = 0;

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "usage: " << argv[0] << " in.dat out.dat\n";
    return 1;
  }
  reader read(argv[1]);

  std::ofstream out(argv[2]);
  auto write = [&out](const auto& x){
    out.write(reinterpret_cast<const char*>(&x),sizeof(x));
  };

  char dm;
  const bool is_mc = read(dm) == 'm';
  TEST(dm)
  write(dm);
  if (!is_mc) {
    read(lumi);
    TEST(lumi);
    write(lumi);
    weight = 1;
  }
  read(nevents_total);
  TEST(nevents_total);

  const auto nevents_pos = out.tellp();
  write(nevents);

  mom_t ph[2];
  { ivanp::timed_counter<> ent;
    for (; read; ++ent) {
      if (is_mc) read(weight);
      read(ph);
      read(njets);

      const vec4 diph
        = vec4(ph[0],vec4::PtEtaPhiM_t{})
        + vec4(ph[1],vec4::PtEtaPhiM_t{});
      const double myy = diph.m();

      if (121<myy && myy<129) {
        ++nevents;
        if (is_mc) write(weight);
        write(ph[0]);
        write(ph[1]);
        write(njets);

        if (njets>4) njets = 4;
        for (decltype(njets) i=0; i<njets; ++i)
          write(read(mom));
      } else {
        if (njets>4) njets = 4;
        read.skip(sizeof(mom_t)*njets);
      }
    }
    if (ent!=nevents_total) {
      cerr << "\033[31m " << nevents_total << " expected, "
        << ent << " events read" << endl;
    }
  }
  TEST(nevents)

  out.flush();
  out.seekp(nevents_pos);
  write(nevents);
  out.flush();
}
