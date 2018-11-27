#include <iostream>
#include <iomanip>
#include <array>

#include "ivanp/io/mem_file.hh"
#include "ivanp/scribe.hh"
#include "ivanp/math/vec4.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;

using namespace ivanp;

int main(int argc, char* argv[]) {
  if (argc!=3 && argc!=2) {
    cout << "usage: " << argv[0] << " file.dat[.xz]\n";
    return 1;
  }

  const mem_file file = (
    ends_with(argv[1],".xz") || ends_with(argv[1],".lzma")
    ? mem_file::pipe(cat("unxz -c ",argv[1]).c_str())
    : mem_file::mmap(argv[1])
  );

  scribe::reader sr(file.mem(),file.size());
  // const auto& head = sr.head();
  // TEST(head);

  cout << "[["
    "\"runNumber\",\"eventNumber\","
    "\"γγ pT [GeV]\",\"γγ m [GeV]\",\"pT γ1/pT γ2\""
  "]";

  const auto event_t = sr.get_type()[0][0];
  const auto photons_i = event_t.index("photons");
  const auto runNumber_i = event_t.index("runNumber");
  const auto eventNumber_i = event_t.index("eventNumber");

  using vec = vec4<double>;
  std::array<vec,2> y;

  for (auto event : sr[0]) {
    const auto& photons = event[photons_i].cast<float(&)[2][4]>();
    y = {{
      { photons[0], vec::PtEtaPhiM },
      { photons[1], vec::PtEtaPhiM }
    }};

    const auto yy = y[0] + y[1];
    const auto pT_yy = yy.pt();
    const auto m_yy = yy.m();

    if (pT_yy < 250) continue;
    if (m_yy < 121 || 129 < m_yy) continue;

    cout << ",\n["
      << event[runNumber_i].cast<uint32_t>() << ','
      << event[eventNumber_i].cast<uint64_t>() << ','
      << pT_yy << ','
      << m_yy << ','
      << (photons[0][0]/photons[1][0]) << ']';
  }
  cout << "]";
}
