#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <string>

#include <TFile.h>

#include "ivanp/timed_counter.hh"
#include "ivanp/error.hh"
#include "ivanp/root/branch_reader.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;
using ivanp::error;

using size_type = uint32_t;
using float_t = float;

template <typename C, typename T>
inline unsigned find(const C& c, const T& x) {
  const auto _begin = begin(c);
  const auto _end = end(c);
  const auto it = std::find(_begin,_end,x);
  if (it==_end) throw error(x," not found");
  return std::distance(_begin,it);
}

template <typename T>
std::string type_name() {
  return (std::is_floating_point<T>::value ? "f" :
      (std::is_signed<T>::value ? "i" : "u")
    ) + std::to_string(sizeof(T));
}

int main(int argc, char* argv[]) {
  std::stringstream out;
  auto write = [&out](const auto& x){
    out.write(reinterpret_cast<const char*>(&x),sizeof(x));
  };

  size_type n_events = 0;
  std::vector<unsigned> ph_i(2), jet_i;

  std::ifstream fnames("mxaod.txt");
  for (std::string fname; getline(fnames,fname); ) {
    TEST(fname);
    TFile fin(fname.c_str());

    TTreeReader reader("CollectionTree",&fin);
    branch_reader<Char_t> isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
    branch_reader<UInt_t> runNumber(reader,"EventInfoAux.runNumber");
    branch_reader<ULong64_t> eventNumber(reader,"EventInfoAux.eventNumber");
    std::array<branch_reader<std::vector<float_t>>,4> _photons {{
      {reader,"HGamPhotonsAuxDyn.pt"},
      {reader,"HGamPhotonsAuxDyn.eta"},
      {reader,"HGamPhotonsAuxDyn.phi"},
      {reader,"HGamPhotonsAuxDyn.m"}
    }};
    std::array<branch_reader<std::vector<float_t>>,4> _jets {{
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.pt"},
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.eta"},
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.phi"},
      {reader,"HGamAntiKt4EMTopoJetsAuxDyn.m"}
    }};
    std::array<branch_reader<float_t>,2> _pT_y {{
      {reader,"HGamEventInfoAuxDyn.pT_y1"},
      {reader,"HGamEventInfoAuxDyn.pT_y2"}
    }};

    for ( ivanp::timed_counter<Long64_t> ent(reader.GetEntries(true));
          reader.Next(); ++ent)
    {
      if (!*isPassed) continue;
      ++n_events;

      const auto& photon_pt = *_photons[0];

      auto pT_y1 = *_pT_y[0], pT_y2 = *_pT_y[1];
      if (pT_y1 < pT_y2) std::swap(pT_y1,pT_y2);
      try {
        ph_i = { find(photon_pt,pT_y1), find(photon_pt,pT_y2) };
        if (ph_i[0]==ph_i[1]) throw error("same index");
      } catch (const std::exception& e) {
        cerr << "Photons: " << e << '\n';
        return 1;
      }

      const auto& jet_pt = *_jets[0];

      jet_i.resize(jet_pt.size());
      if (jet_i.size()) {
        std::iota(jet_i.begin(), jet_i.end(), 0);
        std::sort(jet_i.begin(), jet_i.end(), [&](auto a, auto b){
          return jet_pt[a] > jet_pt[b];
        });
        while (jet_i.size() && jet_pt[jet_i.back()] < 30e3)
          jet_i.pop_back();
      }

      write(*runNumber);
      write(*eventNumber);

      for (auto i : ph_i) {
        write(float_t((*_photons[0])[i]*1e-3));
        write((*_photons[1])[i]);
        write((*_photons[2])[i]);
        write(float_t((*_photons[3])[i]*1e-3));
      }

      const size_type njets = jet_i.size();
      write(njets);
      if (njets) for (auto i : jet_i) {
        write(float_t((*_jets[0])[i]*1e-3));
        write((*_jets[1])[i]);
        write((*_jets[2])[i]);
        write(float_t((*_jets[3])[i]*1e-3));
      }
    }
  }
  TEST(n_events)

  std::ofstream("data.dat")
    << R"({"root":[["event#)"
    << n_events
    << R"(","events"]],"types":{"event":[)"
    << "[\"" << type_name<UInt_t>() << "\",\"runNumber\"],"
    << "[\"" << type_name<ULong64_t>() << "\",\"eventNumber\"],"
    << R"(["4vec#2","photons"],["4vec#","jets"]],"4vec":[[")"
    << type_name<float_t>()
    << R"(","pt","eta","phi","m"]]}})"
    << out.rdbuf();
}
