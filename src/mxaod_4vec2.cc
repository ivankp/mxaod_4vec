#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <string>
#include <memory>

#include <nlohmann/json.hpp>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>

#include "ivanp/timed_counter.hh"
#include "ivanp/error.hh"
#include "ivanp/pcre_wrapper.hh"
#include "ivanp/root/branch_reader.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;
using std::string;
using std::vector;
using ivanp::error;
using ivanp::cat;

using float_t = float;

template <typename C, typename T>
unsigned find(const C& c, const T& x) {
  const auto _begin = begin(c);
  const auto _end = end(c);
  const auto it = std::find(_begin,_end,x);
  if (it==_end) throw error(x," not found");
  return std::distance(_begin,it);
}

template <typename T, typename... Args>
void make(std::unique_ptr<T>& p, Args&&... args) {
  p = std::make_unique<T>(std::forward<Args>(args)...);
}

struct set {
  string name;
  vector<string> data, mc;
  double lumi = 0;
};

int main(int argc, char* argv[]) {
  vector<set> sets;
  { nlohmann::json cfg;
    std::ifstream("mxaods2.json") >> cfg;
    const string dir = cfg["dir"];
    for (const auto& s : cfg["sets"].items()) {
      sets.emplace_back();
      auto& _s = sets.back();
      _s.name = s.key();
      for (const auto& xs : s.value().items()) {
        const bool is_data = xs.key()=="data";
        auto& _xs = is_data ? _s.data : _s.mc;
        for (const auto& xs : xs.value()) {
          const string subdir = xs[0];
          for (const string name : xs[1]) {
            _xs.emplace_back(dir+subdir+name);
            if (is_data) {
              static ivanp::pcre::regex
                re("(?:.*/)?data.*_(\\d+)ipb\\..*\\.root$");
              ivanp::pcre::match m;
              if (re(m,name)) {
                _s.lumi += atoi(m[1][0]);
              } else throw error(
                "filename \"",name,"\" didn\'t match lumi regex"
              );
            }
          }
        }
      }
    }
  }

  double total_lumi = 0;
  for (const auto& s : sets) total_lumi += s.lumi;
  TEST(total_lumi)

  /*
  for (const auto& s : sets) {
    cout << s.name << " " << s.lumi << endl;
    cout << "data" << endl;
    for (const auto& x : s.data)
      cout << "  " << x << endl;
    cout << "mc" << endl;
    for (const auto& x : s.mc)
      cout << "  " << x << endl;
    cout << endl;
  }
  */

  double mc_factor = 0;
  std::vector<unsigned> ph_i(2), jet_i;

  for (const bool is_mc : {false,true}) {
    uint32_t nevents = 0;

    std::ofstream out(cat(
      ( argc>1 && strlen(argv[1]) ? argv[1] : "." ),
      "/hgam_", (is_mc?"mc":"data"), ".dat"
    ));
    auto write = [&out](const auto& x){
      out.write(reinterpret_cast<const char*>(&x),sizeof(x));
    };

    if (is_mc) {
      out << 'm';
    } else {
      out << 'd';
      write(float_t(total_lumi*1e-3));
    }

    const auto nevents_pos = out.tellp();
    write(nevents);

    for (const auto& s : sets) {
    for (const auto& fname : (is_mc ? s.mc : s.data)) {
      TEST(fname);
      TFile fin(fname.c_str());

      if (is_mc) { // MC
        for (auto* key : *fin.GetListOfKeys()) {
          const char* name = key->GetName();
          if (!ivanp::starts_with(name,"CutFlow_") ||
              !ivanp::ends_with(name,"_noDalitz_weighted")) continue;
          TH1 *h = static_cast<TH1*>(static_cast<TKey*>(key)->ReadObj());
          cout << name << endl;
          const double n_all = h->GetBinContent(3);
          cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << endl;
          mc_factor = 1e3*(s.lumi/total_lumi)/n_all;
          break;
        }
      }

      TTreeReader reader("CollectionTree",&fin);
      branch_reader<Char_t> isPassed(reader,"HGamEventInfoAuxDyn.isPassed");
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

      using float_branch = branch_reader<float_t>;

      std::array<float_branch,2> _pT_y {{
        {reader,"HGamEventInfoAuxDyn.pT_y1"},
        {reader,"HGamEventInfoAuxDyn.pT_y2"}
      }};
      float_branch _m_yy(reader,"HGamEventInfoAuxDyn.m_yy");

      // MC
      std::unique_ptr<float_branch> _weight, _cs_br_fe;
      if (is_mc) {
        make(_weight,reader,"HGamEventInfoAuxDyn.weight");
        make(_cs_br_fe,reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      }

      for ( ivanp::timed_counter<Long64_t> ent(reader.GetEntries(true));
            reader.Next(); ++ent)
      {
        if (!*isPassed) continue;

        // diphoton mass cut
        const double m_yy = *_m_yy*1e-3;
        if (m_yy<105. || 160.<m_yy) continue;

        ++nevents; // number of events after cuts

        if (is_mc) {
          write( float_t( // weight
            double(**_weight) * double(**_cs_br_fe) * mc_factor
          ));
        }

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
          // jet pT cut
          while (jet_i.size() && jet_pt[jet_i.back()] < 30e3)
            jet_i.pop_back();
          // jet eta cut
          jet_i.erase(std::remove_if( jet_i.begin(), jet_i.end(),
            [&](const auto& i) { return std::abs((*_jets[1])[i]) > 4.4; }),
            jet_i.end());
        }

        for (auto i : ph_i) {
          write(float_t((*_photons[0])[i]*1e-3));
          write((*_photons[1])[i]);
          write((*_photons[2])[i]);
          write(float_t((*_photons[3])[i]*1e-3));
        }

        const uint8_t njets = jet_i.size();
        write(njets);
        if (njets > 4) jet_i.resize(4);
        for (auto i : jet_i) {
          write(float_t((*_jets[0])[i]*1e-3));
          write((*_jets[1])[i]);
          write((*_jets[2])[i]);
          write(float_t((*_jets[3])[i]*1e-3));
        }
      }
    }}
    TEST(nevents)

    out.flush();
    out.seekp(nevents_pos);
    write(nevents);
    out.flush();
  }
}
