#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <numeric>
#include <vector>
#include <string>

#include <TFile.h>
#include <TH1.h>

#include "ivanp/timed_counter.hh"
#include "ivanp/error.hh"
#include "ivanp/root/branch_reader.hh"
#include "ivanp/math/vec4.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;
using namespace ivanp;

using size_type = uint32_t;
using float_t = float;

template <typename C, typename T>
unsigned find(const C& c, const T& x) {
  const auto _begin = begin(c);
  const auto _end = end(c);
  const auto it = std::find(_begin,_end,x);
  if (it==_end) throw error(x," not found");
  return std::distance(_begin,it);
}

#define POWOFFSET 8

double ratcmp(double x) {
  const bool neg = x < 1;
  if (neg) x = 1./x;
  x = std::log10(x-1) + POWOFFSET;
  if (x < 0) return 0.;
  return neg ? -x : x;
}

TH1D* mkhist(const char* name) {
  TH1D* h = new TH1D(name,name,200,-10,10);
  h->SetXTitle("my / mxaod");
  TAxis* a = h->GetXaxis();
  a->SetNdivisions(10, 0, 0, true);
  a->ChangeLabel(6,-1,-1,-1,-1,-1,"1");
  char str[16];
  for (int i=1; i<6; ++i) {
    const double d = std::pow(10,i*2-POWOFFSET) + 1;
    sprintf(str,"%.7g",d);
    a->ChangeLabel(6+i,-1,-1,-1,-1,-1,str);
    sprintf(str,"1/%.7g",d);
    a->ChangeLabel(6-i,-1,-1,-1,-1,-1,str);
  }
  return h;
}

int main(int argc, char* argv[]) {
  size_type n_events = 0;
  size_type n_events_more_nj = 0;
  size_type n_events_less_nj = 0;

  std::vector<unsigned> ph_i(2), jet_i;

  TFile fout("varcmp_mc.root","recreate");
  if (fout.IsZombie()) return 1;

#define VAR(NAME,...) TH1D* h_##NAME = mkhist(#NAME);

#define ALL_VARS \
  VAR(pT_yy) \
  VAR(m_yy) \
  VAR(pT_y1) \
  VAR(pT_y2) \
  VAR(yAbs_yy) \
  VAR(N_j,Int_t,"N_j_30") \
  VAR(pT_j1,Float_t,"pT_j1_30") \
  VAR(pT_j2,Float_t,"pT_j2_30") \
  VAR(pT_j3,Float_t,"pT_j3_30") \
  VAR(m_jj,Float_t,"m_jj_30") \

  ALL_VARS

  std::ifstream fnames("mxaod_mc.txt");
  for (std::string fname; getline(fnames,fname); ) {
    TEST(fname);
    TFile fin(fname.c_str());

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
    std::array<branch_reader<float_t>,2> _pT_y {{
      {reader,"HGamEventInfoAuxDyn.pT_y1"},
      {reader,"HGamEventInfoAuxDyn.pT_y2"}
    }};

#define VAR_PREF "HGamEventInfoAuxDyn."
#define VAR_PREF_TRUTH "HGamTruthEventInfoAuxDyn."

#define GET_MACRO(_1,_2,_3,NAME,...) NAME

#undef VAR
#define VAR(...) GET_MACRO(__VA_ARGS__, VAR_3, VAR_2, VAR_1)(__VA_ARGS__)
#define VAR_3(NAME,TYPE,ACTUAL) \
    branch_reader<TYPE> _##NAME(reader,VAR_PREF ACTUAL);
#define VAR_2(NAME,TYPE) VAR_3(NAME,TYPE,#NAME)
#define VAR_1(NAME) VAR_2(NAME,Float_t)

    ALL_VARS

    for ( ivanp::timed_counter<Long64_t> ent(reader.GetEntries(true));
          reader.Next(); ++ent)
    {
      if (!*isPassed) continue;
      ++n_events;

      const auto& photon_pt = *_photons[0];

      { auto pT_y1 = *_pT_y[0], pT_y2 = *_pT_y[1];
        if (pT_y1 < pT_y2) std::swap(pT_y1,pT_y2);
        try {
          ph_i = { find(photon_pt,pT_y1), find(photon_pt,pT_y2) };
          if (ph_i[0]==ph_i[1]) throw error("same index");
        } catch (const std::exception& e) {
          cerr << "Photons: " << e << '\n';
          return 1;
        }
      }

      const std::array<vec4<>,2> photons {{
        { (*_photons[0])[ph_i[0]],
          (*_photons[1])[ph_i[0]],
          (*_photons[2])[ph_i[0]],
          (*_photons[3])[ph_i[0]], vec4<>::PtEtaPhiM },
        { (*_photons[0])[ph_i[1]],
          (*_photons[1])[ph_i[1]],
          (*_photons[2])[ph_i[1]],
          (*_photons[3])[ph_i[1]], vec4<>::PtEtaPhiM }
      }};
      const auto yy = photons[0] + photons[1];

      { const auto& jet_pt = *_jets[0];
        jet_i.resize(jet_pt.size());
        if (jet_i.size()) {
          std::iota(jet_i.begin(), jet_i.end(), 0);
          std::sort(jet_i.begin(), jet_i.end(), [&](auto a, auto b){
            return jet_pt[a] > jet_pt[b];
          });
          while (jet_i.size() && jet_pt[jet_i.back()] < 30e3)
            jet_i.pop_back();
        }
      }

      const size_type njets = jet_i.size();
      std::vector<vec4<>> jets(njets);
      for (size_type i=0; i<njets; ++i) {
        jets[i] = {
          (*_jets[0])[jet_i[i]],
          (*_jets[1])[jet_i[i]],
          (*_jets[2])[jet_i[i]],
          (*_jets[3])[jet_i[i]], vec4<>::PtEtaPhiM
        };
      }

      // compare ====================================================

      double pT_yy = yy.pt();
      double m_yy  = yy.m();
      double pT_y1 = photons[0].pt();
      double pT_y2 = photons[1].pt();
      double yAbs_yy = std::abs(yy.rap());

      double N_j = njets;

      if (njets < (unsigned)*_N_j) ++n_events_less_nj; else
      if (njets > (unsigned)*_N_j) ++n_events_more_nj;

#undef VAR
#define VAR(NAME,...) h_##NAME->Fill(ratcmp(NAME / *_##NAME));

#undef ALL_VARS
#define ALL_VARS \
  VAR(pT_yy) \
  VAR(m_yy) \
  VAR(pT_y1) \
  VAR(pT_y2) \
  VAR(yAbs_yy) \
  VAR(N_j) \

      ALL_VARS

      // TEST(njets)

      if (njets < 1) continue; // ===================================

      double pT_j1 = jets[0].pt();

#undef ALL_VARS
#define ALL_VARS \
  VAR(pT_j1) \

      ALL_VARS

      if (njets < 2) continue; // ===================================

      double pT_j2 = jets[1].pt();
      double m_jj = (jets[0]+jets[1]).m();

#undef ALL_VARS
#define ALL_VARS \
  VAR(pT_j2) \
  VAR(m_jj) \

      ALL_VARS

      if (njets < 3) continue; // ===================================

      double pT_j3 = jets[2].pt();

#undef ALL_VARS
#define ALL_VARS \
  VAR(pT_j3) \

      ALL_VARS

    }
  }

  TEST(n_events)
  TEST(n_events_less_nj)
  TEST(n_events_more_nj)

  fout.Write();
}
