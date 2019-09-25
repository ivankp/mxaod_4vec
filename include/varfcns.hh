#ifndef VARFCNS_HH
#define VARFCNS_HH

#include <map>
#include <cmath>

struct less_str {
  bool operator()(const char* a, const char* b) const noexcept {
    return strcmp(a,b) < 0;
  }
};

// ==================================================================
bool nj(unsigned n) noexcept { return nj() >= n; }

double f_HT_jets() noexcept {
  double HT = 0;
  for (const auto& jet : jets) HT += jet.pt();
  return HT;
}
double f_HT_jets_yy() noexcept { return f_HT_jets() + yy.pt(); }
// ==================================================================

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
  { "pT_j1", { []{ return nj(1) ? jets[0].pt() : 0; }, true } },
  { "pT_j2", { []{ return nj(2) ? jets[1].pt() : 0; }, true } },
  { "pT_j3", { []{ return nj(3) ? jets[2].pt() : 0; }, true } },
  { "eta_j1", { []{ return nj(1) ? jets[0].eta() : NAN; }, true } },
  { "eta_j2", { []{ return nj(2) ? jets[1].eta() : NAN; }, true } },
  { "eta_j3", { []{ return nj(3) ? jets[2].eta() : NAN; }, true } },
  { "y_j1", { []{ return nj(1) ? jets[0].rap() : NAN; }, true } },
  { "y_j2", { []{ return nj(2) ? jets[1].rap() : NAN; }, true } },
  { "y_j3", { []{ return nj(3) ? jets[2].rap() : NAN; }, true } },
  { "dy_j1_j2", {
    []{ return nj(2) ? std::abs(jets[0].rap()-jets[1].rap()) : NAN; }, true } },
  { "m_jj", { []{ return nj(2) ? (jets[0]+jets[1]).m() : NAN; }, true } },

  { "Hj_mass", { []{ return nj(1) ? (yy+jets[0]).m() : NAN; }, true } },

  { "HT_jets", { f_HT_jets, true } },
  { "HT_jets_yy", { f_HT_jets_yy, true } },

  { "x_yy", { []{ return yy.pt()/f_HT_jets_yy(); }, true } },
  { "x_j1", { []{ return nj(1) ? jets[0].pt()/f_HT_jets_yy() : 0; }, true } },
  { "x_j2", { []{ return nj(2) ? jets[1].pt()/f_HT_jets_yy() : 0; }, true } },
  { "x_j3", { []{ return nj(3) ? jets[2].pt()/f_HT_jets_yy() : 0; }, true } },

  { "pT_miss", { []{
      auto p = yy;
      for (const auto& jet : jets) p += jet;
      return p.pt();
    }, true } },
  { "s2", { []{
      auto p = yy;
      for (const auto& jet : jets) p += jet;
      return p.m();
    }, true } },
};

#endif
