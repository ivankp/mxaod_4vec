#ifndef VARFCNS_HH
#define VARFCNS_HH

#include <map>

struct less_str {
  bool operator()(const char* a, const char* b) const noexcept {
    return strcmp(a,b) < 0;
  }
};

// ==================================================================
bool nj(unsigned n) noexcept { return nj() >= n; }

double f_HT_jets() noexcept {
  double HT = 0;
  for (const auto& j : j) HT += j.pt();
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
  { "pT_j1", { []{ return nj(1) ? j[0].pt() : 0; }, true } },
  { "pT_j2", { []{ return nj(2) ? j[1].pt() : 0; }, true } },
  { "pT_j3", { []{ return nj(3) ? j[2].pt() : 0; }, true } },
  { "eta_j1", { []{ return nj(1) ? j[0].eta() : 0; }, true } },
  { "eta_j2", { []{ return nj(2) ? j[1].eta() : 0; }, true } },
  { "eta_j3", { []{ return nj(3) ? j[2].eta() : 0; }, true } },
  { "y_j1", { []{ return nj(1) ? j[0].rap() : 0; }, true } },
  { "y_j2", { []{ return nj(2) ? j[1].rap() : 0; }, true } },
  { "y_j3", { []{ return nj(3) ? j[2].rap() : 0; }, true } },
  { "dy_j1_j2", { []{ return nj(2) ? std::abs(j[0].rap()-j[1].rap()) : 0; }, true } },
  { "m_jj", { []{ return nj(2) ? (j[0]+j[1]).m() : 0; }, true } },

  { "Hj_mass", { []{ return nj(1) ? (yy+j[0]).m() : 0; }, true } },

  { "HT_jets", { f_HT_jets, true } },
  { "HT_jets_yy", { f_HT_jets_yy, true } },

  { "x_yy", { []{ return yy.pt()/f_HT_jets_yy(); }, true } },
  { "x_j1", { []{ return nj(1) ? j[0].pt()/f_HT_jets_yy() : 0; }, true } },
  { "x_j2", { []{ return nj(2) ? j[1].pt()/f_HT_jets_yy() : 0; }, true } },
  { "x_j3", { []{ return nj(3) ? j[2].pt()/f_HT_jets_yy() : 0; }, true } },

  { "pT_miss", { []{
      auto p = yy;
      for (const auto& j : j) p += j;
      return p.pt();
    }, true } },
  { "s2", { []{
      auto p = yy;
      for (const auto& j : j) p += j;
      return p.m();
    }, true } },
};

#endif
