#include "TLine.h"
#include "../include/Utilities.h"

namespace util {

  //--------------------------------------------
  TH1F *TH1Flims(const std::string &name, const std::string &title, const std::vector<double> &l) {
    return new TH1F(name.c_str(), title.c_str(), int(l[0]), l[1], l[2]);
  }
  //--------------------------------------------
  TH2F *TH2Flims(const std::string &name, const std::string &title,
                 const std::vector<double> &x, const std::vector<double> &y) {
    return new TH2F(name.c_str(), title.c_str(), int(x[0]), x[1], x[2], int(y[0]), y[1], y[2]);
  }
  //--------------------------------------------
  void DrawArea(const std::vector<double> &a, int lcolor, int lwidth, int lstyle) {
    // plotted as (y, x): horizontal axis = detector y, vertical axis = detector x
    double xmin = a[0], xmax = a[1], ymin = a[2], ymax = a[3];
    TLine *l[4] = {new TLine(ymin, xmin, ymax, xmin), new TLine(ymin, xmax, ymax, xmax),
                   new TLine(ymin, xmin, ymin, xmax), new TLine(ymax, xmin, ymax, xmax)};
    for (auto *li : l) {
      li->SetLineColor(lcolor); li->SetLineWidth(lwidth); li->SetLineStyle(lstyle); li->Draw();
    }
  }

}
