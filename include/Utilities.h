#ifndef GEP_UTILITIES_H
#define GEP_UTILITIES_H

#include <string>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"

namespace util {

  // TH1F from {nbins, min, max} (the "*_lims" arrays in the config)
  TH1F *TH1Flims(const std::string &name, const std::string &title, const std::vector<double> &lims);
  TH2F *TH2Flims(const std::string &name, const std::string &title,
                 const std::vector<double> &xlims, const std::vector<double> &ylims);

  // draw a rectangle {xmin, xmax, ymin, ymax}; detector x is vertical so it is drawn on the y axis
  void DrawArea(const std::vector<double> &area, int lcolor = 2, int lwidth = 2, int lstyle = 9);

}

#endif
