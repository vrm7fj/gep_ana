#ifndef GEP_CUT_H
#define GEP_CUT_H

#include <vector>

namespace cut {

  // lo < val < hi ; window = {lo, hi}
  bool InWindow(double val, const std::vector<double> &window);

  // elliptical cut: ((x-mux)/(n*sx))^2 + ((y-muy)/(n*sy))^2 <= 1
  // peak = {mean, sigma}
  bool InEllipse(double x, double y, const std::vector<double> &peak_x,
                 const std::vector<double> &peak_y, double nsig_x, double nsig_y);

  // rectangular area = {xmin, xmax, ymin, ymax} in detector coordinates
  bool InArea(double x, double y, const std::vector<double> &area);
  // shrink an area by margins (x and y) on all sides -> fiducial region
  std::vector<double> Shrink(const std::vector<double> &area, double mx, double my);

}

#endif
