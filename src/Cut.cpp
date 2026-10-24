#include <cmath>
#include "../include/Cut.h"

namespace cut {

  //--------------------------------------------
  bool InWindow(double val, const std::vector<double> &w) {
    return val > w[0] && val < w[1];
  }
  //--------------------------------------------
  bool InEllipse(double x, double y, const std::vector<double> &px,
                 const std::vector<double> &py, double nx, double ny) {
    double ax = (x - px[0]) / (px[1] * nx);
    double ay = (y - py[0]) / (py[1] * ny);
    return ax * ax + ay * ay <= 1.;
  }
  //--------------------------------------------
  bool InArea(double x, double y, const std::vector<double> &a) {
    return x > a[0] && x < a[1] && y > a[2] && y < a[3];
  }
  //--------------------------------------------
  std::vector<double> Shrink(const std::vector<double> &a, double mx, double my) {
    return {a[0] + mx, a[1] - mx, a[2] + my, a[3] - my};
  }

}
