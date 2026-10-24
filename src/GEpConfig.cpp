#include <iostream>
#include "../include/Constants.h"
#include "../include/GEpConfig.h"

GEpConfig::GEpConfig(const JSONManager &jmgr) {
  fLabel     = jmgr.GetValueFromKey<std::string>("kinematic");
  fEbeam     = jmgr.GetValueFromKey<double>("ebeam");
  fSBStheta  = jmgr.GetValueFromKey<double>("sbs_theta_deg")  * constant::deg2rad;
  fECALtheta = jmgr.GetValueFromKey<double>("ecal_theta_deg") * constant::deg2rad;
  fECALdist  = jmgr.GetValueFromKey<double>("ecal_dist");
  fECALvoff  = jmgr.GetValueFromKey<double>("ecal_voffset", 0.);
  fECALhoff  = jmgr.GetValueFromKey<double>("ecal_hoffset", 0.);
  fSBSfield  = jmgr.GetValueFromKey<double>("sbs_field", 0.);
}

void GEpConfig::Print() const {
  std::cout << "------ GEp configuration: " << fLabel << " ------" << std::endl
            << " Ebeam (GeV)      : " << fEbeam << std::endl
            << " SBS theta (deg)  : " << fSBStheta / constant::deg2rad << std::endl
            << " ECAL theta (deg) : " << fECALtheta / constant::deg2rad << std::endl
            << " ECAL dist (m)    : " << fECALdist << std::endl
            << " ECAL offsets (m) : v = " << fECALvoff << ", h = " << fECALhoff << std::endl
            << " SBS field        : " << fSBSfield << std::endl
            << "-------------------------------------------" << std::endl;
}
