#ifndef GEP_CONFIG_H
#define GEP_CONFIG_H

#include <string>
#include "JSONManager.h"

/*
   GEpConfig: everything that defines one GEp kinematic setting.
   Unlike gmn-ana's SBSconfig (values hard-coded per config number), this is
   read from the JSON config so a new kinematic only needs a new config file.
   Beam energy is the nominal value from the config (no EPICS averaging here).
*/
class GEpConfig {
public:
  GEpConfig(const JSONManager &jmgr);

  std::string GetLabel()       const { return fLabel; }
  double GetEbeam()            const { return fEbeam; }        // GeV
  double GetSBStheta_rad()     const { return fSBStheta; }     // proton arm (beam right)
  double GetECALtheta_rad()    const { return fECALtheta; }    // electron arm (beam left)
  double GetECALdist()         const { return fECALdist; }     // m, target -> ECAL face
  double GetECALvoffset()      const { return fECALvoff; }     // m, along ECAL x (vertical)
  double GetECALhoffset()      const { return fECALhoff; }     // m, along ECAL y (horizontal)
  double GetSBSfield()         const { return fSBSfield; }     // % or T, whatever you put in config

  void Print() const;

private:
  std::string fLabel;
  double fEbeam, fSBStheta, fECALtheta, fECALdist, fECALvoff, fECALhoff, fSBSfield;
};

#endif
