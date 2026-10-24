#include <algorithm>
#include <iostream>
#include "../include/SetROOTVar.h"

namespace setrootvar {

  //--------------------------------------------
  static int attach(TTree *T, const std::string &branchname, void *memory) {
    if (!T->GetBranch(branchname.c_str())) {
      std::cerr << "[setrootvar] WARNING: branch '" << branchname
                << "' not found, skipping (variable keeps its default value)." << std::endl;
      return 0;
    }
    T->SetBranchStatus(branchname.c_str(), 1);
    T->SetBranchAddress(branchname.c_str(), memory);
    return 1;
  }
  //--------------------------------------------
  int setbranch(TTree *T, std::string prefix, std::string suffix, void *memory) {
    T->SetMakeClass(1);   // allows access to individual sub-branches
    std::string name = suffix.empty() ? prefix : prefix + "." + suffix;
    return attach(T, name, memory);
  }
  //--------------------------------------------
  int setbranch(TTree *T, std::string prefix,
                std::vector<std::string> suffix, std::vector<void*> memory) {
    return setbranch(T, prefix, suffix, memory, std::vector<int>{});
  }
  //--------------------------------------------
  int setbranch(TTree *T, std::string prefix,
                std::vector<std::string> suffix, std::vector<void*> memory,
                std::vector<int> ndatapos) {
    if (suffix.size() != memory.size()) {
      std::cerr << "[setrootvar] suffix and memory containers differ in size **!**" << std::endl;
      throw std::runtime_error("setrootvar: size mismatch");
    }
    T->SetMakeClass(1);
    int nok = 0;
    for (size_t i = 0; i < suffix.size(); i++) {
      bool isNdata = std::find(ndatapos.begin(), ndatapos.end(), int(i)) != ndatapos.end();
      std::string name = (isNdata ? "Ndata." : "") + prefix + "." + suffix[i];
      nok += attach(T, name, memory[i]);
    }
    return nok;
  }

}
