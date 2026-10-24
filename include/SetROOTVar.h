#ifndef GEP_SET_ROOT_VAR_H
#define GEP_SET_ROOT_VAR_H

#include <string>
#include <vector>
#include "TTree.h"

/*
   Branch-address helpers, same interface as gmn-ana's setrootvar.
   Difference: a missing branch prints a warning and is skipped instead of
   crashing, since GEp replay branch names are still being settled.
   All return the number of branches successfully attached.
*/
namespace setrootvar {

  // one variable: prefix.suffix (or just prefix if suffix is empty)
  int setbranch(TTree *T, std::string prefix, std::string suffix, void *memory);

  // several variables sharing one prefix
  int setbranch(TTree *T, std::string prefix,
                std::vector<std::string> suffix, std::vector<void*> memory);

  // same, but the entries at positions ndatapos are Ndata.prefix.suffix counters
  int setbranch(TTree *T, std::string prefix,
                std::vector<std::string> suffix, std::vector<void*> memory,
                std::vector<int> ndatapos);

}

#endif
