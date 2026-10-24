#ifndef GEP_JSON_MANAGER_H
#define GEP_JSON_MANAGER_H

/*
   Minimal header-only JSON config reader (same interface as the JSONManager
   gmn-ana uses from D. Flay's SBS-ANA, but without the hard-coded ifarm path).
   Needs nlohmann/json.hpp, which ships with ROOT >= 6.24 (and conda-forge ROOT).
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

class JSONManager {
public:
  JSONManager(const std::string &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
      std::cerr << "[JSONManager] Cannot open config file " << filename << " **!**" << std::endl;
      throw std::runtime_error("JSONManager: missing config");
    }
    ifs >> fObj;
  }

  bool HasKey(const std::string &key) const { return fObj.contains(key); }

  template <typename T>
  T GetValueFromKey(const std::string &key) const {
    if (!HasKey(key)) {
      std::cerr << "[JSONManager] Key '" << key << "' not found in config **!**" << std::endl;
      throw std::runtime_error("JSONManager: missing key " + key);
    }
    return fObj.at(key).get<T>();
  }

  // returns default value if the key is absent
  template <typename T>
  T GetValueFromKey(const std::string &key, T def) const {
    return HasKey(key) ? fObj.at(key).get<T>() : def;
  }

  std::string GetValueFromKey_str(const std::string &key) const {
    return GetValueFromKey<std::string>(key);
  }

  template <typename T>
  void GetVectorFromKey(const std::string &key, std::vector<T> &vec) const {
    vec = GetValueFromKey<std::vector<T>>(key);
  }

private:
  nlohmann::json fObj;
};

#endif
