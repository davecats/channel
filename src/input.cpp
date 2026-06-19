#include "channel/input.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace channel {
namespace {

using EntryMap = std::map<std::string, std::map<std::string, std::string>>;

std::string trim(std::string value) {
  auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
  value.erase(value.begin(), std::find_if(value.begin(), value.end(), [&](unsigned char c) { return !is_space(c); }));
  value.erase(std::find_if(value.rbegin(), value.rend(), [&](unsigned char c) { return !is_space(c); }).base(),
              value.end());
  return value;
}

std::string lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

void strip_comment(std::string& line) {
  const auto bang = line.find('!');
  const auto hash = line.find('#');
  const auto semi = line.find(';');
  auto pos = std::string::npos;
  for (const auto candidate : {bang, hash, semi}) {
    if (candidate != std::string::npos && (pos == std::string::npos || candidate < pos)) pos = candidate;
  }
  if (pos != std::string::npos) line.erase(pos);
}

EntryMap read_entries(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("input file not found: " + path);
  }

  EntryMap entries;
  std::string section;
  std::string line;
  while (std::getline(input, line)) {
    strip_comment(line);
    line = trim(line);
    if (line.empty()) continue;
    if (line.front() == '[' && line.back() == ']') {
      section = lower(trim(line.substr(1, line.size() - 2)));
      continue;
    }
    const auto equals = line.find('=');
    if (equals == std::string::npos) continue;
    const auto key = lower(trim(line.substr(0, equals)));
    const auto value = trim(line.substr(equals + 1));
    if (!section.empty() && !key.empty()) entries[section][key] = value;
  }
  return entries;
}

std::string require_value(const EntryMap& entries, const std::string& section, const std::string& key) {
  const auto section_it = entries.find(section);
  if (section_it == entries.end()) {
    throw std::runtime_error("missing input section [" + section + "]");
  }
  const auto value_it = section_it->second.find(key);
  if (value_it == section_it->second.end()) {
    throw std::runtime_error("missing input key " + section + "." + key);
  }
  return value_it->second;
}

bool optional_value(const EntryMap& entries, const std::string& section, const std::string& key, std::string& value) {
  const auto section_it = entries.find(section);
  if (section_it == entries.end()) return false;
  const auto value_it = section_it->second.find(key);
  if (value_it == section_it->second.end()) return false;
  value = value_it->second;
  return true;
}

int parse_int(const std::string& value, const std::string& name) {
  std::istringstream stream(value);
  int parsed = 0;
  stream >> parsed;
  if (stream.fail()) {
    throw std::runtime_error("invalid integer value for " + name + ": " + value);
  }
  stream >> std::ws;
  if (!stream.eof()) {
    throw std::runtime_error("invalid integer value for " + name + ": " + value);
  }
  return parsed;
}

double parse_double(const std::string& value, const std::string& name) {
  std::istringstream stream(value);
  double parsed = 0.0;
  stream >> parsed;
  if (stream.fail()) {
    throw std::runtime_error("invalid real value for " + name + ": " + value);
  }
  stream >> std::ws;
  if (!stream.eof()) {
    throw std::runtime_error("invalid real value for " + name + ": " + value);
  }
  return parsed;
}

std::vector<double> parse_double_vector(const std::string& value, const std::string& name, int count) {
  std::istringstream stream(value);
  std::vector<double> parsed(static_cast<std::size_t>(count), 0.0);
  for (int i = 0; i < count; ++i) {
    stream >> parsed[static_cast<std::size_t>(i)];
    if (!stream) {
      throw std::runtime_error("invalid real vector value for " + name + ": " + value);
    }
  }
  return parsed;
}

int require_int(const EntryMap& entries, const std::string& section, const std::string& key) {
  return parse_int(require_value(entries, section, key), section + "." + key);
}

double require_double(const EntryMap& entries, const std::string& section, const std::string& key) {
  return parse_double(require_value(entries, section, key), section + "." + key);
}

} // namespace

ChannelInput read_channel_input(const std::string& path) {
  const auto entries = read_entries(path);
  ChannelInput cfg;
  cfg.mesh.nx = require_int(entries, "mesh", "nx");
  cfg.mesh.ny = require_int(entries, "mesh", "ny");
  cfg.mesh.nz = require_int(entries, "mesh", "nz");
  cfg.mesh.alfa0 = require_double(entries, "mesh", "alfa0");
  cfg.mesh.beta0 = require_double(entries, "mesh", "beta0");
  cfg.mesh.stretching = require_double(entries, "mesh", "stretching");
  cfg.mesh.ymin = require_double(entries, "mesh", "ymin");
  cfg.mesh.ymax = require_double(entries, "mesh", "ymax");

  cfg.velocity.reynolds = require_double(entries, "velocity", "ni");
  if (cfg.velocity.reynolds == 0.0) {
    throw std::runtime_error("velocity.ni must be nonzero");
  }
  cfg.velocity.viscosity = 1.0 / cfg.velocity.reynolds;
  cfg.velocity.meanpx = require_double(entries, "velocity", "meanpx");
  cfg.velocity.meanpz = require_double(entries, "velocity", "meanpz");
  cfg.velocity.meanflowx = require_double(entries, "velocity", "meanflowx");
  cfg.velocity.meanflowz = require_double(entries, "velocity", "meanflowz");
  cfg.velocity.u0 = require_double(entries, "velocity", "u0");
  cfg.velocity.uN = require_double(entries, "velocity", "un");

  cfg.scalars.nphi = require_int(entries, "scalars", "nphi");
  if (cfg.scalars.nphi < 0) {
    throw std::runtime_error("scalars.nPhi must be nonnegative");
  }
  if (cfg.scalars.nphi > 0) {
    const auto pr = require_value(entries, "scalars", "pr");
    cfg.scalars.inverse_prandtl = parse_double_vector(pr, "scalars.pr", cfg.scalars.nphi);
    for (auto& value : cfg.scalars.inverse_prandtl) {
      if (value == 0.0) throw std::runtime_error("scalars.pr entries must be nonzero");
      value = 1.0 / value;
    }
  }

  cfg.timestepping.dt = require_double(entries, "timestepping", "deltat");
  cfg.timestepping.cflmax = require_double(entries, "timestepping", "cflmax");
  cfg.timestepping.time = require_double(entries, "timestepping", "time");
  cfg.timestepping.t_max = require_double(entries, "timestepping", "t_max");
  cfg.timestepping.nstep = require_int(entries, "timestepping", "nstep");

  std::string npy;
  if (optional_value(entries, "parallel", "npy", npy) || optional_value(entries, "mesh", "npy", npy)) {
    cfg.parallel.npy = parse_int(npy, "parallel.npy");
    cfg.parallel.npy_was_set = true;
  }
  if (cfg.parallel.npy < 1) {
    throw std::runtime_error("parallel.npy must be positive");
  }
  return cfg;
}

} // namespace channel
