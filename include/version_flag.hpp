/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <args.hxx>
#include <string>
#include <utility>
#include <vector>


#pragma once

namespace sf::args {

class output_version : public ::args::Error {
 public:
  output_version(std::string const &flag) : Error(flag) {}
  virtual ~output_version() {}
};

class version_flag : public ::args::Flag {
 public:
  version_flag(::args::Group &group_, std::string const &name_,
               std::string const &help_, ::args::Matcher &&matcher_,
               ::args::Options options_ = {})
      : Flag(group_, name_, help_, std::move(matcher_), options_) {}

  virtual ~version_flag() {}

  virtual void ParseValue(std::vector<std::string> const &) {
    throw output_version(Name());
  }

  bool Get() const noexcept { return Matched(); }
};

}  // namespace sf::args
