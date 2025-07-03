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


// Taywee/args uses plain names (i.e. no e.g. trailing underscore) for
// non-static data members.
template <typename t_value, typename t_reader = ::args::ValueReader>
class value_flag : public ::args::ValueFlagBase {
 protected:
  t_value &value;
  t_value default_value{};

 private:
  t_reader reader;

 protected:
  virtual std::string GetDefaultString(
      ::args::HelpParams const &) const override {
    return ::args::detail::ToString(default_value);
  }

 public:
  value_flag(::args::Group &group_, std::string const &name_,
             std::string const &help_, ::args::Matcher &&matcher_,
             t_value &value_, t_value const &default_value_,
             ::args::Options options_)
      : ::args::ValueFlagBase(name_, help_, std::move(matcher_), options_),
        value(value_),
        default_value(default_value_) {
    group_.Add(*this);
  }

  value_flag(::args::Group &group_, std::string const &name_,
             std::string const &help_, ::args::Matcher &&matcher_,
             t_value &value_, ::args::Options options_)
      : value_flag(group_, name_, help_, std::move(matcher_), value_, value_,
                   options_) {}

  value_flag(::args::Group &group_, std::string const &name_,
             std::string const &help_, ::args::Matcher &&matcher_,
             t_value &value_)
      : value_flag(group_, name_, help_, std::move(matcher_), value_, value_,
                   ::args::Options::None) {}

  virtual ~value_flag() {}
  virtual void ParseValue(std::vector<std::string> const &values) override;
  virtual void Reset() noexcept override;
  t_value &Get() noexcept { return value; }
  t_value &operator*() noexcept { return value; }
  t_value const &operator*() const noexcept { return value; }
  t_value *operator->() noexcept { return &value; }
  t_value const *operator->() const noexcept { return &value; }
  t_value const &GetDefault() const noexcept { return default_value; }
};


template <typename t_value, typename t_reader>
void value_flag<t_value, t_reader>::Reset() noexcept {
	ValueFlagBase::Reset();
	value = default_value;
}


template <typename t_value, typename t_reader>
void value_flag<t_value, t_reader>::ParseValue(
    std::vector<std::string> const &values) {
  auto const &value_ = values.at(0);
  reader(name, value_, value);
}

}  // namespace sf::args
