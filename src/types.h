/*
  Copyright 2019, Zihua Wu <wuzihua@picb.ac.cn>

  This file is part of ACGTN.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

typedef std::int8_t int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;

typedef std::uint8_t uint8, uchar, byte;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

typedef std::vector<uint32> Vuint32;
typedef std::vector<uint64> Vuint64;

typedef std::vector<std::string> Vstring;

class uint40 {
 public:
  inline uint40() {}

  inline uint40(uint32 l, uint8 h) : low(l), high(h) {}

  inline uint40(const uint40& a) : low(a.low), high(a.high) {}

  inline uint40(const uint64& a) : low(a & 0xFFFFFFFF), high((a >> 32) & 0xFF) {
    assert(a <= 0xFFFFFFFFFFULL);
  }

  inline uint64 ull() const { return ((uint64)high) << 32 | (uint64)low; }

  inline operator uint64() const { return ull(); }

  inline uint64 u64() const { return ((uint64)high) << 32 | (uint64)low; }

  inline uint40& operator++() {
    if (low == std::numeric_limits<uint32>::max())
      ++high, low = 0;
    else
      ++low;
    return *this;
  }

  inline uint40& operator--() {
    if (low == 0)
      --high, low = std::numeric_limits<uint32>::max();
    else
      --low;
    return *this;
  }

  inline uint40& operator+=(const uint40& b) {
    uint64 add = low + b.low;
    low = add & 0xFFFFFFFF;
    high += b.high + ((add >> 32) & 0xFF);
    return *this;
  }

  inline bool operator==(const uint40& b) const {
    return (low == b.low) && (high == b.high);
  }

  inline bool operator!=(const uint40& b) const {
    return (low != b.low) || (high != b.high);
  }

  inline bool operator<(const uint40& b) const {
    return (high < b.high) || (high == b.high && low < b.low);
  }

  inline bool operator<=(const uint40& b) const {
    return (high < b.high) || (high == b.high && low <= b.low);
  }

  inline bool operator>(const uint40& b) const {
    return (high > b.high) || (high == b.high && low > b.low);
  }

  inline bool operator>=(const uint40& b) const {
    return (high > b.high) || (high == b.high && low >= b.low);
  }

  friend std::ostream& operator<<(std::ostream& os, const uint40& a) {
    return os << a.ull();
  }

 private:
  uint32 low;
  uint8 high;

} __attribute__((packed));

namespace std {

template <>
class numeric_limits<uint40> {
 public:
  static uint40 min() {
    return uint40(std::numeric_limits<uint32_t>::min(),
                  std::numeric_limits<uint8_t>::min());
  }

  static uint40 max() {
    return uint40(std::numeric_limits<uint32_t>::max(),
                  std::numeric_limits<uint8_t>::max());
  }
};

}  // namespace std
