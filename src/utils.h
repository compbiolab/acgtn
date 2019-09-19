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

#include "types.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#define FreePtr(p) \
  if (p) {         \
    delete p;      \
    p = nullptr;   \
  }

#define FreeArrPtr(p) \
  if (p) {            \
    delete[] p;       \
    p = nullptr;      \
  }

#define MIN(x, y) ((x) <= (y) ? (x) : (y))
#define MAX(x, y) ((x) >= (y) ? (x) : (y))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define BITMASK(n) (((n) == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << (n)) - 1))

#define base2rev(c) \
  ((c) < 'G' ? ((c) == 'A' ? 'T' : 'G') : ((c) == 'G' ? 'C' : 'A'))
#define base2int(c) \
  ((c) < 'G' ? ((c) == 'A' ? 0ULL : 1ULL) : ((c) == 'G' ? 2ULL : 3ULL))
#define int2base(c) (1413956417ULL >> (8 * (c)) & 0xFF)

#ifndef BITS_INDEX_32
#define uintx uint40
#else
#define uintx uint32
#endif

const uintx MaxUintx = std::numeric_limits<uintx>::max();

using Range = std::pair<uint64, uint64>;

extern const uint64 BitMask[65];

constexpr uint64 KILOBYTE = 1024;
constexpr uint64 MEGABYTE = KILOBYTE * KILOBYTE;
constexpr uint64 GIGABYTE = KILOBYTE * MEGABYTE;

inline bool is_even(uint64 x) { return (0 == (x & 1UL)); }
inline bool is_odd(uint64 x) { return (0 != (x & 1UL)); }

// Return floor((x+y)/2)
inline uint64 average(uint64 x, uint64 y) { return (x & y) + ((x ^ y) >> 1); }

// Return ceil((x+y)/2)
inline uint64 ceil_average(uint64 x, uint64 y) {
  return (x | y) - ((x ^ y) >> 1);
}

// Return random number in the range [0, 1, ..., m).
// Must have m > 0.
inline uint64 rand_idx(uint64 m) {
  if (m == 1) return 0;  // could also use % 1
  uint64 x = (uint64)rand();
  x ^= x >> 16;  // avoid using low bits of rand() alone
  return x % m;
}

template <typename Type>
inline Type *resize(Type *p, uint64 n) {
  return (Type *)realloc((void *)p, n * sizeof(Type));
}

template <typename Type>
inline void memzeros(Type *dst, uint64 n) {
  memset(dst, 0, n * sizeof(Type));
}

template <typename Type>
inline void memcopy(const Type *src, Type *dst, uint64 n) {
  memcpy(dst, src, n * sizeof(Type));
}

template <typename Type>
inline void read_member(Type &t, std::ifstream &in) {
  in.read((char *)&t, sizeof(t));
}

template <typename Type>
inline uint64 write_member(const Type &t, std::ofstream &out) {
  out.write((char *)&t, sizeof(t));
  return sizeof(t);
}

int remove_dir(const char *dir);

std::vector<Range> jobBounds(uint64 num, uint64 thds);
