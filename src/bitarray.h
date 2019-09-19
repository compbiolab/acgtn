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

#include "bits.h"

#include <assert.h>
#include <cstring>
#include <fstream>

#ifndef BITS_USE_ASM
#define DIVMOD(n, d, bm)        \
  uint64 d = n / BITS_PER_WORD; \
  uint64 bm = 1UL << (n % BITS_PER_WORD);

#define DIVMOD_TEST(n, d, bm)             \
  uint64 d = n / BITS_PER_WORD;           \
  uint64 bm = 1UL << (n % BITS_PER_WORD); \
  uint64 t = bm & f_[d];
#endif

class bitarray {
 public:
  explicit bitarray(uint64 nbits, uint64 *f = nullptr) {
    uint64 nw = ctor_core(nbits);
    if (f) {
      f_ = f;
      myfq_ = false;
    } else {
      f_ = new uint64[nw]();
      myfq_ = true;
    }
  }

  ~bitarray() {
    if (myfq_) delete[] f_;
  }

  bitarray(const bitarray &) = delete;
  bitarray &operator=(const bitarray &) = delete;

  uint64 size() const { return n_; }

  // Test whether n-th bit set.
  uint64 test(uint64 n) const {
    assert(n < n_);
#ifdef BITS_USE_ASM
    return asm_bt(f_, n);
#else
    DIVMOD_TEST(n, d, bm);
    return t;
#endif
  }

  // Set n-th bit.
  void set(uint64 n) {
    assert(n < n_);
#ifdef BITS_USE_ASM
    asm_b_s(f_, n);
#else
    DIVMOD(n, d, bm);
    f_[d] |= bm;
#endif
  }

  // Clear n-th bit.
  void clear(uint64 n) {
    assert(n < n_);
#ifdef BITS_USE_ASM
    asm_b_r(f_, n);
#else
    DIVMOD(n, d, bm);
    f_[d] &= ~bm;
#endif
  }

  // Toggle n-th bit.
  void change(uint64 n) {
    assert(n < n_);
#ifdef BITS_USE_ASM
    asm_b_c(f_, n);
#else
    DIVMOD(n, d, bm);
    f_[d] ^= bm;
#endif
  }

  // Test whether n-th bit is set and set it.
  uint64 test_set(uint64 n) {
    assert(n < n_);
#ifdef BITS_USE_ASM
    return asm_bts(f_, n);
#else
    DIVMOD_TEST(n, d, bm);
    f_[d] |= bm;
    return t;
#endif
  }

  // Test whether n-th bit is set and clear it.
  uint64 test_clear(uint64 n) {
    assert(n < n_);
#ifdef BITS_USE_ASM
    return asm_btr(f_, n);
#else
    DIVMOD_TEST(n, d, bm);
    f_[d] &= ~bm;
    return t;
#endif
  }

  // Test whether n-th bit is set and toggle it.
  uint64 test_change(uint64 n) {
    assert(n < n_);
#ifdef BITS_USE_ASM
    return asm_btc(f_, n);
#else
    DIVMOD_TEST(n, d, bm);
    f_[d] ^= bm;
    return t;
#endif
  }

  // Clear all bits.
  void clear_all() {
    for (uint64 k = 0; k < nfw_; ++k) f_[k] = 0;
    if (mp_) f_[nfw_] = 0;
  }

  // Set all bits.
  void set_all() {
    for (uint64 k = 0; k < nfw_; ++k) f_[k] = ~0UL;
    if (mp_) f_[nfw_] = mp_;
  }

  // Return whether all bits are set.
  bool all_set_q() const {
    for (uint64 k = 0; k < nfw_; ++k)
      if (~f_[k]) return false;
    if (mp_) {
      uint64 z = f_[nfw_] & mp_;
      if (z != mp_) return false;
    }
    return true;
  }

  // Return whether all bits are clear.
  uint64 all_clear_q() const {
    for (uint64 k = 0; k < nfw_; ++k)
      if (f_[k]) return false;
    if (mp_) {
      uint64 z = f_[nfw_] & mp_;
      if (z != 0) return false;
    }
    return true;
  }

  // Return index of next set bit or value beyond end.
  // Note: the given index n is included in the search
  uint64 next_set(uint64 n) const {
    while ((n < n_) && (!test(n))) ++n;
    return n;
  }

  // Return index of next clear bit or value beyond end.
  // Note: the given index n is included in the search
  uint64 next_clear(uint64 n) const {
    while ((n < n_) && (test(n))) ++n;
    return n;
  }

  // Return index of first set bit or value beyond end.
  uint64 first_set() {
    uint64 k = 0;
    while (k < nfw_) {
      if (f_[k]) break;
      ++k;
    }
    return next_set(k * BITS_PER_WORD);
  }

  // Return index of first clear bit or value beyond end.
  uint64 first_clear() {
    uint64 k = 0;
    while (k < nfw_) {
      if (~f_[k]) break;
      ++k;
    }
    return next_clear(k * BITS_PER_WORD);
  }

  // Return the number of set bits.
  uint64 count_ones() const {
    uint64 ct = 0;
    for (uint64 j = 0; j < nfw_; ++j) ct += bit_count(f_[j]);
    if (mp_) ct += bit_count(mp_ & f_[nfw_]);
    return ct;
  }

  // Return the number of clear bits.
  uint64 count_zeros() const { return n_ - count_ones(); }

  void load(std::ifstream &in) {
    uint64 nfw = nfw_;
    if (mp_) ++nfw;
    assert(f_ && nfw);
    in.read((char *)f_, nfw * sizeof(uint64));
  }

  uint64 serialize(std::ofstream &out) const {
    uint64 nfw = nfw_;
    if (mp_) ++nfw;
    out.write((char *)f_, nfw * sizeof(uint64));
    return nfw * sizeof(uint64);
  }

  void resize(uint64 nbits) {
    assert(nbits > 0);
    uint64 n = nbits;
    uint64 nw = n / BITS_PER_WORD;
    uint64 nbl = n - nw * BITS_PER_WORD;
    uint64 nfw = nw;
    uint64 mp = 0ULL;
    if (0 != nbl) {
      ++nw;
      mp = ~0ULL;
      mp >>= (BITS_PER_WORD - nbl);
    }

    uint64 *f = new uint64[nw]();
    uint64 ncp = mp_ ? nfw_ + 1 : nfw_;
    ncp = ncp < nw ? ncp : nw;
    memcpy(f, f_, ncp * sizeof(uint64));
    delete[] f_;

    f_ = f;
    n_ = n;
    nfw_ = nfw;
    mp_ = mp;
  }

 private:
  uint64 ctor_core(uint64 nbits) {
    assert(nbits > 0);
    n_ = nbits;
    uint64 nw = n_ / BITS_PER_WORD;
    uint64 nbl = n_ - nw * BITS_PER_WORD;
    nfw_ = nw;
    mp_ = 0ULL;
    if (0 != nbl) {
      ++nw;
      mp_ = ~0ULL;
      mp_ >>= (BITS_PER_WORD - nbl);
    }
    return nw;
  }

  uint64 *f_;   // bit bucket
  uint64 n_;    // number of bits
  uint64 nfw_;  // number of words where all bits are used, may be zero
  uint64 mp_;   // mask for partially used word if there is one, else zero
  bool myfq_;   // whether f[] was allocated by class
};

#undef CHECK
#undef DIVMOD
#undef DIVMOD_TEST
