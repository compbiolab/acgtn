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

#include "utils.h"

#include <cassert>
#include <initializer_list>
#include <iostream>

template <typename t_intvector>
class intvector_reference;

template <uint64 t_width = 0>
class intvector {
 public:
  using intvector_type = intvector<t_width>;
  using reference = intvector_reference<intvector_type>;

  intvector(uint64 size = 0, uint64 width = t_width)
      : data_(nullptr),
        size_(size),
        width_(width),
        nbits_(size * width),
        nwords_(0) {
    bit_resize(nbits_);
  }

  intvector(const intvector& v);
  intvector& operator=(const intvector& v);

  ~intvector() { FreeArrPtr(data_); }

  inline reference operator[](const uint64 idx);
  inline uint64 operator[](const uint64 idx) const;

  inline reference at(const uint64 idx);
  inline uint64 at(const uint64 idx) const;

  uint64 size() const { return size_; }
  uint64 width() const { return width_; }
  uint64* data() const { return data_; }

  void set_size(const uint64 size);
  void set_width(const uint64 width);
  void set_size_width(const uint64 size, const uint64 width);

  void load(std::ifstream& in);
  uint64 serialize(std::ofstream& out) const;

 private:
  static_assert(t_width <= 64, "intvector: width of must be at most 64bits.");

  void bit_resize(const uint64 nbits);

  uint64* data_;
  uint64 size_;
  uint64 width_;
  uint64 nbits_;
  uint64 nwords_;
};

template <uint64 t_width>
intvector<t_width>::intvector(const intvector<t_width>& v)
    : data_(nullptr), size_(v.size()), width_(v.width()), nwords_(0) {
  size_ = v.size();
  width_ = v.width();
  bit_resize(size_ * width_);
  memcpy(data_, v.data(), nwords_ * sizeof(uint64));
}

template <uint64 t_width>
intvector<t_width>& intvector<t_width>::operator=(const intvector<t_width>& v) {
  if (this != &v) {
    size_ = v.size();
    width_ = v.width();
    bit_resize(size_ * width_);
    memcpy(data_, v.data(), nwords_ * sizeof(uint64));
  }
  return *this;
}

static inline void bits_write_int(uint64* word, uint64 val, uint64 offset,
                                  const uint64 len);

static inline uint64 bits_read_int(const uint64* word, uint64 offset,
                                   const uint64 len);

template <uint64 t_width>
inline auto intvector<t_width>::operator[](const uint64 idx) -> reference {
  assert(idx < this->size() && width_);
  uint64 i = idx * width_;
  return reference(data_ + (i >> 6), i & 0x3F, width_);
}

template <uint64 t_width>
inline uint64 intvector<t_width>::operator[](const uint64 idx) const {
  assert(idx < this->size() && width_);
  uint64 i = idx * width_;
  return bits_read_int(data_ + (i >> 6), idx & 0x3F, width_);
}

template <uint64 t_width>
inline auto intvector<t_width>::at(const uint64 idx) -> reference {
  assert(idx < this->size() && width_);
  uint64 i = idx * width_;
  return reference(data_ + (i >> 6), i & 0x3F, width_);
}

template <uint64 t_width>
inline uint64 intvector<t_width>::at(const uint64 idx) const {
  assert(idx < this->size() && width_);
  uint64 i = idx * width_;
  return bits_read_int(data_ + (i >> 6), idx & 0x3F, width_);
}

template <uint64 t_width>
void intvector<t_width>::set_size(const uint64 size) {
  if (size == size_) return;
  bit_resize(size * width_);
  size_ = size;
}

template <uint64 t_width>
void intvector<t_width>::set_width(const uint64 width) {
  assert(width <= 64);
  if (width == width_) return;
  bit_resize(size_ * width);
  width_ = width;
}

template <uint64 t_width>
void intvector<t_width>::set_size_width(const uint64 size, const uint64 width) {
  set_size(size);
  set_width(width);
}

template <uint64 t_width>
void intvector<t_width>::bit_resize(const uint64 nbits) {
  if ((0 == nbits) || (nbits_ == nbits && nwords_ != 0)) return;
  uint64 nw = nbits / 64;
  if (0 != (nbits - nw * 64)) nw++;
  if (data_) delete[] data_;
  data_ = new uint64[nw]();
  nbits_ = nbits;
  nwords_ = nw;
}

template <uint64 t_width>
void intvector<t_width>::load(std::ifstream& in) {
  uint64 size, width;
  read_member(size, in);
  read_member(width, in);
  set_size_width(size, width);
  in.read((char*)data_, nwords_ * sizeof(uint64));
}

template <uint64 t_width>
uint64 intvector<t_width>::serialize(std::ofstream& out) const {
  uint64 written_bytes = 0;
  written_bytes += write_member(size_, out);
  written_bytes += write_member(width_, out);
  uint64 nbytes = nwords_ * sizeof(uint64);
  out.write((char*)data_, nbytes);
  return written_bytes + nbytes;
}

template <typename t_intvector>
class intvector_reference {
 public:
  intvector_reference(uint64* word, uint64 offset, uint64 len)
      : word_(word), offset_(offset), len_(len){};

  intvector_reference& operator=(uint64 val) {
    bits_write_int(word_, val, offset_, len_);
    return *this;
  };

  intvector_reference& operator=(const intvector_reference& x) {
    return *this = uint64(x);
  };

  operator uint64() const { return bits_read_int(word_, offset_, len_); }

  intvector_reference& operator++() {
    uint64 x = bits_read_int(word_, offset_, len_);
    bits_write_int(word_, x + 1, offset_, len_);
    return *this;
  }

  uint64 operator++(int) {
    uint64 val = (uint64) * this;
    ++(*this);
    return val;
  }

  intvector_reference& operator--() {
    uint64 x = bits_read_int(word_, offset_, len_);
    bits_write_int(word_, x - 1, offset_, len_);
    return *this;
  }

  uint64 operator--(int) {
    uint64 val = (uint64) * this;
    --(*this);
    return val;
  }

  intvector_reference& operator+=(const uint64 val) {
    uint64 x = bits_read_int(word_, offset_, len_);
    bits_write_int(word_, x + val, offset_, len_);
    return *this;
  }

  intvector_reference& operator-=(const uint64 val) {
    uint64 x = bits_read_int(word_, offset_, len_);
    bits_write_int(word_, x - val, offset_, len_);
    return *this;
  }

 private:
  uint64* word_;
  const uint64 offset_;
  const uint64 len_;
};

uint64 bits_read_int(const uint64* word, uint64 offset, const uint64 len) {
  const uint64 w = (*word) >> offset;
  if (offset + len > 64) {
    return w |
           ((*(word + 1) & BitMask[(offset + len) & 0x3F]) << (64 - offset));
  } else {
    return w & BitMask[len];
  }
}

void bits_write_int(uint64* word, uint64 val, uint64 offset, const uint64 len) {
  val &= BitMask[len];
  if (offset + len < 64) {
    *word &= (~(BitMask[len] << offset));
    *word |= (val << offset);
  } else {
    *word &= BitMask[offset];
    *word |= (val << offset);
    if ((offset = (offset + len) & 0x3F)) {
      *(word + 1) &= (~BitMask[offset]);
      *(word + 1) |= (val >> (len - offset));
    }
  }
}
