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

#include <type_traits>

#include "types.h"

template <typename TData, typename TUint>
class Array;

using Seq32 = Array<uchar, uint32>;
using Seq64 = Array<uchar, uint64>;

template <typename TData, typename TUint = uint64>
class Array {
  static_assert(std::is_unsigned<TUint>::value, "Array: require unsigned int!");

 public:
  explicit Array() : data_(nullptr), size_(0), cap_(0) {}

  explicit Array(TUint cap) : data_(nullptr), size_(0), cap_(cap) {
    reserve(cap_);
  }

  Array(const Array<TData, TUint>& o) : data_(nullptr), size_(0), cap_(1) {
    *this = o;
  }

  ~Array() { freeData(); }

  inline TData* data() { return data_; }

  inline const TData* data() const { return data_; }

  inline TUint size() const { return size_; }

  inline TUint cap() const { return cap_; }

  void clear() { size_ = 0; }

  void push_back(const TData& val) {
    if (!data_) reserve(1024);
    if (size_ == cap_) recapCopy(cap_ + 1);
    data_[size_++] = val;
  }

  void pop_back() {
    assert(size_);
    --size_;
  }

  inline TData& operator[](TUint idx) {
    assert(idx < size_);
    return data_[idx];
  }

  inline const TData& operator[](TUint idx) const {
    assert(idx < size_);
    return data_[idx];
  }

  inline TData& at(TUint idx) { return operator[](idx); }

  inline const TData& at(TUint idx) const { return operator[](idx); }

  inline TData& front() {
    assert(size_);
    return data_[0];
  }

  inline const TData& front() const { return front(); }

  inline TData& back() {
    assert(size_);
    return data_[size_ - 1];
  }

  inline const TData& back() const {
    assert(size_);
    return data_[size_ - 1];
  }

  inline TData* begin() {
    assert(size_);
    return &data_[0];
  }

  inline const TData* begin() const { return begin(); }

  inline TData* end() { return begin() + size_; }

  inline const TData* end() const { return begin() + size_; }

  Array<TData, TUint>& operator=(const Array<TData, TUint>& o) {
    if (this == &o) return *this;
    if (!o.size()) {
      size_ = 0;
      return *this;
    }
    resize(o.size());
    for (TUint i = 0; i < size_; ++i) data_[i] = o[i];
    return *this;
  }

  Array<TData, TUint>& operator+=(const Array<TData, TUint>& o) {
    for (TUint i = 0; i < o.size(); ++i) push_back(o[i]);
    return *this;
  }

  bool operator==(const Array<TData, TUint>& o) const {
    if (size_ != o.size()) {
      return false;
    }
    for (TUint i = 0; i < size_; ++i) {
      if (data_[i] != o[i]) return false;
    }
    return true;
  }

  void fill(TUint begin, TUint end, const TData& val) {
    assert(size_ && begin < end && end <= size_);
    for (TUint i = begin; i < end; ++i) data_[i] = val;
  }

  void fill(const TData& val) {
    assert(size_);
    for (TUint i = 0; i < size_; ++i) data_[i] = val;
  }

  void zeros(TUint begin, TUint end) {
    assert(size_ && begin < end && end <= size_);
    memset(&data_[begin], 0, sizeof(TData) * (end - begin));
  }

  void zeros() {
    assert(size_);
    memset(data_, 0, sizeof(TData) * size_);
  }

  void reserve(TUint cap) {
    assert(cap);
    if (!data_) lazyInit(cap);
    if (cap <= cap_) {
      if (cap <= size_) size_ = cap;
      cap_ = cap;
      return;
    }
    recapCopy(cap);
  }

  void resize(TUint size) {
    assert(size);
    if (size > cap_) reserve(size);
    size_ = size;
  }

  void erase(TUint idx) {
    assert(idx < size_);
    for (TUint i = idx; i < size_ - 1; ++i) data_[i] = data_[i + 1];
    --size_;
  }

  void erase(TUint idx, TUint len) {
    assert(len && idx + len - 1 < size_);
    for (TUint i = idx; i < size_ - len; ++i) data_[i] = data_[i + len];
    size_ -= len;
  }

  void insert(const TData& val, TUint idx) {
    if (!data_) lazyInit();
    assert(idx <= size_);
    if (size_ == cap_) recapCopy(cap_ + 1);
    for (TUint i = size_; i > idx; --i) data_[i] = data_[i - 1];
    data_[idx] = val;
    ++size_;
  }

 private:
  void lazyInit(TUint cap) {
    assert(!data_);
    data_ = allocMem(cap > cap_ ? cap : cap_);
  }

  TData* allocMem(TUint size) const {
    TData* tmp = new TData[size];
    if (!tmp) {
      std::cerr << "Array: allocate memory failed" << std::endl;
      exit(EXIT_FAILURE);
    }
    return tmp;
  }

  void freeData() {
    if (data_) {
      delete[] data_;
      data_ = nullptr;
    }
  }

  void recapCopy(TUint cap) {
    if (cap <= cap_) return;
    TUint newcap = growCap(cap);
    TData* tmp = allocMem(newcap);
    if (data_) {
      for (TUint i = 0; i < size_; ++i) tmp[i] = data_[i];
      freeData();
    }
    data_ = tmp;
    cap_ = newcap;
  }

  void recapNoCopy(TUint cap) {
    if (cap <= cap_) return;
    TUint newcap = growCap(cap);
    freeData();
    data_ = allocMem(newcap);
    cap_ = newcap;
  }

  TUint growCap(TUint cap) const {
    assert(cap >= cap_);
    if (cap == cap_) return cap;
    TUint newcap = cap_;
    TUint doublecap = newcap * 2;

    if (cap > doublecap) {
      newcap = cap;
    } else {
      if (size_ < 1024) {
        newcap = doublecap;
      } else {
        while (newcap < cap) newcap += newcap / 4;
      }
    }

    if (newcap < cap || newcap >= std::numeric_limits<TUint>::max()) {
      std::cerr << "Array: grow capacity overflow" << std::endl;
      exit(EXIT_FAILURE);
    }

    return newcap;
  }

  TData* data_;
  TUint size_;
  TUint cap_;
};
