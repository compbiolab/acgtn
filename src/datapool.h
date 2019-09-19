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

#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

template <class DataType>
class DataPool {
  using Pool = std::vector<DataType*>;

 public:
  static const uint64 DefaultPoolCapacity = 32;

  DataPool(uint64 capacity = DefaultPoolCapacity) : cap_(capacity), size_(0) {
    items_.resize(cap_);
    addrs_.reserve(cap_);
  }

  ~DataPool() {
    for (auto p = addrs_.begin(); p != addrs_.end(); ++p) delete *p;
  }

  template <class... ArgsType>
  DataType* Acquire(ArgsType&&... args) {
    std::unique_lock<std::mutex> lock(mutex_);

    while (size_ >= cap_) full_.wait(lock);

    assert(items_.size());

    DataType*& tmp = items_.back();
    items_.pop_back();
    if (!tmp) {
      tmp = new DataType(std::forward<ArgsType>(args)...);
      addrs_.push_back(tmp);
    }

    size_++;
    return tmp;
  }

  void Release(const DataType* data) {
    std::lock_guard<std::mutex> lock(mutex_);

    assert(data && size_ && size_ <= cap_);
    assert(std::find(addrs_.begin(), addrs_.end(), data) != addrs_.end());

    items_.push_back((DataType*)data);
    size_--;

    full_.notify_one();
  }

 private:
  const uint64 cap_;
  uint64 size_;

  Pool items_;
  Pool addrs_;

  std::mutex mutex_;
  std::condition_variable full_;
};
