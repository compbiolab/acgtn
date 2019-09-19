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
#include <queue>
#include <thread>

template <class DataType>
class DataQueue {
  using Item = std::pair<uint64, DataType*>;
  using Queue = std::deque<Item>;

 public:
  static const uint64 DefaultQueueCapacity = 32;

  DataQueue(uint64 capacity = DefaultQueueCapacity, uint64 njobs = 1)
      : cap_(capacity), njobs_(njobs), ndone_(0) {}

  ~DataQueue() {}

  void Done() {
    std::lock_guard<std::mutex> lock(mutex_);
    assert(ndone_ < njobs_);
    ndone_++;
    empty_.notify_all();
  }

  void Push(uint64 id, const DataType* data) {
    std::unique_lock<std::mutex> lock(mutex_);

    while (items_.size() > cap_ && id > items_.front().first) full_.wait(lock);

    items_.push_back(std::make_pair(id, (DataType*)data));
    if (items_.size() > 1)
      std::sort(items_.begin(), items_.end(),
                [](const Item& a, const Item& b) { return a.first < b.first; });

    empty_.notify_one();
  }

  bool Pop(uint64& id, DataType*& data) {
    std::unique_lock<std::mutex> lock(mutex_);

    while (!items_.size() && ndone_ != njobs_) empty_.wait(lock);

    if (items_.size()) {
      id = items_.front().first;
      data = items_.front().second;
      items_.pop_front();
      full_.notify_one();
      return true;
    }

    assert(ndone_ == njobs_ && !items_.size());
    return false;
  }

 private:
  const uint64 cap_;
  const uint64 njobs_;
  uint64 ndone_;
  Queue items_;

  std::mutex mutex_;
  std::condition_variable full_;
  std::condition_variable empty_;
};
