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

#include <rocksdb/db.h>

#include "fqreader.h"
#include "header.h"
#include "index.h"
#include "params.h"
#include "roaringcpp.h"

class Compressor {
 public:
  Compressor();
  Compressor(const CompParam& param);
  ~Compressor();

  int Run();

 private:
  struct FastqInfo;

  void ReorderReads();
  void EncodeReads();

  static void AssignKeysWorker(void* param);
  static void ReadInputsWorker(void* param);
  static void ReloadReadsWorker(void* param);
  static void EncodeReadsWorker(void* param);

  static void PutFastq(rocksdb::WriteBatch* batch, const FastqInfo& info);
  static void GetFastq(rocksdb::Iterator* it, FastqInfo& info);

  rocksdb::Options GetOptions();
  void OpenDb(const std::string& dir, bool readOnly = false);

  CompParam param_;
  Header header_;
  Roaring isrc_;

  Vstring fns_;
  Vuint64 codes_;
  
#ifndef NDEBUG
  uint64 nsuc_;
  uint64 nrds_;
  uint64 nfld_;
#endif

  rocksdb::DB* db_;
  std::string dbdir_;
};
