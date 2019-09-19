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
#include "params.h"
#include "roaringcpp.h"
#include "strmcoder.h"

class Decompressor {
 public:
  Decompressor();
  Decompressor(const DecompParam& param);
  ~Decompressor();

  int Run();

 private:
  void ReorderReads();
  void AssignReads();

  static void ReadInputWorker(void* param);
  static void DecodeReadsWorker(void* param);
  static void ReloadReadsWorker(void* param);

  static void PutFastq(rocksdb::WriteBatch* batch, const Fastq& fastq);
  static void GetFastq(rocksdb::Iterator* it, Fastq& fastq);

  rocksdb::Options GetOptions();
  void OpenDb(const std::string& dir, bool readOnly = false);

  DecompParam param_;
  Header header_;
  Roaring isrc_;
  uint64 bufsz_;

  rocksdb::DB* db_;
  std::string dbdir_;
};
