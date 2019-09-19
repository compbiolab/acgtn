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

#include <zlib.h>
#include <vector>

#include "crc64.h"
#include "header.h"
#include "kseq.h"
#include "myarray.h"
#include "types.h"
#include "utils.h"

KSEQ_INIT(gzFile, gzread);

struct Fastq {
  Fastq() : key(0), seq() {}

  Fastq& operator=(const Fastq& other) {
    if (this != &other) {
      key = other.key;
      seq = other.seq;
    }
    return *this;
  }

  uint64 key;
  Seq32 seq;
};

using Fastqs = Array<Fastq>;

class Fqreader {
  static const uint64 BufferSize = 256 << 20;
  static const uint64 SwapSize = 1024;

 public:
  Fqreader(const std::string fn, uint64 bufsz = BufferSize);
  Fqreader(const std::vector<std::string>& fns, uint64 bufsz = BufferSize);
  ~Fqreader();

  bool ReadFastqs(Fastqs& fqs);

  Header GetHeader() const { return hdr_; }

 private:
  bool OpenNext();
  void CloseNext();

  Header hdr_;
  Vstring fns_;
  uint64 bufsz_;
  uint64 idx_;
  uint64 rid_;
  gzFile fp_;
  kseq_t* ks_;
  std::vector<crc64> crcs_;
};
