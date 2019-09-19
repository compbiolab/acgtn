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

#include <tsl/robin_map.h>
#include <string>
#include <vector>

#include "bitarray.h"
#include "intvector.h"
#include "params.h"

struct Dnaseq {
  uchar* dna;   // input
  uint64 len;   // input
  uint64 ntry;  // input

  uint64 idx;  // output
  uint64 pos;  // output
  bool isfwd;  // output

  Dnaseq() : dna(nullptr), len(0), ntry(2), idx(0), pos(0), isfwd(true) {}
};

class Index {
 public:
  static constexpr uint64 K = 14;
  static constexpr uint64 N = 5;
  static constexpr uint64 L = 54;  // K + 8 * N;

  uint64 PosBits;
  uint64 SeqBits;
  uint64 PosMask;
  uint64 RefSize;

  Index() {}
  Index(const IndexParam& param);
  ~Index();

  int build();
  bool load(const std::string file);

  bool encode(Dnaseq& read);

  inline uchar* pos2seq(const uint64 pos, uchar seq[L]) const;

  uchar* refseq() const { return seq_; }
  uint64 refsize() const { return gsz_; }

  void printInfo();

 private:
  void loadGenome(tsl::robin_map<std::string, Range>& chrmap);
  void loadSNP(const tsl::robin_map<std::string, Range>& chrmap);
  void computeSA();
  void computeKmerMap();

  bool extendSeq(Dnaseq& read);
  uint64 estimateNSA(uint64* bucketCnts, const uint64 k) const;

  static void sortSAWorker(void* param);
  static void kmerMapWorker(void* param);
  static void estimateNSAWorker(void* param);

  IndexParam param_;
  uchar* seq_;
  bitarray* kfwd_;
  // bitarray* kbwd_;
  intvector<> sa_;
  intvector<> kmap_;
  tsl::robin_map<uint64, uchar> snpmap_;
  uint64 gsz_;
  uint64 nsa_;
};
