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

#include <string>

#include "myarray.h"

struct Kmer {
  static const uchar id2base[5];
  static const uint64 base2id[256];
  static const uchar base2rev[256];

  static bool isbase(const uchar base);
  static bool isnbase(const uchar base);
  static bool hasN(const uchar* kmer, uint64 k);
  static bool hasN(const std::string& kmer);

  static uint64 toId(const uchar* kmer, uint64 k);
  static uint64 nextId(const uint64 id, const uchar b, uint64 k);
  static uint64 prevId(const uint64 id, const uchar b, uint64 k);

  static uint64 toId(const std::string& kmer);
  static void toSeq(const uint64 id, uint64 k, uchar* kmer);
  static std::string toStr(uint64 id, uint64 k);

  static uint64 rcId(const uchar* kmer, uint64 k);
  static uint64 rcId(const uint64 id, uint64 k);

  static void rcSeq(Seq32& kmer);
  static void rcSeq(uchar* kmer, uint64 k);
  static void rcSeq(const uchar* kmer, uint64 k, uchar* rcseq);

  static std::string rcStr(const uchar* kmer, uint64 k);
  static std::string rcStr(const std::string& kmer);

  static uint64 revId(const uchar* kmer, uint64 k);
  static uint64 minId(const uchar* kmer, uint64 k);
};
