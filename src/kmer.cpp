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

#include <cassert>
#include <iostream>

#include "kmer.h"
#include "utils.h"

const uchar Kmer::id2base[5] = {'A', 'C', 'G', 'T', 'N'};

const uint64 Kmer::base2id[256] = {
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 3, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

const uchar Kmer::base2rev[256] = {
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', 'T',  '\0', 'G',  '\0', '\0', '\0', 'C',
    '\0', '\0', '\0', '\0', '\0', '\0', 'N',  '\0', '\0', '\0', '\0', '\0',
    'A',  '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', 'T',  '\0', 'G',  '\0', '\0', '\0', 'C',  '\0', '\0', '\0', '\0',
    '\0', '\0', 'N',  '\0', '\0', '\0', '\0', '\0', 'A',  '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0'};

bool Kmer::isbase(const uchar base) { return 4ULL > base2id[base]; }
bool Kmer::isnbase(const uchar base) { return 5ULL != base2id[base]; }

bool Kmer::hasN(const uchar* kmer, uint64 k) {
  for (uint64 i = 0; i < k; ++i) {
    if (4ULL == base2id[kmer[i]]) return true;
  }
  return false;
}

bool Kmer::hasN(const std::string& kmer) {
  for (const uchar b : kmer) {
    if (4ULL == base2id[b]) return true;
  }
  return false;
}

uint64 Kmer::toId(const uchar* kmer, uint64 k) {
  uint64 id = 0;
  for (uint64 i = 0; i < k; ++i) {
    id <<= 2;
    id |= base2id[kmer[i]];
  }
  return id;
}

uint64 Kmer::nextId(const uint64 id, const uchar b, uint64 k) {
  assert(k <= 32);
  return ((id << 2) | base2id[b]) & BitMask[2 * k];
}

uint64 Kmer::prevId(const uint64 id, const uchar b, uint64 k) {
  return ((id >> 2) & BitMask[2 * k - 2]) | (base2id[b] << (2 * k - 2));
}

uint64 Kmer::toId(const std::string& kmer) {
  uint64 id = 0;
  for (uint64 i = 0; i < kmer.size(); ++i) {
    id <<= 2;
    id |= base2id[(uchar)kmer[i]];
  }
  return id;
}

void Kmer::toSeq(const uint64 id, uint64 k, uchar* kmer) {
  for (uint64 i = 0; i < k; ++i) kmer[k - i - 1] = id2base[(id >> 2 * i) & 0x3];
}

std::string Kmer::toStr(uint64 id, uint64 k) {
  std::string kmer(k, '\0');
  for (uint64 i = 0; i < k; ++i) kmer[k - i - 1] = id2base[(id >> 2 * i) & 0x3];
  return kmer;
}

uint64 Kmer::rcId(const uchar* kmer, uint64 k) {
  uint64 id = 0;
  for (uint64 i = 0; i < k; ++i) {
    id <<= 2;
    id |= 3ULL - base2id[kmer[k - i - 1]];
    // id |= base2rc[*(kmer - i)];
  }
  return id;
}

uint64 Kmer::rcId(const uint64 id, uint64 k) {
  uint64 rcid = 0;
  for (uint64 i = 0; i < k; ++i) {
    rcid <<= 2;
    // rcid |= base2id[int2base((id >> 2 * i) & 0x3)];
    rcid |= 3ULL - ((id >> 2 * i) & 0x3);
  }
  return rcid;
}

void Kmer::rcSeq(Seq32& kmer) {
  const uint32 k = kmer.size();
  for (uint32 i = 0; i<k>> 1; ++i) {
    uchar tmp = kmer[k - i - 1];
    kmer[k - i - 1] = base2rev[kmer[i]];
    kmer[i] = base2rev[tmp];
  }
}

void Kmer::rcSeq(uchar* kmer, uint64 k) {
  for (uint64 i = 0; i<k>> 1; ++i) {
    uchar tmp = kmer[k - i - 1];
    kmer[k - i - 1] = base2rev[kmer[i]];
    kmer[i] = base2rev[tmp];
  }
}

void Kmer::rcSeq(const uchar* kmer, uint64 k, uchar* rcseq) {
  for (uint64 i = 0; i < k; ++i) {
    rcseq[k - i - 1] = base2rev[kmer[i]];
  }
}

std::string Kmer::rcStr(const uchar* kmer, uint64 k) {
  std::string rcstr(k, '\0');
  for (uint64 i = 0; i < k; ++i) {
    rcstr[k - i - 1] = base2rev[kmer[i]];
  }
  return rcstr;
}

std::string Kmer::rcStr(const std::string& kmer) {
  const uint64 k = kmer.size();
  std::string rcstr(k, '\0');
  for (uint64 i = 0; i < k; ++i) {
    rcstr[k - i - 1] = base2rev[(uchar)kmer[i]];
  }
  return rcstr;
}

uint64 Kmer::revId(const uchar* kmer, uint64 k) {
  uint64 id = 0;
  for (uint64 i = 0; i < k; ++i) {
    id <<= 2;
    id |= base2id[kmer[k - i - 1]];
  }
  return id;
}

uint64 Kmer::minId(const uchar* kmer, uint64 k) {
  bool isRc = false;
  for (uint64 i = 0; i < k; ++i) {
    if (kmer[i] > base2rev[kmer[k - i - 1]]) {
      isRc = true;
      break;
    }
  }
  return isRc ? rcId(kmer, k) : toId(kmer, k);
}
