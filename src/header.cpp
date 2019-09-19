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

#include "header.h"
#include "logger.h"
#include "strmcoder.h"

Header::Header()
    : tag(TAG), version(VERSION), bufsz(0), bases(0), minlen(-1), maxlen(0) {}

void Header::serialize(FILE* out) const {
  if (reads.size() != inputs.size() || maxlen < minlen || bufsz < MinBufsz)
    Logger::Error("Bad acgtn header");

  std::string fns = "";
  for (uint64 i = 0; i < inputs.size(); ++i) {
    assert(inputs[i] != "");
    fns += inputs[i] + ',';
  }

  fwrite(&tag, sizeof(uint64), 1, out);
  fwrite(&version, sizeof(uint64), 1, out);
  fwrite(&bufsz, sizeof(uint64), 1, out);
  fwrite(&bases, sizeof(uint64), 1, out);
  fwrite(&minlen, sizeof(uint64), 1, out);
  fwrite(&maxlen, sizeof(uint64), 1, out);

  std::string tmp = CompressWithLzma(fns, 6);
  uint64 tmpsz = tmp.size();
  fwrite(&tmpsz, sizeof(uint64), 1, out);
  fwrite(tmp.data(), 1, tmpsz, out);

  for (uint64 num : reads) fwrite(&num, sizeof(uint64), 1, out);
  for (uint64 crc : crc64s) fwrite(&crc, sizeof(uint64), 1, out);
}

void Header::load(FILE* in) {
  assert(1 == fread(&tag, sizeof(uint64), 1, in));
  assert(1 == fread(&version, sizeof(uint64), 1, in));
  assert(1 == fread(&bufsz, sizeof(uint64), 1, in));
  assert(1 == fread(&bases, sizeof(uint64), 1, in));
  assert(1 == fread(&minlen, sizeof(uint64), 1, in));
  assert(1 == fread(&maxlen, sizeof(uint64), 1, in));

  if (TAG != tag || maxlen < minlen || bufsz < MinBufsz)
    Logger::Error("Bad acgtn header");

  reads.clear();
  crc64s.clear();
  inputs.clear();

  uint64 tmpsz = 0;
  fread(&tmpsz, sizeof(uint64), 1, in);
  std::string tmp;
  tmp.resize(tmpsz);
  char* ctmp = (char*)tmp.data();
  fread(ctmp, 1, tmpsz, in);
  std::string fns = DecompressWithLzma(tmp);

  std::string fn = "";
  for (uint64 i = 0; i < fns.size(); ++i) {
    if (fns[i] == ',') {
      inputs.push_back(fn);
      fn = "";
    } else {
      fn += fns[i];
    }
  }

  for (uint64 i = 0; i < inputs.size(); ++i) {
    uint64 num;
    assert(1 == fread(&num, sizeof(uint64), 1, in));
    reads.push_back(num);
  }
  
  for (uint64 i = 0; i < inputs.size(); ++i) {
    uint64 crc;
    assert(1 == fread(&crc, sizeof(uint64), 1, in));
    crc64s.push_back(crc);
  }
}

void Header::print(bool verbose) const {
  std::cout << "raw file name\t->\tdecode to file" << std::endl;
  for (uint64 i = 0; i < inputs.size(); ++i) {
    std::string input = inputs[i];
    auto pos = input.find_last_of('.');
    assert(pos != input.npos);
    std::string outfn = input.substr(0, pos) + ".dna(.gz)";
    std::cout << input << "\t->\t" << outfn << std::endl;
  }

  if (!verbose) return;
  std::cout << std::endl;

  uint64 num = reads.back();
  std::cout << "# of total reads: " << std::to_string(num) << std::endl;
  std::cout << "# of total bases: " << std::to_string(bases) << std::endl;
  std::cout << "min length of reads : " << std::to_string(minlen) << std::endl;
  std::cout << "max length of reads : " << std::to_string(maxlen) << std::endl;
  std::cout << std::endl;

  std::cout << "file name: # of reads" << std::endl;
  uint64 tmp = 0;
  for (uint64 i = 0; i < inputs.size(); ++i) {
    assert(reads[i] > tmp);
    std::cout << inputs[i] << ": " << reads[i] - tmp << std::endl;
    tmp = reads[i];
  }
  std::cout << std::endl;
}
