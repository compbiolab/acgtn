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

#include "utils.h"

struct Header {
  const static uint64 TAG = 0x7F414347544EFE;
  const static uint64 VERSION = 1ULL << 48;  // major:16, minor:16, patch:32
  const static uint64 MinBufsz = 64 << 20;

  uint64 tag;
  uint64 version;
  uint64 bufsz;
  uint64 bases;
  uint64 minlen;
  uint64 maxlen;
  Vuint64 reads;
  Vuint64 crc64s;
  Vstring inputs;

  Header();

  void serialize(FILE* out) const;
  void load(FILE* in);
  void print(bool verbose = false) const;
};
