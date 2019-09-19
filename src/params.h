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

#include <vector>

#include "tinythread.h"
#include "types.h"

const uint64 MinBufferSize = 16 << 20;
const uint64 MaxBufferSize = 1ULL << 32;
const uint64 DefaultBufferSize = 128 << 20;
const uint64 DefaultThreadNumber = tthread::thread::hardware_concurrency();

struct IndexParam {
  // static const uint64 DefaultK = 16;

  IndexParam() : threads(DefaultThreadNumber) {}

  std::vector<std::string> gfns;
  std::string snpfn;
  std::string output;
  uint64 threads;
};

struct CompParam {
  CompParam() : ntry(3), bufsz(DefaultBufferSize), threads(DefaultThreadNumber) {}

  std::string idxfn;
  std::string outfn;
  uint64 ntry;
  uint64 bufsz;
  uint64 preset;
  uint64 threads;
  std::vector<std::string> infns;
};

struct DecompParam {
  DecompParam() : threads(DefaultThreadNumber), usegz(false) {}

  std::string infn;
  std::string outdir;
  uint64 threads;
  bool usegz;
};

struct MainParam {
  enum class Mode { INDEX, ENCODE, DECODE };

  Mode mode;
  IndexParam index;
  CompParam comp;
  DecompParam decomp;
};
