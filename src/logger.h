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

#include <ctime>
#include <iostream>

#define STR_(X) #X
#define STR(X) STR_(X)

#define LOGGER_LEVEL(level)                                    \
  std::time_t t = std::time(nullptr);                          \
  char mbstr[100];                                             \
  if (std::strftime(mbstr, sizeof(mbstr),                      \
                    "[" STR(level) " | %Y-%m-%d %H:%M:%S] | ", \
                    std::localtime(&t))) {                     \
    out << mbstr << msg << std::endl;                          \
  }

struct Logger {
  static bool verbose;

  static void Info(const std::string &msg, std::ostream &out = std::cerr) {
    if (!verbose) return;
    LOGGER_LEVEL(INFO);
  }

  static void Warning(const std::string &msg, std::ostream &out = std::cerr) {
    LOGGER_LEVEL(WARNING);
  }

  static void Error(const std::string &msg, std::ostream &out = std::cerr) {
    LOGGER_LEVEL(ERROR);
    exit(EXIT_FAILURE);
  }

#ifndef NDEBUG
  static void Debug(const std::string &msg, std::ostream &out = std::cerr) {
    LOGGER_LEVEL(DEBUG);
  }
#endif
};

#undef LOGGER_LEVEL
#undef STR
#undef STR_
