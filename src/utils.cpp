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

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "utils.h"
int remove_dir(const char *dir) {
  char cur_dir[] = ".";
  char up_dir[] = "..";
  char dir_name[512];
  DIR *dirp;
  struct dirent *dp;
  struct stat dir_stat;

  if (0 != access(dir, F_OK)) {
    return 0;
  }

  if (0 > stat(dir, &dir_stat)) {
    perror("get directory stat error");
    return -1;
  }

  if (S_ISREG(dir_stat.st_mode)) {
    remove(dir);
  } else if (S_ISDIR(dir_stat.st_mode)) {
    dirp = opendir(dir);
    while ((dp = readdir(dirp)) != NULL) {
      if ((0 == strcmp(cur_dir, dp->d_name)) ||
          (0 == strcmp(up_dir, dp->d_name))) {
        continue;
      }

      sprintf(dir_name, "%s/%s", dir, dp->d_name);
      remove_dir(dir_name);
    }

    closedir(dirp);
    rmdir(dir);
  } else {
    perror("unknow file type!");
  }

  return 0;
}

std::vector<Range> jobBounds(uint64 num, uint64 thds) {
  std::vector<Range> bounds(thds);
  for (uint64 thd = 0, beg = 0; thd < thds; thd++) {
    bounds[thd].first = beg;
    if (beg < num) beg += (num - beg) / (thds - thd);
    bounds[thd].second = beg;
  }
  return bounds;
}

const uint64 BitMask[65] = {0x0ULL,
                            0x1ULL,
                            0x3ULL,
                            0x7ULL,
                            0xFULL,
                            0x1FULL,
                            0x3FULL,
                            0x7FULL,
                            0xFFULL,
                            0x1FFULL,
                            0x3FFULL,
                            0x7FFULL,
                            0xFFFULL,
                            0x1FFFULL,
                            0x3FFFULL,
                            0x7FFFULL,
                            0xFFFFULL,
                            0x1FFFFULL,
                            0x3FFFFULL,
                            0x7FFFFULL,
                            0xFFFFFULL,
                            0x1FFFFFULL,
                            0x3FFFFFULL,
                            0x7FFFFFULL,
                            0xFFFFFFULL,
                            0x1FFFFFFULL,
                            0x3FFFFFFULL,
                            0x7FFFFFFULL,
                            0xFFFFFFFULL,
                            0x1FFFFFFFULL,
                            0x3FFFFFFFULL,
                            0x7FFFFFFFULL,
                            0xFFFFFFFFULL,
                            0x1FFFFFFFFULL,
                            0x3FFFFFFFFULL,
                            0x7FFFFFFFFULL,
                            0xFFFFFFFFFULL,
                            0x1FFFFFFFFFULL,
                            0x3FFFFFFFFFULL,
                            0x7FFFFFFFFFULL,
                            0xFFFFFFFFFFULL,
                            0x1FFFFFFFFFFULL,
                            0x3FFFFFFFFFFULL,
                            0x7FFFFFFFFFFULL,
                            0xFFFFFFFFFFFULL,
                            0x1FFFFFFFFFFFULL,
                            0x3FFFFFFFFFFFULL,
                            0x7FFFFFFFFFFFULL,
                            0xFFFFFFFFFFFFULL,
                            0x1FFFFFFFFFFFFULL,
                            0x3FFFFFFFFFFFFULL,
                            0x7FFFFFFFFFFFFULL,
                            0xFFFFFFFFFFFFFULL,
                            0x1FFFFFFFFFFFFFULL,
                            0x3FFFFFFFFFFFFFULL,
                            0x7FFFFFFFFFFFFFULL,
                            0xFFFFFFFFFFFFFFULL,
                            0x1FFFFFFFFFFFFFFULL,
                            0x3FFFFFFFFFFFFFFULL,
                            0x7FFFFFFFFFFFFFFULL,
                            0xFFFFFFFFFFFFFFFULL,
                            0x1FFFFFFFFFFFFFFFULL,
                            0x3FFFFFFFFFFFFFFFULL,
                            0x7FFFFFFFFFFFFFFFULL,
                            0xFFFFFFFFFFFFFFFFULL};
