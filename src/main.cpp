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

#include "compressor.h"
#include "decompressor.h"
#include "header.h"
#include "logger.h"
#include "main.h"

#include <getopt.h>
#include <set>

bool Logger::verbose = false;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    usage();
    return -1;
  }

  MainParam param;
  if (!parse_arguments(argc, argv, param)) return -1;

  if (param.mode == MainParam::Mode::ENCODE) return compress(param.comp);
  if (param.mode == MainParam::Mode::DECODE) return decompress(param.decomp);
  return index(param.index);
}

void usage() {
  char cwd[PATH_MAX + 1];
  if (!getcwd(cwd, sizeof(cwd)))
    Logger::Error("Cannot get current work directory");

  std::string version = std::to_string(Header::VERSION >> 48) + '.';
  version += std::to_string((Header::VERSION >> 32) & BitMask[16]) + '.';
  version += std::to_string(Header::VERSION & BitMask[32]);

  std::cerr << "acgtn - Archival Compression for hiGh-Throughput sequeNcing "
               "data (v"
            << version << ")\n\n";
  std::cerr << "usage:\n\tacgtn <index|encode|decode> [options] <files|file "
               "list>\n\n";

  std::cerr << "<index> options:\n";
  std::cerr << "\t-o <str>    : index prefix, optional\n";
  std::cerr << "\t-s <str>    : common snp file, optional\n";
  std::cerr << "\t-t <int>    : worker threads number, default: "
            << DefaultThreadNumber << '\n';
  std::cerr << "\t-l          : genome files are provided in list\n";
  std::cerr << "\t-v          : verbose mode, default: false\n";
  std::cerr << std::endl;

  std::cerr << "<encode> options:\n";
  std::cerr << "\t-i <str>    : index file, required\n";
  std::cerr << "\t-o <str>    : output file prefix, optional\n";
  std::cerr << "\t-b <str>    : block size to work, default: 256M\n";
  std::cerr << "\t-n <int>    : trying times to assign reads, default: 3\n";
  std::cerr
      << "\t-p <int>    : compression preset level, must <= 9, default: 6\n";
  std::cerr << "\t-t <int>    : worker threads number, default: "
            << DefaultThreadNumber << '\n';
  std::cerr << "\t-l          : input files are provided in list\n";
  std::cerr << "\t-v          : verbose mode, default: false\n";
  std::cerr << std::endl;

  std::cerr << "<decode> options:\n";
  std::cerr << "\t-o <str>    : output directory, default: " << cwd << "\n";
  std::cerr << "\t-t <int>    : worker threads number, default: "
            << DefaultThreadNumber << '\n';
  std::cerr
      << "\t-g          : compress output files using gzip, default: false\n";
  std::cerr << "\t-v          : verbose mode, default: false\n";
  std::cerr << std::endl;
}

int index(const IndexParam& param) {
  Index index(param);
  return index.build();
}

int compress(const CompParam& param) {
  Compressor compressor(param);
  compressor.Run();
  return 0;
}

int decompress(const DecompParam& param) {
  Decompressor decompressor(param);
  decompressor.Run();
  return 0;
}

bool parse_arguments(int argc, char* argv[], MainParam& param) {
  std::string mode = argv[1];
  if ((mode != "index") && (mode != "encode") && (mode != "decode")) {
    std::cerr << "invalid mode: " << mode
              << ", mode must be <index|encode|decode>\n"
              << std::endl;
    usage();
    return false;
  }

  if (mode == "index") {
    param.mode = MainParam::Mode::INDEX;
    bool inlist = false;

    const char* optstr = "o:s:t:lvh";
    static struct option options[] = {{"output", required_argument, NULL, 'o'},
                                      {"snpfile", optional_argument, NULL, 's'},
                                      {"threads", optional_argument, NULL, 't'},
                                      {"verbose", no_argument, NULL, 'v'},
                                      {"help", no_argument, NULL, 'h'}};
    int c;
    int optidx;
    optind++;
    while (true) {
      c = getopt_long(argc, argv, optstr, options, &optidx);
      if (c == -1) break;
      switch (c) {
        case 0:
          break;
        case 'o':
          param.index.output = optarg;
          break;
        case 's':
          param.index.snpfn = optarg;
          break;
        case 't':
          param.index.threads = std::stoull(optarg);
          break;
        case 'l':
          inlist = true;
          break;
        case 'v': {
          Logger::verbose = true;
          break;
        }
        case 'h':
          usage();
          break;
        default:
          break;
      }
    }

    char actualpath[PATH_MAX + 1];
    if (inlist) {
      if (optind + 1 != argc)
        Logger::Error("List of inputs must be in one file!");
      std::string flist = argv[optind];
      if (!realpath(flist.c_str(), actualpath))
        Logger::Error("Invalid input file list " + flist);
      flist = actualpath;
      std::ifstream istrm(flist, std::ios::binary);
      std::string infn;
      while (istrm >> infn) {
        if (!realpath(infn.c_str(), actualpath))
          Logger::Error("Invalid input file " + infn);
        infn = actualpath;
        param.index.gfns.push_back(infn);
      }
    } else {
      for (int i = optind; i < argc; ++i) {
        if (!realpath(argv[i], actualpath))
          Logger::Error("Invalid input file " + std::string(argv[i]));
        std::string infn = actualpath;
        param.index.gfns.push_back(infn);
      }
    }

    const std::vector<std::string>& gfns = param.index.gfns;

    if (!gfns.size()) Logger::Error("Empty genome reference file");

    std::set<std::string> gfnset(gfns.begin(), gfns.end());
    if (gfns.size() != gfnset.size())
      Logger::Error("Duplicate genome reference files specified");

    if (!param.index.output.size()) Logger::Error("No output prefix specified");

    return true;
  }

  if (mode == "encode") {
    param.mode = MainParam::Mode::ENCODE;
    param.comp.preset = 6;
    bool inlist = false;

    const char* optstr = "i:o:n:b:p:t:lvh";
    static struct option options[] = {{"index", required_argument, NULL, 'i'},
                                      {"output", required_argument, NULL, 'o'},
                                      {"ntry", optional_argument, NULL, 'n'},
                                      {"bufsz", optional_argument, NULL, 'b'},
                                      {"preset", optional_argument, NULL, 'p'},
                                      {"threads", optional_argument, NULL, 't'},
                                      {"list", no_argument, NULL, 'l'},
                                      {"verbose", no_argument, NULL, 'v'},
                                      {"help", no_argument, NULL, 'h'}};
    int c;
    int optidx;
    optind++;
    while (true) {
      c = getopt_long(argc, argv, optstr, options, &optidx);
      if (c == -1) break;
      switch (c) {
        case 0:
          break;
        case 'i':
          param.comp.idxfn = optarg;
          break;
        case 'o':
          param.comp.outfn = optarg;
          break;
        case 'n':
          param.comp.ntry = std::stoull(optarg);
          break;
        case 'b': {
          char* p;
          uint64 sz = strtoul(optarg, &p, 10);
          if (*p == 'M' || *p == 'm')
            sz *= MEGABYTE;
          else if (*p == 'G' || *p == 'g')
            sz *= GIGABYTE;
          else {
            usage();
            exit(EXIT_FAILURE);
          }
          param.comp.bufsz = sz;
          break;
        }
        case 'p':
          param.comp.preset = std::stoull(optarg);
          break;
        case 't':
          param.comp.threads = std::stoull(optarg);
          break;
        case 'l':
          inlist = true;
          break;
        case 'v':
          Logger::verbose = true;
          break;
        case 'h':
          usage();
          break;
        default:
          break;
      }
    }

    char actualpath[PATH_MAX + 1];
    if (inlist) {
      if (optind + 1 != argc)
        Logger::Error("List of inputs must be in one file!");
      std::string flist = argv[optind];
      if (!realpath(flist.c_str(), actualpath))
        Logger::Error("Invalid input file list " + flist);
      flist = actualpath;
      std::ifstream istrm(flist, std::ios::binary);
      std::string infn;
      while (istrm >> infn) {
        if (!realpath(infn.c_str(), actualpath))
          Logger::Error("Invalid input file " + infn);
        infn = actualpath;
        param.comp.infns.push_back(infn);
      }
    } else {
      for (int i = optind; i < argc; ++i) {
        if (!realpath(argv[i], actualpath))
          Logger::Error("Invalid input file " + std::string(argv[i]));
        std::string infn = actualpath;
        param.comp.infns.push_back(infn);
      }
    }

    if (!param.comp.infns.size()) Logger::Error("Empty input file");

    std::vector<std::string> infns;
    for (const auto& fn : param.comp.infns) {
      const auto pos = fn.find_last_of('/');
      assert(pos != fn.npos);
      infns.push_back(fn.substr(pos + 1));
    }

    std::set<std::string> infnset(infns.begin(), infns.end());
    if (infns.size() != infnset.size())
      Logger::Error("Basename of inputs cannot be identical!");

    if (param.comp.outfn == "") {
      srand(time(NULL));
      std::string outfn = "tmp" + std::to_string(rand());
      param.comp.outfn = outfn;
      Logger::Warning("Empty output file, output to " + outfn + ".acgtn");
    }

    if (!param.comp.idxfn.size()) Logger::Error("No index file specified");

    if (!param.comp.ntry) Logger::Error("Try assign read cannot be 0 time");

    if (param.comp.bufsz < MinBufferSize) {
      Logger::Warning("Too small bufsz setting, reset to 256M");
      param.comp.bufsz = DefaultBufferSize;
    }

    if (param.comp.bufsz > MaxBufferSize) {
      Logger::Warning("Too big bufsz setting, reset to 256M");
      param.comp.bufsz = DefaultBufferSize;
    }

    if (param.comp.preset > 9) {
      Logger::Warning("Invalid preset setting, reset to 6");
      param.comp.preset = 6;
    }
  } else {
    param.mode = MainParam::Mode::DECODE;

    char cwd[PATH_MAX + 1];
    if (!getcwd(cwd, sizeof(cwd)))
      Logger::Error("Cannot get current work directory");
    param.decomp.outdir = cwd;

    const char* optstr = "o:t:gvh";
    static struct option options[] = {{"outdir", required_argument, NULL, 'o'},
                                      {"threads", optional_argument, NULL, 't'},
                                      {"gzip", no_argument, NULL, 'g'},
                                      {"verbose", no_argument, NULL, 'v'},
                                      {"help", no_argument, NULL, 'h'}};
    int c;
    int optidx;
    optind++;
    while (true) {
      c = getopt_long(argc, argv, optstr, options, &optidx);
      if (c == -1) break;
      switch (c) {
        case 0:
          break;
        case 'o':
          param.decomp.outdir = optarg;
          break;
        case 't':
          param.decomp.threads = std::stoull(optarg);
          break;
        case 'g':
          param.decomp.usegz = true;
          break;
        case 'v':
          Logger::verbose = true;
          break;
        case 'h':
          usage();
          break;
        default:
          break;
      }
    }

    param.decomp.infn = argv[optind];

    char actualpath[PATH_MAX + 1];
    if (!realpath(param.decomp.infn.c_str(), actualpath))
      Logger::Error("Invalid input file " + param.decomp.infn);

    if (!realpath(param.decomp.outdir.c_str(), actualpath))
      Logger::Error("Invalid output directory " + param.decomp.outdir);
  }
  return true;
}
