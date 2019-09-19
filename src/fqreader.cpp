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

#include "fqreader.h"
#include "logger.h"

Fqreader::Fqreader(const std::string fn, uint64 bufsz)
    : bufsz_(bufsz), idx_(0), rid_(0), fp_(nullptr), ks_(nullptr) {
  fns_.push_back(fn);

#define initHeader                          \
  hdr_.bufsz = bufsz;                       \
  hdr_.reads.clear();                       \
  hdr_.inputs.clear();                      \
                                            \
  for (auto& fn : fns_) {                   \
    const auto pos = fn.find_last_of('/');  \
    std::string input = fn.substr(pos + 1); \
    assert(input != "");                    \
                                            \
    hdr_.reads.push_back(0);                \
    hdr_.crc64s.push_back(0);               \
    hdr_.inputs.push_back(input);           \
    crcs_.push_back(crc64(crc64::cc[15]));  \
  }

  initHeader;
  assert(OpenNext());
}

Fqreader::Fqreader(const Vstring& fns, uint64 bufsz)
    : fns_(fns), bufsz_(bufsz), idx_(0), rid_(0), fp_(nullptr), ks_(nullptr) {
  initHeader;
  assert(OpenNext());
}

#undef initHeader

Fqreader::~Fqreader() { CloseNext(); }

bool Fqreader::ReadFastqs(Fastqs& fqs) {
  fqs.clear();
  uint64 bufsz = 0;
  while (bufsz < bufsz_) {
    int ret = kseq_read(ks_);
    if (ret == -1) {
      hdr_.crc64s[idx_ - 1] = crcs_[idx_ - 1].get_a();
      if (!OpenNext()) return fqs.size() ? true : false;
      ret = kseq_read(ks_);
    }
    assert(idx_ && (rid_ < (1ULL << 40)));
    if (ret == -2) Logger::Error("Bad fastq file " + fns_[idx_ - 1]);
    if (ks_->seq.l == 0) continue;
    Fastq fq;
    fq.key = rid_;
    uint64 len = ks_->seq.l;
    fq.seq.resize(len);
    memcpy(fq.seq.data(), ks_->seq.s, len);
    bufsz += len;
    fqs.push_back(fq);
    ++rid_;
    if (len < hdr_.minlen) hdr_.minlen = len;
    if (len > hdr_.maxlen) hdr_.maxlen = len;
    hdr_.bases += len;
    for (uint64 i = 0; i < len; ++i) crcs_[idx_ - 1].byte_in(fq.seq[i]);
  }
  return true;
}

void Fqreader::CloseNext() {
  if (!idx_ && !fp_) return;

  if (fp_) gzclose(fp_);
  if (ks_) kseq_destroy(ks_);
  fp_ = nullptr;
  ks_ = nullptr;
}

bool Fqreader::OpenNext() {
  if (idx_) hdr_.reads[idx_ - 1] = rid_;
  if (idx_ >= fns_.size()) return false;
  CloseNext();

  const auto fn = fns_[idx_];
  if (!(fp_ = gzopen(fn.c_str(), "rb"))) Logger::Error("Failed to open " + fn);
  ks_ = kseq_init(fp_);

  ++idx_;
  return true;
}
