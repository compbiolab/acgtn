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

#include <zlib.h>
#include <algorithm>
#include <deque>

#include "index.h"
#include "kmer.h"
#include "kseq.h"
#include "logger.h"

KSEQ_INIT(gzFile, gzread)

static uchar* SEQ = nullptr;
static Index* INDEX = nullptr;
static tsl::robin_map<uint64, uchar>* SNPMAP = nullptr;

static const std::vector<std::string> int2snp{
    "",  "A",  "C",  "AC",  "G",  "AG",  "CG",  "ACG",
    "T", "AT", "CT", "ACT", "GT", "AGT", "CGT", "ACGT"};

static const tsl::robin_map<std::string, uint8> snp2int{
    {"A", 1},    {"C", 2},   {"AC", 3},   {"G", 4},    {"AG", 5},
    {"CG", 6},   {"ACG", 7}, {"T", 8},    {"AT", 9},   {"CT", 10},
    {"ACT", 11}, {"GT", 12}, {"AGT", 13}, {"CGT", 14}, {"ACGT", 15}};

struct sortSAWorkerParam {
  uint64* cnts;
  uint64 k;
  uint64 thd;
  uint64 beg;
  uint64 end;
};

struct kmerMapWorkerParam {
  uint64 thd;
  uint64 beg;
  uint64 end;
};

struct estiNSAWorkerParam {
  uint64 k;
  uint64 thd;
  uint64 beg;
  uint64 end;
  uint64* tnsa;
  uint64* cnts;
};

inline int seqCompartor(const uint64* vseq, const uint64* vmid,
                        const uint64 step);

inline void getseq(const uint64 pos, uchar seq[56]);

Index::Index(const IndexParam& param)
    : param_(param), seq_(nullptr), kfwd_(nullptr), gsz_(0), nsa_(0) {}

Index::~Index() {
  FreeArrPtr(seq_);
  FreePtr(kfwd_);
}

int Index::build() {
  INDEX = this;
  SNPMAP = &snpmap_;

  tsl::robin_map<std::string, Range> chrmap;
  loadGenome(chrmap);
  loadSNP(chrmap);
  computeSA();
  computeKmerMap();

  const std::string outfn = param_.output + ".idx";
  std::ofstream out(outfn, std::ios::binary);
  if (!out) {
    std::cerr << "Error: Failed to open file `" << outfn << "`" << std::endl;
    exit(EXIT_FAILURE);
  }

  write_member(gsz_, out);
  out.write((char*)seq_, gsz_ + 1);
  kfwd_->serialize(out);
  sa_.serialize(out);
  kmap_.serialize(out);

  if (snpmap_.size()) {
    intvector<> snps;
    snps.set_size_width(snpmap_.size(), PosBits + 2);
    uint64 i = 0;
    for (const auto& snp : snpmap_)
      snps[i++] = ((Kmer::base2id[snp.second]) << PosBits) | (uint64)snp.first;
    snps.serialize(out);
  }

  out.close();

  /*
  uchar tmpseq[L];
  for (uint64 i = 0; i < nsa_; i++) {
    const uint64 pos = sa_[i];
    const uchar* seq = pos2seq(pos, tmpseq);
    std::cout << std::string((char*)seq, L) << '\t' << (pos & PosMask) << '\t'
              << (pos >> PosBits) << std::endl;
  }
  */

  RefSize = gsz_;
  Logger::Info("Finish build index");
  return 0;
}

bool Index::load(const std::string file) {
  std::ifstream in(file, std::ios::binary);
  if (!in) {
    std::cerr << "Error: Failed to open file `" << file << "`" << std::endl;
    exit(EXIT_FAILURE);
  }

  read_member(gsz_, in);
  seq_ = new uchar[gsz_ + 1];
  in.read((char*)seq_, gsz_ + 1);
  const uint64 nkmer = 1ULL << (K * 2);
  kfwd_ = new bitarray(nkmer);
  kfwd_->load(in);
  sa_.load(in);
  kmap_.load(in);

  nsa_ = sa_.size();
  SeqBits = sa_.width();
  PosBits = ld(gsz_) + 1;
  PosMask = BitMask[PosBits];
  assert(SeqBits == PosBits || SeqBits == PosBits + 1);

  if (SeqBits != PosBits) {
    intvector<> snps;
    snps.load(in);
    assert(snps.width() == PosBits + 2);

    for (uint64 i = 0; i < snps.size(); ++i) {
      const uint64 pos = snps[i] & PosMask;
      const uchar snp = int2base(snps[i] >> PosBits);
      assert(Kmer::isbase(snp));
      snpmap_.emplace(std::make_pair(pos, snp));
    }
  }

  in.close();
  /*
  for (uint64 i = 0; i < kmermap_.size(); i++) {
    if (!(kmermap_[i] & absentmask_)) {
      std::cout << Kmer::toStr(i, k_) << '\t' << kmermap_[i] << std::endl;
    }
  }
  */

  RefSize = gsz_;
  Logger::Info("Finish load index");
  return true;
}

bool Index::encode(Dnaseq& read) {
  read.isfwd = true;
  if (read.len < 68) return false;
  uint64 ntry = MIN(read.ntry, read.len - 67);
  read.idx = 0;
  for (uint64 i = 0; i < ntry; ++i) {
    if (extendSeq(read)) return true;
    read.idx++;
  }

  Kmer::rcSeq(read.dna, read.len);
  read.isfwd = false;
  read.idx = 0;
  for (uint64 i = 0; i < ntry; ++i) {
    if (extendSeq(read)) return true;
    read.idx++;
  }
  /*
  for (read.idx = 0; read.idx < 3; ++read.idx) {
    const uchar* seq = read.dna + E - read.idx * 9;
    const uint64 begid = Kmer::revId(seq - K, K);

    const bool hasfwd = kfwd_->test(begid);
    const bool hasbwd = kbwd_->test(Kmer::rcId(begid, K));

    read.isfwd = true;
    if (!(hasfwd || hasbwd)) {
      continue;
    } else if (hasfwd && !hasbwd) {
      if (extendSeq(read)) return true;
    } else if (!hasfwd && hasbwd) {
      Kmer::rcSeq(read.dna, read.len);
      read.isfwd = false;
      if (extendSeq(read)) return true;
    } else {
      if (extendSeq(read)) return true;
      Kmer::rcSeq(read.dna, read.len);
      read.isfwd = false;
      if (extendSeq(read)) return true;
    }
  }
  */
  return false;
}

bool Index::extendSeq(Dnaseq& read) {
  bool success = false;
  uint64 skip = read.len / 10;
  uint64 extsz = 0;
  uint64 ntry = MIN(read.ntry, read.len - 67);
  if (ntry > 1) extsz = (read.len - L - 2 * skip) / (ntry - 1);
  const uchar* seq = read.dna + read.idx * extsz + skip;
  const uint64 kmerid = Kmer::revId(seq + L - K, K);
  if (!kfwd_->test(kmerid)) return false;
  uint64 low = kmap_[kmerid];
  uint64 high = kmap_[kmerid + 1];
  assert(high > low);
  uint64 count = high - low;
  uint64 curr, step;

  const uint64* vseq = (uint64*)(seq + L - K - 8);
  // MIN(E / B - read.idx * 2, (read.len - read.idx * 9 - K - 12) / B);
  // if (read.step < step) step = read.step;
  uchar tmpseq[L];
  while (count) {
    step = count / 2;
    curr = low + step;
    /*
    std::fprintf(stderr, "idx=%lu, low=%lu, high=%lu, curr=%lu\n", read.idx,
                 low, high, curr);
    */
    const uchar* currseq = pos2seq(sa_[curr], tmpseq);
    uint64 nerr = 0;
    for (uint64 i = 0; i < L; ++i)
      if (seq[i] != currseq[i]) ++nerr;
    if (nerr < 5) {
      success = true;
      break;
    }
    const uint64* vcurr = (uint64*)(currseq + L - K - 8);
    const int result = seqCompartor(vseq, vcurr, N);
    if (result > 0) {
      low = curr + 1;
      count -= step + 1;
    } else if (result < 0) {
      // high = curr;
      count = step;
    } else {
      success = true;
      break;
    }
  }

  /*
  do {
    assert(low < high);
    mid = average(low, high);
    std::fprintf(stderr, "low=%lu, high=%lu, mid=%lu, ncomp=%lu\n", low, high,
                 mid, ncomp);
    const uchar* midseq = pos2seq(sa_[mid], tmpseq);
    const uint64* vmid = (uint64*)(midseq + E - B);
    const int result = seqCompartor(vseq, vmid, ncomp);
    if (result > 0) {
      low = mid;
    } else if (result < 0) {
      high = mid;
    } else {
      success = true;
    }
    if (low + 1 >= high) {
      // if (!seqCompartor(vseq, vmid, 1)) success = true;
      break;
    }
  } while (!success);
  */

  if (!success) return false;
  // read.pos = (sa_[mid] & PosMask) + read.cur + L - read.len;
  read.pos = (sa_[curr] & PosMask) - read.idx * extsz - skip;
  return true;
}

inline uchar* Index::pos2seq(const uint64 pos, uchar seq[L]) const {
  const uint64 n = pos >> PosBits;
  if (!n) return seq_ + pos;

  for (uint64 i = 0, p = pos & PosMask; i < L; ++i) {
    if (snpmap_.find(p) != snpmap_.end()) {
      seq[i] = snpmap_.at(p++);
    } else {
      seq[i] = seq_[p++];
    }
  }
  return seq;
}

void Index::loadGenome(tsl::robin_map<std::string, Range>& chrmap) {
  const std::vector<std::string> files = param_.gfns;
  Array<uchar> seq;
  Range range;
  for (const auto& fn : files) {
    gzFile fp = gzopen(fn.c_str(), "rb");
    kseq_t* ks = kseq_init(fp);
    int l = -1;
    while ((l = kseq_read(ks)) >= 0) {
      if (ks->seq.l == 0) continue;
      std::string chrom(ks->name.s, ks->name.l);
      assert(chrom != "");
      range.first = seq.size();
      for (uint64 i = 0; i < ks->seq.l; ++i) seq.push_back(*(ks->seq.s + i));
      range.second = seq.size();
      chrmap.emplace(std::make_pair(chrom, range));
    }
    kseq_destroy(ks);
    gzclose(fp);
  }

  for (uint64 i = 0; i < seq.size(); ++i) seq[i] = std::toupper(seq[i]);
  gsz_ = seq.size();
  seq_ = new uchar[gsz_ + 1];
  memcpy(seq_, seq.data(), seq.size());
  seq_[gsz_] = '\0';
  SEQ = seq_;

  PosBits = ld(gsz_) + 1;
  SeqBits = PosBits;
  assert(SeqBits <= 40);
  PosMask = BitMask[PosBits];
  Logger::Info("Finish load genome");
}

void Index::loadSNP(const tsl::robin_map<std::string, Range>& chrmap) {
  const std::string snpfn = param_.snpfn;
  if (snpfn == "") return;

  std::ifstream in(snpfn, std::ios::binary);
  if (!in) {
    std::cerr << "Error: Failed to open file `" << snpfn << "`" << std::endl;
    exit(EXIT_FAILURE);
  }

  uint64 pos;
  uchar allele;
  std::string rsid, chrom;
  while (in >> rsid >> chrom >> pos >> allele) {
    if (chrmap.find(chrom) == chrmap.end()) continue;
    pos += chrmap.at(chrom).first;
    assert(pos < chrmap.at(chrom).second && Kmer::isbase(allele));
    snpmap_.emplace(std::make_pair(pos, allele));
  }

  if (snpmap_.size()) SeqBits = PosBits + 1;
  assert(SeqBits <= 40);
  in.close();
  Logger::Info("Finish load snp");
}

void Index::computeSA() {
  const uint64 k = 8;
  const uint64 bucketN = 1ULL << (2 * k);
  uint64* bucketCnts = new uint64[bucketN]();
  nsa_ = estimateNSA(bucketCnts, k);

  const uint64 thds = param_.threads;
  assert(thds > 1);
  std::vector<Range> ranges(thds);
  const uint64 maxCnts = nsa_ / (thds - 1);
  uint64* cnts = new uint64[thds];

  uint64 currCnt = 0;
  ranges[0].first = 0;
  for (uint64 i = 0, thd = 0; i < bucketN; ++i) {
    const uint64 bcnt = bucketCnts[i];
    if ((currCnt + bcnt > maxCnts) && (thd + 1 < thds)) {
      cnts[thd] = currCnt;
      ranges[thd].second = i;
      ranges[++thd].first = i;
      currCnt = 0;
    }
    currCnt += bcnt;
  }
  cnts[thds - 1] = currCnt;
  ranges[thds - 1].second = bucketN;
  delete[] bucketCnts;

  std::vector<tthread::thread*> threads(thds);
  for (uint64 thd = 0; thd < thds; ++thd) {
    uint64 beg = ranges[thd].first;
    uint64 end = ranges[thd].second;

    sortSAWorkerParam* param = new sortSAWorkerParam;
    param->cnts = cnts;
    param->k = k;
    param->thd = thd;
    param->beg = beg;
    param->end = end;

    threads[thd] = new tthread::thread(sortSAWorker, (void*)param);
  }

  for (uint64 thd = 0; thd < thds; ++thd) {
    tthread::thread* t = threads[thd];
    t->join();
    delete t;
  }

  nsa_ = 0;
  for (uint64 thd = 0; thd < thds; thd++) nsa_ += cnts[thd];
  sa_.set_size_width(nsa_, SeqBits);
  delete[] cnts;

  for (uint64 thd = 0, n = 0; thd < thds; ++thd) {
    const std::string infn = param_.output + std::to_string(thd) + ".sa";
    std::ifstream in(infn, std::ios::binary);
    if (!in) {
      std::cerr << "Error: Failed to open file `" << infn << "`" << std::endl;
      exit(EXIT_FAILURE);
    }
    intvector<> bsa;
    bsa.load(in);
    in.close();
    const uint64 size = bsa.size();
    for (uint64 i = 0; i < size; i++) {
      assert(n < nsa_);
      sa_[n++] = (uint64)bsa[i];
    }
    std::remove(infn.c_str());
  }
  Logger::Info("Finish compute SA");
}

void Index::sortSAWorker(void* param) {
  sortSAWorkerParam* param0 = (sortSAWorkerParam*)param;
  const tsl::robin_map<uint64, uchar>& snpmap = INDEX->snpmap_;
  const uchar* const seq = INDEX->refseq();

  const uint64 k = param0->k;
  const uint64 thd = param0->thd;
  const uint64 beg = param0->beg;
  const uint64 end = param0->end;
  uint64 cnts = param0->cnts[thd];
  const uint64 gsz = INDEX->gsz_;
  const std::string output = INDEX->param_.output;

  /*
  std::string info = "sortSABucket";
  info += std::to_string(thd);
  RUNINFO(info);
  */

  uintx* saBucket = new uintx[cnts]();
  if (!saBucket) {
    std::cerr << "Error: Failed to allocate memory" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::deque<uchar> snpque;
  uint64 snpidx = -1;
  uchar tmpseq[L];
  for (uint64 p = 0, j = 0; p < gsz - L - 1; ++p) {
    if (!Kmer::hasN(seq + p, L)) {
      uint64 bucket = Kmer::revId(seq + p + L - k, k);
      if (bucket >= beg && bucket < end) saBucket[j++] = p;
    }

    if (p == 0) {
      for (uint64 i = 0; i < L; ++i) {
        if (snpmap.find(i) != snpmap.end()) {
          assert(Kmer::isbase(snpmap.at(p + i)) && snpidx == -1);
          snpque.push_back(snpmap.at(i));
          snpidx = i;
        } else {
          snpque.push_back('\0');
        }
      }
    } else {
      if (snpque.front() != '\0') snpidx = -1;
      snpque.pop_front();
      if (snpmap.find(p + L - 1) != snpmap.end()) {
        assert(Kmer::isbase(snpmap.at(p + L - 1)) && snpidx == -1);
        snpque.push_back(snpmap.at(p + L - 1));
        snpidx = L;
      } else {
        snpque.push_back('\0');
      }
      if (snpidx != -1) --snpidx;
    }

    assert(snpque.size() == L && (snpidx == -1 || snpidx < L));
    if (snpidx == -1) continue;
    for (uint64 i = 0; i < L; i++) tmpseq[i] = seq[p + i];
    assert(snpque[snpidx]);
    tmpseq[snpidx] = snpque[snpidx];
    if (!Kmer::hasN(tmpseq, L)) {
      const uint64 bucket = Kmer::revId(tmpseq + L - k, k);
      if (bucket >= beg && bucket < end)
        saBucket[j++] = p | (1ULL << INDEX->PosBits);
    }
  }

  auto sufffixComparator = [](const void* a, const void* b) -> int {
    uchar seqa[56], seqb[56];
    getseq(*(uintx*)a, seqa);
    getseq(*(uintx*)b, seqb);
    const uint64* ga = (uint64*)(seqa + 48);
    const uint64* gb = (uint64*)(seqb + 48);

    uint64 va = 0, vb = 0;
    for (uint64 i = 0; i < 7; ++i) {
      va = *(ga - i);
      vb = *(gb - i);

      if (va > vb) {
        return 1;
      } else if (va < vb) {
        return -1;
      };
    }

    return 0;
  };

  auto sufffixEqual = [&sufffixComparator](const uintx a,
                                           const uintx b) -> bool {
    return sufffixComparator((void*)&a, (void*)&b) == 0;
  };

  qsort(saBucket, cnts, sizeof(uintx), sufffixComparator);
  uintx* last = std::unique(saBucket, saBucket + cnts, sufffixEqual);
  cnts = std::distance(saBucket, last) - 1;
  *(param0->cnts + thd) = cnts;

  intvector<> sa;
  sa.set_size_width(cnts, INDEX->SeqBits);
  for (uint64 i = 0; i < cnts; ++i) {
    sa[i] = (uint64)(saBucket[i]);
  }
  delete[] saBucket;

  const std::string outfn = output + std::to_string(thd) + ".sa";
  std::ofstream out(outfn, std::ios::binary);
  if (!out) {
    std::cerr << "Error: Failed to open file `" << outfn << "`" << std::endl;
    exit(EXIT_FAILURE);
  }
  sa.serialize(out);
  out.close();
  delete param0;
}

void Index::computeKmerMap() {
  const uint64 nkmer = 1ULL << (K * 2);
  kfwd_ = new bitarray(nkmer);

  const uint64 nbits = ld(nsa_) + 1;
  kmap_.set_size_width(nkmer + 1, nbits);

  uchar tmpseq[L];
  const uint64 thds = param_.threads;
  auto bounds = jobBounds(nsa_, thds);
  uint64 prevEnd = 0;
  for (uint64 i = 0; i < thds - 1; ++i) {
    if (i) bounds[i].first = prevEnd;
    uint64 currEnd = bounds[i].second;
    const uchar* seq = pos2seq(sa_[currEnd], tmpseq);
    const uint64 currId = Kmer::revId(seq + L - K, K);
    for (uint64 j = currEnd + 1; j < bounds[i + 1].second; ++j) {
      ++currEnd;
      uchar tmpseq[L];
      const uchar* seq = pos2seq(sa_[j], tmpseq);
      const uint64 nextId = Kmer::revId(seq + L - K, K);
      assert(currId <= nextId);
      if (currId < nextId) break;
    }
    assert(currEnd >= prevEnd);
    bounds[i].second = currEnd;
    prevEnd = currEnd;
  }
  assert(nsa_ >= prevEnd);
  bounds[thds - 1].first = prevEnd;

  std::vector<tthread::thread*> threads(thds);
  for (uint64 thd = 0; thd < thds; ++thd) {
    const uint64 beg = bounds[thd].first;
    const uint64 end = bounds[thd].second;

    kmerMapWorkerParam* param = new kmerMapWorkerParam;
    param->thd = thd;
    param->beg = beg;
    param->end = end;

    threads[thd] = new tthread::thread(kmerMapWorker, (void*)param);
  }

  for (uint64 thd = 0; thd < thds; ++thd) {
    tthread::thread* t = threads[thd];
    t->join();
    delete t;
  }
  Logger::Info("Finish compute kmer map");
}

void Index::kmerMapWorker(void* param) {
  const kmerMapWorkerParam* param0 = (kmerMapWorkerParam*)param;
  const uint64 thd = param0->thd;
  const uint64 beg = param0->beg;
  const uint64 end = param0->end;

  /*
  std::string info = "kmerMapWorker";
  info += std::to_string(thd);
  RUNINFO(info);
  */

  intvector<>* sa = &INDEX->sa_;
  intvector<>* kmap = &INDEX->kmap_;
  bitarray* kfwd = INDEX->kfwd_;

  uchar tmpseq[L];
  uint64 prevId = -1;
  for (uint64 i = beg; i < end; ++i) {
    const uchar* seq = INDEX->pos2seq(sa->at(i), tmpseq);
    const uint64 currId = Kmer::revId(seq + L - K, K);

    assert(currId + 1 >= prevId + 1);
    if (currId == prevId) continue;
    kfwd->set(currId);
    kmap->at(currId) = i;
    if (prevId + 1 != currId && prevId != -1) kmap->at(prevId + 1) = i;
    prevId = currId;
  }

  kmap->at(prevId + 1) = end == INDEX->nsa_ ? end - 1 : end;
  delete param0;
}

uint64 Index::estimateNSA(uint64* bucketCnts, const uint64 k) const {
  const uint64 bucketN = 1ULL << (2 * k);
  const uint64 thds = param_.threads;
  uint64* tnsa = new uint64[thds];
  const auto bounds = jobBounds(gsz_ - L - 1, thds);

  std::vector<tthread::thread*> threads(thds);
  std::vector<uint64*> bcnts(thds);
  for (uint64 thd = 0; thd < thds; ++thd) {
    const uint64 beg = bounds[thd].first;
    const uint64 end = bounds[thd].second;

    uint64* cnts = new uint64[bucketN]();
    bcnts[thd] = cnts;

    estiNSAWorkerParam* param = new estiNSAWorkerParam;
    param->k = k;
    param->thd = thd;
    param->beg = beg;
    param->end = end;
    param->cnts = cnts;
    param->tnsa = tnsa;

    threads[thd] = new tthread::thread(estimateNSAWorker, (void*)param);
  }

  for (uint64 thd = 0; thd < thds; ++thd) {
    tthread::thread* t = threads[thd];
    t->join();
    delete t;
  }

  uint64 nsa = 0;
  for (uint64 thd = 0; thd < thds; ++thd) {
    nsa += tnsa[thd];
    const uint64* cnts = bcnts[thd];
    for (uint64 i = 0; i < bucketN; ++i) bucketCnts[i] += cnts[i];
    delete[] cnts;
  }

  delete[] tnsa;
  Logger::Info("Finish estimate # of SA");
  return nsa;
}

void Index::estimateNSAWorker(void* param) {
  estiNSAWorkerParam* param0 = (estiNSAWorkerParam*)param;
  const uint64 k = param0->k;
  const uint64 thd = param0->thd;
  const uint64 beg = param0->beg;
  const uint64 end = param0->end;
  uint64* tnsa = param0->tnsa;
  uint64* cnts = param0->cnts;
  const tsl::robin_map<uint64, uchar>& snpmap = INDEX->snpmap_;

  /*
  std::string info = "estimateNSAWorker";
  info += std::to_string(thd);
  RUNINFO(info);
  */

  const uchar* const seq = INDEX->seq_;

  uint64 nsa = 0, snpidx = -1;
  std::deque<uchar> snpque;
  uchar tmpseq[L];
  for (uint64 p = beg; p < end; ++p) {
    if (!Kmer::hasN(seq + p, L)) {
      const uint64 bucket = Kmer::revId(seq + p + L - k, k);
      cnts[bucket]++;
      ++nsa;
    }

    if (p == beg) {
      assert(p + L <= end);
      for (uint64 i = 0; i < L; ++i) {
        if (snpmap.find(p + i) != snpmap.end()) {
          assert(Kmer::isbase(snpmap.at(p + i)) && snpidx == -1);
          snpque.push_back(snpmap.at(p + i));
          snpidx = i;
        } else {
          snpque.push_back('\0');
        }
      }
    } else {
      if (snpque.front() != '\0') snpidx = -1;
      snpque.pop_front();
      if (snpmap.find(p + L - 1) != snpmap.end()) {
        assert(Kmer::isbase(snpmap.at(p + L - 1)) && snpidx == -1);
        snpque.push_back(snpmap.at(p + L - 1));
        snpidx = L;
      } else {
        snpque.push_back('\0');
      }
      if (snpidx != -1) --snpidx;
    }

    assert(snpque.size() == L && (snpidx == -1 || snpidx < L));
    if (snpidx == -1) continue;
    for (uint64 i = 0; i < L; ++i) tmpseq[i] = seq[p + i];
    assert(snpque[snpidx]);
    tmpseq[snpidx] = snpque[snpidx];
    if (!Kmer::hasN(tmpseq, L)) {
      const uint64 bucket = Kmer::revId(tmpseq + L - k, k);
      cnts[bucket]++;
      ++nsa;
    }
  }

  tnsa[thd] = nsa;
  delete param0;
}

void Index::printInfo() {
  const uint64 nkmer = 1ULL << (K * 2);
  std::vector<uint64> cnts(257, 0);
  for (uint64 i = 0; i < nkmer; ++i) {
    if (kfwd_->test(i)) {
      uint64 block = kmap_[i + 1] - kmap_[i];
      assert(block);
      if (block > 256) {
        cnts[256]++;
      } else {
        cnts[block - 1]++;
      }
    }
  }

  std::cout << "#kmer block" << std::endl;
  for (uint64 i = 0; i < 257; ++i)
    std::cout << i + 1 << '\t' << cnts[i] << std::endl;
}

inline int seqCompartor(const uint64* vseq, const uint64* vmid,
                        const uint64 step) {
  uint64 vs = 0, vm = 0;
  for (uint64 i = 0; i < step; ++i) {
    vs = *(vseq - i);
    vm = *(vmid - i);

    if (vs > vm) {
      return 1;
    } else if (vs < vm) {
      return -1;
    };
  }

  return 0;
};

inline void getseq(const uint64 pos, uchar seq[56]) {
  seq[0] = seq[1] = 0;
  const uint64 n = pos >> INDEX->PosBits;
  if (!n) {
    memcpy(seq + 2, SEQ + pos, 54);
    return;
  }

  for (uint64 i = 2, p = pos & INDEX->PosMask; i < 56; ++i) {
    if (SNPMAP->find(p) != SNPMAP->end()) {
      seq[i] = SNPMAP->at(p++);
    } else {
      seq[i] = SEQ[p++];
    }
  }
}
