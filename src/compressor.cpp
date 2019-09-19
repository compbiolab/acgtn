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

#include <rocksdb/table.h>
#include "vp4.h"

#include "compressor.h"
#include "datapool.h"
#include "dataqueue.h"
#include "kmer.h"
#include "logger.h"
#include "strmcoder.h"

static rocksdb::WriteOptions DBWOPT = rocksdb::WriteOptions();

using FastqsPool = DataPool<Fastqs>;
using FastqsQueue = DataQueue<Fastqs>;

using BlockPool = DataPool<Seq64>;
using BlockQueue = DataQueue<Seq64>;

struct AssignKeysParam {
  uint64 beg;
  uint64 end;
  Fastqs* fqs;
  uint64* keys;
  uint64 ntry;
  Index* index;
};

struct ReadFastqParam {
  Fqreader* reader;
  FastqsPool* fqsPool;
  FastqsQueue* fqsQueue;
};

struct ReloadReadsParam {
  uint64 bufsz;
  rocksdb::DB* db;
  Roaring* isrc;
  FastqsPool* fqsPool;
  FastqsQueue* fqsQueue;
};

struct EncodeReadsParam {
  uint64 bufsz;
  uint64 preset;
  FastqsPool* fqsPool;
  FastqsQueue* fqsQueue;
  BlockPool* outPool;
  BlockQueue* outQueue;
};

struct Compressor::FastqInfo {
  Fastq fastq;
  uint64 pos;
  uint64 rid;
  bool isfwd;
};

Compressor::Compressor() : db_(nullptr), dbdir_("") {}

Compressor::Compressor(const CompParam& param)
    : param_(param), db_(nullptr), dbdir_("") {}

Compressor::~Compressor() {
  if (db_) delete db_;
}

int Compressor::Run() {
#ifndef NDEBUG
  nsuc_ = 0;
  nrds_ = 0;
  nfld_ = 0;
#endif

  DBWOPT.disableWAL = true;
  fns_ = param_.infns;
  ReorderReads();
  EncodeReads();
  return 0;
}

void Compressor::ReorderReads() {
  char cwd[PATH_MAX + 1];
  if (!getcwd(cwd, sizeof(cwd)))
    Logger::Error("Cannot get current work directory");
  std::string tmpdb = cwd;
  // tmpdb += "/tmp_acgtn_encode_XXXXXX";
  std::srand(std::time(nullptr));
  tmpdb += "/tmp_acgtn_encode_" + std::to_string(std::rand());
  // char* tmpdir = mkdtemp((char*)tmpdb.c_str());
  if (mkdir(tmpdb.c_str(), 0750) == -1)
    Logger::Error("Failed to create directory for fastq db");
  dbdir_ = tmpdb;
  OpenDb(dbdir_, false);

  Index index;
  index.load(param_.idxfn);

  Fqreader reader(param_.infns, param_.bufsz);
  FastqsPool fqsPool(param_.threads);
  FastqsQueue fqsQueue(param_.threads, 1);

  ReadFastqParam readerParam;
  readerParam.reader = &reader;
  readerParam.fqsPool = &fqsPool;
  readerParam.fqsQueue = &fqsQueue;

  tthread::thread readerThread(ReadInputsWorker, (void*)&readerParam);

  uint64 id = 0;
  Fastqs* fastqs = nullptr;

  FastqInfo info;
  while (fqsQueue.Pop(id, fastqs)) {
    Fastqs fqs;
    for (auto it = fastqs->begin(); it != fastqs->end(); ++it)
      fqs.push_back(*it);
    fqsPool.Release(fastqs);

    uint64 nfq = fqs.size();
    uint64* keys = new uint64[nfq];

    uint64 thds = param_.threads;
    const auto bounds = jobBounds(nfq, thds);
    std::vector<tthread::thread*> akThds(thds);
    for (uint64 thd = 0; thd < thds; ++thd) {
      uint64 beg = bounds[thd].first;
      uint64 end = bounds[thd].second;

      AssignKeysParam* param = new AssignKeysParam;
      param->beg = beg;
      param->end = end;
      param->fqs = &fqs;
      param->keys = keys;
      param->ntry = param_.ntry;
      param->index = &index;

      akThds[thd] = new tthread::thread(AssignKeysWorker, (void*)param);
    }

    for (uint64 thd = 0; thd < thds; ++thd) {
      tthread::thread* t = akThds[thd];
      t->join();
      delete t;
    }

    rocksdb::WriteBatch batch;
    for (uint64 i = 0; i < nfq; ++i) {
      Fastq& fq = fqs[i];

      info.fastq = fq;
      info.rid = fq.key;
      info.isfwd = true;

      if (keys[i] == (uint64)-1) {
        info.pos = -1;
        PutFastq(&batch, info);
#ifndef NDEBUG
        ++nrds_;
#endif
        continue;
      }

      uint64 key = keys[i] & BitMask[40];
      bool success = key != BitMask[40];
      info.pos = success ? key : -2;
      if (success && (keys[i] >> 40)) {
        Kmer::rcSeq(fq.seq);
        info.fastq.seq = fq.seq;
        info.isfwd = false;
      }

      PutFastq(&batch, info);
#ifndef NDEBUG
      if (success) {
        ++nsuc_;
      } else {
        ++nfld_;
      }
#endif
    }

    delete[] keys;
    rocksdb::Status s = db_->Write(DBWOPT, &batch);
    if (!s.ok()) Logger::Error("Write fastq batch");
  }

  header_ = reader.GetHeader();
  readerThread.join();
  header_.print(true);

  db_->Flush(rocksdb::FlushOptions());
  delete db_;
  db_ = nullptr;

#ifndef NDEBUG
  std::cerr << "hasn=" << nrds_ << "\tsuccess=" << nsuc_ << "\tfailed=" << nfld_
            << std::endl;
#endif

  Logger::Info("Finish reorder reads");
}

void Compressor::ReadInputsWorker(void* param) {
  ReadFastqParam* param0 = (ReadFastqParam*)param;
  Fqreader* reader = param0->reader;
  FastqsPool* fqsPool = param0->fqsPool;
  FastqsQueue* fqsQueue = param0->fqsQueue;

  uint64 id = 0;
  Fastqs* fastqs = fqsPool->Acquire();

  while (reader->ReadFastqs(*fastqs)) {
    fqsQueue->Push(id++, fastqs);
    fastqs = fqsPool->Acquire();
  }

  fqsPool->Release(fastqs);
  fqsQueue->Done();
}

void Compressor::AssignKeysWorker(void* param) {
  AssignKeysParam* param0 = (AssignKeysParam*)param;
  uint64 beg = param0->beg;
  uint64 end = param0->end;
  Fastqs* fqs = param0->fqs;
  uint64* keys = param0->keys;
  Index* index = param0->index;

  Seq64 seq;
  Dnaseq read;
  read.ntry = param0->ntry;

  for (uint64 i = beg; i < end; ++i) {
    Fastq fq = fqs->at(i);
    if (Kmer::hasN(fq.seq.data(), fq.seq.size())) {
      keys[i] = -1;
      continue;
    }

    if (fq.seq.size() > seq.size()) seq.reserve(fq.seq.size());
    memcpy(seq.data(), fq.seq.data(), fq.seq.size());
    read.dna = seq.data();
    read.len = fq.seq.size();

    bool success = index->encode(read);
    uint64 key = read.pos;
    if (success && !read.isfwd) key |= (1ULL << 40);
    keys[i] = success ? key : BitMask[40];
  }
}

void Compressor::EncodeReads() {
  assert(dbdir_ != "");
  OpenDb(dbdir_, true);
  assert(db_);

  std::string fn = param_.outfn + ".acgtn";
  FILE* out = fopen(fn.c_str(), "wb");
  if (!out) Logger::Error("Failed to open " + fn);
  header_.serialize(out);

  uint64 threads = param_.threads;
  FastqsPool fqsPool(param_.threads);
  FastqsQueue fqsQueue(param_.threads, 1);

  ReloadReadsParam readerParam;
  readerParam.bufsz = param_.bufsz;
  readerParam.db = db_;
  readerParam.isrc = &isrc_;
  readerParam.fqsPool = &fqsPool;
  readerParam.fqsQueue = &fqsQueue;

  BlockPool outPool(threads);
  BlockQueue outQueue(threads, threads);

  EncodeReadsParam encodeParam;
  encodeParam.bufsz = param_.bufsz;
  encodeParam.preset = param_.preset;
  encodeParam.fqsPool = &fqsPool;
  encodeParam.fqsQueue = &fqsQueue;
  encodeParam.outPool = &outPool;
  encodeParam.outQueue = &outQueue;

  tthread::thread readerThread(ReloadReadsWorker, (void*)&readerParam);

  std::vector<tthread::thread*> encodeThreads(threads);
  for (uint64 i = 0; i < threads; ++i)
    encodeThreads[i] =
        new tthread::thread(EncodeReadsWorker, (void*)&encodeParam);

  uint64 blockId = 0;
  Seq64* buff = nullptr;

  while (outQueue.Pop(blockId, buff)) {
    uint64 outsz = buff->size();
    assert(outsz);
    fwrite(&outsz, sizeof(uint64), 1, out);
    fwrite(buff->data(), 1, buff->size(), out);
    outPool.Release(buff);
    buff = nullptr;
  }

  readerThread.join();

  for (uint64 i = 0; i < threads; ++i) {
    tthread::thread* t = encodeThreads[i];
    t->join();
    delete t;
  }

  uint64 endstrm = 0;
  fwrite(&endstrm, sizeof(uint64), 1, out);
  uint64 rcsz = isrc_.getSizeInBytes();
  char* isrcOut = new char[rcsz];
  isrc_.write(isrcOut);
  fwrite(isrcOut, 1, rcsz, out);
  fwrite(&rcsz, sizeof(uint64), 1, out);
  delete[] isrcOut;

  fclose(out);
  remove_dir(dbdir_.c_str());

  Logger::Info("Finish encode reads");
}

void Compressor::ReloadReadsWorker(void* param) {
  ReloadReadsParam* param0 = (ReloadReadsParam*)param;
  uint64 bufsz = param0->bufsz;
  rocksdb::DB* db = param0->db;
  Roaring* isrc = param0->isrc;
  FastqsPool* fqsPool = param0->fqsPool;
  FastqsQueue* fqsQueue = param0->fqsQueue;

  uint64 tmpsz = 0;
  uint64 blockId = 0;
  Fastqs* fastqs = nullptr;
  fastqs = fqsPool->Acquire();
  fastqs->clear();

  rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    assert(it->status().ok());
    FastqInfo info;
    GetFastq(it, info);
    if (!info.isfwd) isrc->add(info.fastq.key);
    if (tmpsz + info.fastq.seq.size() >= bufsz) {
      fqsQueue->Push(blockId++, fastqs);
      fastqs = fqsPool->Acquire();
      fastqs->clear();
      tmpsz = 0;
    }
    fastqs->push_back(info.fastq);
    tmpsz += info.fastq.seq.size();
  }

  if (fastqs->size()) {
    fqsQueue->Push(blockId++, fastqs);
    fastqs = fqsPool->Acquire();
  }

  fqsPool->Release(fastqs);
  fqsQueue->Done();

  isrc->runOptimize();
  isrc->setCopyOnWrite(true);

  delete it;
}

void Compressor::EncodeReadsWorker(void* param) {
  EncodeReadsParam* param0 = (EncodeReadsParam*)param;
  uint64 bufsz = param0->bufsz;
  uint64 preset = param0->preset;
  FastqsPool* fqsPool = param0->fqsPool;
  FastqsQueue* fqsQueue = param0->fqsQueue;
  BlockPool* outPool = param0->outPool;
  BlockQueue* outQueue = param0->outQueue;

  uint64 buffId = 0;
  Fastqs* fastqs = nullptr;
  Seq64* outBuff = nullptr;

  while (fqsQueue->Pop(buffId, fastqs)) {
    Seq32 seqs;
    uint64 fqnum = fastqs->size();
    uint64* ids = new uint64[fqnum];
    for (uint64 i = 0; i < fqnum; ++i) {
      const Fastq& fq = fastqs->at(i);
      ids[i] = fq.key;
      seqs += fq.seq;
      seqs.push_back('\n');
    }

    assert(seqs.size() < BitMask[32]);

    fqsPool->Release(fastqs);

    outBuff = outPool->Acquire(bufsz);
    outBuff->clear();

    uint64 idBound = fqnum * 5 * 5 / 3 + 1024;
    uchar* idOut = new uchar[idBound];
    if (!idOut) Logger::Error("Failed to allocate memory");
    uint64 idsz = p4nenc64(ids, fqnum, idOut);
    assert(idsz < BitMask[32]);

    uint64 tmpsz = (idsz << 32) | (fqnum & BitMask[32]);
    outBuff->resize(idsz + 8);
    memcpy(outBuff->data(), &tmpsz, 8);
    memcpy(outBuff->data() + 8, idOut, idsz);
    delete[] idOut;
    delete[] ids;

    Seq32 comps;
    EncodeLzma(seqs, preset, comps);
    assert(comps.size() < BitMask[32]);
    uint32 compsz = comps.size() & BitMask[32];
#ifndef NDEBUG
    std::cerr << "fqnum=" << fqnum << " idsz=" << idsz
              << " seqsz=" << seqs.size() << " compsz=" << compsz << std::endl;
#endif
    tmpsz = outBuff->size();
    outBuff->resize(tmpsz + compsz + 4);
    memcpy(outBuff->data() + tmpsz, &compsz, 4);
    memcpy(outBuff->data() + tmpsz + 4, comps.data(), compsz);

    outQueue->Push(buffId, outBuff);
  }

  outQueue->Done();
}

void Compressor::PutFastq(rocksdb::WriteBatch* batch,
                          const Compressor::FastqInfo& info) {
  uint64 pos = info.pos;
  uint64 rid = info.rid;
  assert(rid < BitMask[40]);

  std::string key;
  key.resize(10);

  if (pos == -1 || pos == -2) {
    assert(info.isfwd);
    key[0] = pos == -1 ? 0xFF : 0x7F;  // hasN:-1, unknown:-2
    key[1] = 0xFF;
    key[2] = 0xFF;
    key[3] = 0xFF;
    key[4] = 0xFE;  // keep fwd

  } else {
    pos = pos << 1;
    if (!info.isfwd) ++pos;  // fwd:0, bwd:1
    key[0] = (pos >> 32) & 0xFF;
    key[1] = (pos >> 24) & 0xFF;
    key[2] = (pos >> 16) & 0xFF;
    key[3] = (pos >> 8) & 0xFF;
    key[4] = pos & 0xFF;
  }

  key[5] = (rid >> 32) & 0xFF;
  key[6] = (rid >> 24) & 0xFF;
  key[7] = (rid >> 16) & 0xFF;
  key[8] = (rid >> 8) & 0xFF;
  key[9] = rid & 0xFF;

  std::string value =
      std::string((char*)info.fastq.seq.data(), info.fastq.seq.size());

  batch->Put(key, value);
}

void Compressor::GetFastq(rocksdb::Iterator* it, FastqInfo& info) {
  const std::string value = it->value().ToString();
  info.fastq.seq.clear();
  for (const uchar c : value) info.fastq.seq.push_back(c);

  uint64 rid = it->key()[5] & 0xFF;
  rid = (rid << 8) | (it->key()[6] & 0xFF);
  rid = (rid << 8) | (it->key()[7] & 0xFF);
  rid = (rid << 8) | (it->key()[8] & 0xFF);
  rid = (rid << 8) | (it->key()[9] & 0xFF);
  assert(rid < BitMask[40]);
  info.fastq.key = rid;

  info.isfwd = true;
  if (it->key()[4] & 0x1U) info.isfwd = false;
}

rocksdb::Options Compressor::GetOptions() {
  rocksdb::Options options;

  options.max_open_files = -1;
  options.max_background_jobs = param_.threads;
  options.max_background_flushes = param_.threads;
  options.max_background_compactions = param_.threads;
  options.max_subcompactions = param_.threads;
  options.allow_mmap_reads = true;

  options.IncreaseParallelism(param_.threads);

  options.create_if_missing = true;
  options.error_if_exists = false;
  options.db_log_dir = "/tmp";

  options.access_hint_on_compaction_start =
      rocksdb::Options::AccessHint::SEQUENTIAL;

  options.allow_concurrent_memtable_write = true;
  options.enable_write_thread_adaptive_yield = true;

  const size_t block_cache_bytes = 1 << 30;
  const size_t memtable_bytes = 4 * size_t(1 << 30);

  rocksdb::BlockBasedTableOptions topt;
  topt.format_version = 4;
  topt.block_size = 4 << 20;
  topt.block_cache = rocksdb::NewLRUCache(block_cache_bytes);
  options.table_factory.reset(NewBlockBasedTableFactory(topt));

  options.OptimizeLevelStyleCompaction(memtable_bytes);
  options.max_write_buffer_number = 4;
  options.write_buffer_size = (memtable_bytes / 4);
  options.min_write_buffer_number_to_merge = 1;

  options.num_levels = 10;
  options.compression = rocksdb::kLZ4Compression;
  options.level_compaction_dynamic_level_bytes = true;
  options.bottommost_compression = rocksdb::kZSTD;

  return options;
}

void Compressor::OpenDb(const std::string& dir, bool readOnly) {
  rocksdb::Options options = GetOptions();

  rocksdb::Status s;
  if (readOnly) {
    s = rocksdb::DB::OpenForReadOnly(options, dir, &db_);
  } else {
    s = rocksdb::DB::Open(options, dir, &db_);
  }

  if (!s.ok()) {
    if (db_) {
      delete db_;
    }
    db_ = nullptr;
    std::cerr << "Error: Can't open" << dir << std::endl;
    exit(EXIT_FAILURE);
  }
}
