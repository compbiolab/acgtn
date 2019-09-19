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

#include "datapool.h"
#include "dataqueue.h"
#include "decompressor.h"
#include "kmer.h"
#include "logger.h"
#include "myarray.h"

static rocksdb::WriteOptions DBWOPT = rocksdb::WriteOptions();

using FastqsPool = DataPool<Fastqs>;
using FastqsQueue = DataQueue<Fastqs>;

using BinaryPool = DataPool<Seq64>;
using BinaryQueue = DataQueue<Seq64>;

struct ReloadReadsParam {
  uint64 bufsz;
  rocksdb::DB* db;
  Header* header;
  FastqsPool* fqsPool;
  FastqsQueue* fqsQueue;
};

struct ReadInputParam {
  FILE* in;
  uint64 bufsz;
  BinaryPool* buffPool;
  BinaryQueue* buffQueue;
};

struct DecodeReadsParam {
  Roaring* isrc;
  BinaryPool* buffPool;
  BinaryQueue* buffQueue;
  FastqsPool* fqsPool;
  FastqsQueue* fqsQueue;
};

Decompressor::Decompressor() : header_(), db_(nullptr), dbdir_("") {}

Decompressor::Decompressor(const DecompParam& param)
    : param_(param), header_(), db_(nullptr), dbdir_("") {}

Decompressor::~Decompressor() {
  if (db_) delete db_;
}

int Decompressor::Run() {
  DBWOPT.disableWAL = true;
  ReorderReads();
  AssignReads();
  return 0;
}

void Decompressor::AssignReads() {
  assert(dbdir_ != "");
  OpenDb(dbdir_, true);
  assert(db_);

  Vstring fns = header_.inputs;
  std::vector<crc64> crcs;
  std::vector<FILE*> outs(fns.size(), nullptr);
  for (uint64 i = 0; i < fns.size(); ++i) {
    std::string fn = fns[i];
    assert(fn != "");
    auto pos = fn.find_first_of('.');
    assert(pos != fn.npos);
    std::string outfn = fn.substr(0, pos) + ".dna";
    FILE* out = fopen(outfn.c_str(), "wb");
    if (!out) Logger::Error("Failed to open " + outfn);
    outs[i] = out;
    crcs.push_back(crc64(crc64::cc[15]));
  }

  FastqsPool fqsPool(param_.threads);
  FastqsQueue fqsQueue(param_.threads, 1);

  ReloadReadsParam readerParam;
  readerParam.bufsz = header_.bufsz;
  readerParam.db = db_;
  readerParam.header = &header_;
  readerParam.fqsPool = &fqsPool;
  readerParam.fqsQueue = &fqsQueue;

  tthread::thread readerThread(ReloadReadsWorker, (void*)&readerParam);

  uint64 id = 0;
  Fastqs* fastqs = nullptr;
  while (fqsQueue.Pop(id, fastqs)) {
    Fastqs fqs;
    for (auto it = fastqs->begin(); it != fastqs->end(); ++it)
      fqs.push_back(*it);
    fqsPool.Release(fastqs);

    for (Fastq& fq : fqs) {
      assert(fq.key < fns.size());
      for (uint64 i = 0; i < fq.seq.size() - 1; ++i)
        crcs[fq.key].byte_in(fq.seq[i]);
      fwrite(fq.seq.data(), 1, fq.seq.size(), outs[fq.key]);
    }
  }

  readerThread.join();
  for (auto f : outs) fclose(f);
  for (uint64 i = 0; i < fns.size(); ++i) {
    if (header_.crc64s[i] != crcs[i].get_a())
      Logger::Warning("CRC64 check failed in " + fns[i]);
  }
  remove_dir(dbdir_.c_str());
}

void Decompressor::ReloadReadsWorker(void* param) {
  ReloadReadsParam* param0 = (ReloadReadsParam*)param;
  uint64 bufsz = param0->bufsz;
  rocksdb::DB* db = param0->db;
  Header* header = param0->header;
  FastqsPool* fqsPool = param0->fqsPool;
  FastqsQueue* fqsQueue = param0->fqsQueue;

  uint64 tmpsz = 0;
  uint64 blockId = 0;
  Fastqs* fastqs = nullptr;
  fastqs = fqsPool->Acquire();
  fastqs->clear();

  uint64 fid = 0, rid = 0;
  rocksdb::Iterator* it = db->NewIterator(rocksdb::ReadOptions());
  for (it->SeekToFirst(); it->Valid(); it->Next()) {
    assert(it->status().ok());
    Fastq fastq;
    GetFastq(it, fastq);
    fastq.seq.push_back('\n');
    if (rid == header->reads[fid]) ++fid;
    fastq.key = fid;
    if (tmpsz + fastq.seq.size() >= bufsz) {
      fqsQueue->Push(blockId++, fastqs);
      fastqs = fqsPool->Acquire();
      fastqs->clear();
      tmpsz = 0;
    }
    fastqs->push_back(fastq);
    tmpsz += fastq.seq.size();
    ++rid;
  }

  if (fastqs->size()) {
    fqsQueue->Push(blockId++, fastqs);
    fastqs = fqsPool->Acquire();
  }

  fqsPool->Release(fastqs);
  fqsQueue->Done();

  delete it;
}

void Decompressor::ReorderReads() {
  char cwd[PATH_MAX + 1];
  if (!getcwd(cwd, sizeof(cwd)))
    Logger::Error("Cannot get current work directory");
  std::string tmpdb = cwd;
  // tmpdb += "/tmp_acgtn_decode_XXXXXX";
  std::srand(std::time(nullptr));
  tmpdb += "/tmp_acgtn_decode_" + std::to_string(std::rand());
  // char* tmpdir = mkdtemp((char*)tmpdb.c_str());
  // if (!tmpdir) Logger::Error("Failed to create directory for fastq db");
  if (mkdir(tmpdb.c_str(), 0750) == -1)
    Logger::Error("Failed to create directory for fastq db");
  dbdir_ = tmpdb;
  OpenDb(dbdir_, false);

  std::string infn = param_.infn;
  FILE* in = fopen(infn.c_str(), "rb");
  if (!in) Logger::Error("Failed to open " + infn);

  header_.load(in);
  bufsz_ = header_.bufsz;
  header_.print(true);

  uint64 rcsz = 0;
  int64 offset = -sizeof(uint64);
  fseek(in, offset, SEEK_END);
  fread(&rcsz, sizeof(uint64), 1, in);
  offset = -(rcsz + sizeof(uint64));
  fseek(in, offset, SEEK_END);
  char* isrcIn = new char[rcsz];
  fread(isrcIn, 1, rcsz, in);
  isrc_ = Roaring::read(isrcIn);
  delete[] isrcIn;

  fseek(in, 0, SEEK_SET);
  header_.load(in);

  uint64 threads = param_.threads;
  BinaryPool buffPool(threads);
  BinaryQueue buffQueue(threads, 1);

  FastqsPool fqsPool(threads);
  FastqsQueue fqsQueue(threads, threads);

  ReadInputParam readerParam;
  readerParam.in = in;
  readerParam.bufsz = bufsz_;
  readerParam.buffPool = &buffPool;
  readerParam.buffQueue = &buffQueue;

  tthread::thread readerThread(ReadInputWorker, (void*)&readerParam);

  DecodeReadsParam decodeParam;
  decodeParam.isrc = &isrc_;
  decodeParam.buffPool = &buffPool;
  decodeParam.buffQueue = &buffQueue;
  decodeParam.fqsPool = &fqsPool;
  decodeParam.fqsQueue = &fqsQueue;

  std::vector<tthread::thread*> decodeThreads(threads);
  for (uint64 i = 0; i < threads; ++i)
    decodeThreads[i] =
        new tthread::thread(DecodeReadsWorker, (void*)&decodeParam);

  uint64 buffId = 0;
  Fastqs* fastqs = nullptr;

  while (fqsQueue.Pop(buffId, fastqs)) {
    Fastqs fqs;
    for (auto it = fastqs->begin(); it != fastqs->end(); ++it)
      fqs.push_back(*it);
    fqsPool.Release(fastqs);
    rocksdb::WriteBatch batch;
    for (Fastq& fq : fqs) PutFastq(&batch, fq);
    rocksdb::Status s = db_->Write(DBWOPT, &batch);
    if (!s.ok()) Logger::Error("Write fastq batch");
  }

  readerThread.join();

  for (uint64 i = 0; i < threads; ++i) {
    tthread::thread* t = decodeThreads[i];
    t->join();
    delete t;
  }

  fclose(in);

  db_->Flush(rocksdb::FlushOptions());
  delete db_;
  db_ = nullptr;

  Logger::Info("Finish reorder reads");
}

void Decompressor::ReadInputWorker(void* param) {
  ReadInputParam* param0 = (ReadInputParam*)param;
  FILE* in = param0->in;
  uint64 bufsz = param0->bufsz;
  BinaryPool* buffPool = param0->buffPool;
  BinaryQueue* buffQueue = param0->buffQueue;

  uint64 id = 0;
  Seq64* buff = buffPool->Acquire(bufsz);
  buff->clear();

  while (true) {
    uint64 insz, ret;
    ret = fread(&insz, sizeof(uint64), 1, in);
    if (!insz) {
      buffPool->Release(buff);
      buffQueue->Done();
      Logger::Info("Finish read worker");
      return;
    }
    buff->resize(insz);
    ret = fread(buff->data(), 1, insz, in);
    if (ret != insz) Logger::Error("Read input data");
    buffQueue->Push(id++, buff);
    buff = buffPool->Acquire(bufsz);
    buff->clear();
  }
}

void Decompressor::DecodeReadsWorker(void* param) {
  DecodeReadsParam* param0 = (DecodeReadsParam*)param;
  Roaring* isrc = param0->isrc;
  BinaryPool* buffPool = param0->buffPool;
  BinaryQueue* buffQueue = param0->buffQueue;
  FastqsPool* fqsPool = param0->fqsPool;
  FastqsQueue* fqsQueue = param0->fqsQueue;

  uint64 buffId = 0;
  Seq64 inBuff;
  Seq64* buff = nullptr;
  Fastqs* fastqs = nullptr;

  while (buffQueue->Pop(buffId, buff)) {
    inBuff.resize(buff->size());
    memcpy(inBuff.data(), buff->data(), buff->size());
    buffPool->Release(buff);

    uint64 tmpsz = 0;
    memcpy(&tmpsz, inBuff.data(), 8);
    uint64 idsz = tmpsz >> 32;
    uint64 fqnum = tmpsz & BitMask[32];
    uchar* inIds = new uchar[idsz];
    uint64* ids = new uint64[fqnum + 64];
    memcpy(inIds, inBuff.data() + 8, idsz);
    uint64 ret = p4ndec64(inIds, fqnum, ids);
    // uint64 ret = p4ndec64(inBuff.data() + 8, fqnum, ids);
    assert(ret == idsz);
    delete[] inIds;
    tmpsz = idsz + 8;

    Seq32 rawseq;
    uint32 rawsz = 0;
    memcpy(&rawsz, inBuff.data() + tmpsz, 4);
    assert(rawsz);
    rawseq.resize(rawsz);
    assert(tmpsz + 4 + rawsz == inBuff.size());
    memcpy(rawseq.data(), inBuff.data() + tmpsz + 4, rawsz);

    Seq32 seqs;
    DecodeLzma(rawseq, seqs);
    assert(seqs.size() < BitMask[32]);

    fastqs = fqsPool->Acquire();
    fastqs->clear();

    assert(fqnum == std::count(seqs.begin(), seqs.end(), '\n'));
    Seq32 seq;
    for (uint64 i = 0, j = 0; i < seqs.size(); ++i) {
      uchar base = seqs[i];
      if (base == '\n') {
        uint64 rid = ids[j];
        assert(rid < BitMask[40] && seq.size());
        if (isrc->contains(rid)) Kmer::rcSeq(seq);
        Fastq fq;
        fq.key = rid;
        fq.seq = seq;
        fastqs->push_back(fq);
        seq.clear();
        ++j;
      } else {
        seq.push_back(seqs[i]);
      }
    }
    delete[] ids;

    fqsQueue->Push(buffId, fastqs);
  }

  fqsQueue->Done();
  Logger::Info("Finish decode worker");
}

void Decompressor::PutFastq(rocksdb::WriteBatch* batch, const Fastq& fastq) {
  uint64 rid = fastq.key;

  std::string key;
  key.resize(5);

  key[0] = (rid >> 32) & 0xFF;
  key[1] = (rid >> 24) & 0xFF;
  key[2] = (rid >> 16) & 0xFF;
  key[3] = (rid >> 8) & 0xFF;
  key[4] = rid & 0xFF;

  std::string value = std::string((char*)fastq.seq.data(), fastq.seq.size());

  batch->Put(key, value);
}

void Decompressor::GetFastq(rocksdb::Iterator* it, Fastq& fastq) {
  const std::string value = it->value().ToString();
  fastq.seq.clear();
  for (const uchar c : value) fastq.seq.push_back(c);
}

rocksdb::Options Decompressor::GetOptions() {
  rocksdb::Options options;

  options.max_open_files = -1;
  options.max_background_jobs = param_.threads;
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

void Decompressor::OpenDb(const std::string& dir, bool readOnly) {
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
