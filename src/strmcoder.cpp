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

#include "strmcoder.h"

void EncodeLzma(const Seq32& in, int level, Seq32& out) {
  uint64 bound = lzma_stream_buffer_bound(in.size());
  if (out.cap() < bound) out.reserve(bound);
  size_t out_pos = 0;
  if (LZMA_OK != lzma_easy_buffer_encode(
                     level, LZMA_CHECK_CRC32, NULL,
                     reinterpret_cast<uint8_t*>(const_cast<uchar*>(in.data())),
                     in.size(), reinterpret_cast<uint8_t*>(out.data()),
                     &out_pos, out.cap()))
    abort();
  out.resize(out_pos);
}

void DecodeLzma(const Seq32& in, Seq32& out) {
  static const size_t kMemLimit = 1 << 30;  // 1 GB.
  lzma_stream strm = LZMA_STREAM_INIT;
  out.resize(8192);
  size_t result_used = 0;
  lzma_ret ret;
  ret = lzma_stream_decoder(&strm, kMemLimit, LZMA_CONCATENATED);
  if (ret != LZMA_OK) abort();
  size_t avail = out.size();
  strm.next_in = reinterpret_cast<const uint8_t*>(in.data());
  strm.avail_in = in.size();
  strm.next_out = reinterpret_cast<uint8_t*>(&out[0]);
  strm.avail_out = avail;
  while (true) {
    ret = lzma_code(&strm, strm.avail_in == 0 ? LZMA_FINISH : LZMA_RUN);
    if (ret == LZMA_STREAM_END) {
      result_used += avail - strm.avail_out;
      if (strm.avail_in) abort();
      out.resize(result_used);
      lzma_end(&strm);
      return;
    }
    if (ret != LZMA_OK) abort();
    if (strm.avail_out == 0) {
      result_used += avail - strm.avail_out;
      out.resize(out.size() << 1);
      strm.next_out = reinterpret_cast<uint8_t*>(&out[0] + result_used);
      strm.avail_out = avail = out.size() - result_used;
    }
  }
}

std::string CompressWithLzma(const std::string& in, int level) {
  std::string result;
  result.resize(in.size() + (in.size() >> 2) + 128);
  size_t out_pos = 0;
  if (LZMA_OK != lzma_easy_buffer_encode(
                     level, LZMA_CHECK_CRC32, NULL,
                     reinterpret_cast<uint8_t*>(const_cast<char*>(in.data())),
                     in.size(), reinterpret_cast<uint8_t*>(&result[0]),
                     &out_pos, result.size()))
    abort();
  result.resize(out_pos);
  return result;
}

std::string DecompressWithLzma(const std::string& in) {
  static const size_t kMemLimit = 1 << 30;  // 1 GB.
  lzma_stream strm = LZMA_STREAM_INIT;
  std::string result;
  result.resize(8192);
  size_t result_used = 0;
  lzma_ret ret;
  ret = lzma_stream_decoder(&strm, kMemLimit, LZMA_CONCATENATED);
  if (ret != LZMA_OK) abort();
  size_t avail = result.size();
  strm.next_in = reinterpret_cast<const uint8_t*>(in.data());
  strm.avail_in = in.size();
  strm.next_out = reinterpret_cast<uint8_t*>(&result[0]);
  strm.avail_out = avail;
  while (true) {
    ret = lzma_code(&strm, strm.avail_in == 0 ? LZMA_FINISH : LZMA_RUN);
    if (ret == LZMA_STREAM_END) {
      result_used += avail - strm.avail_out;
      if (strm.avail_in) abort();
      result.resize(result_used);
      lzma_end(&strm);
      return result;
    }
    if (ret != LZMA_OK) abort();
    if (strm.avail_out == 0) {
      result_used += avail - strm.avail_out;
      result.resize(result.size() << 1);
      strm.next_out = reinterpret_cast<uint8_t*>(&result[0] + result_used);
      strm.avail_out = avail = result.size() - result_used;
    }
  }
}
