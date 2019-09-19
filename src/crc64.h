#pragma once

#include "types.h"

class crc64 {
 public:
  static const uint64 cc[];
  
  explicit crc64(uint64 c = 0) {
    if (0 == c) c = 0x1bULL;
    init(c);
  }

  ~crc64() {}

  uint64 get_a() const { return a_; }

  uint64 byte_in(uchar b) {
    a_ ^= b;
    shift();
    shift();
    shift();
    shift();
    shift();
    shift();
    shift();
    shift();
    return a_;
  }

 private:
  void reset() { set_a(~0ULL); }
  void set_a(uint64 a) { a_ = a; }

  void init(uint64 c) {
    c_ = c;
    c_ >>= 1;
    uint64 h = 1ULL << 63;
    c_ |= h;
    reset();
  }

  void shift() {
    bool s = (a_ & 1);
    a_ >>= 1;
    if (0 != s) a_ ^= c_;
  }

  uint64 a_;
  uint64 c_;
};
