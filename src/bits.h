
#pragma once

#include "bitasm.h"

#define BITS_PER_WORD 64

// Return zero if bit[i] is zero, else return one-bit word with bit[i] set.
static inline uint64 test_bit(uint64 a, uint64 i) { return (a & (1UL << i)); }

// Return a with bit[i] set.
static inline uint64 set_bit(uint64 a, uint64 i) { return (a | (1UL << i)); }

// Return a with bit[i] cleared.
static inline uint64 clear_bit(uint64 a, uint64 i) { return (a & ~(1UL << i)); }

// Return a with bit[i] changed.
static inline uint64 change_bit(uint64 a, uint64 i) { return (a ^ (1UL << i)); }

// Copy bit at [isrc] to position [idst]. Return the modified word.
static inline uint64 copy_bit(uint64 a, uint64 isrc, uint64 idst) {
  uint64 x = ((a >> isrc) ^ (a >> idst)) & 1;
  a ^= (x << idst);
  return a;
}

// Return word with length-n bit block starting at bit p set.
static inline uint64 bit_block(uint64 p, uint64 n) {
  uint64 x = (1UL << n) - 1;
  return x << p;
}

// Return number of set bits.
static inline uint64 bit_count(uint64 x) {
#if defined BITS_USE_ASM
  return asm_bit_count(x);
#else
  x -= (x >> 1) & 0x5555555555555555UL;
  x = ((x >> 2) & 0x3333333333333333UL) + (x & 0x3333333333333333UL);
  x = ((x >> 4) + x) & 0x0f0f0f0f0f0f0f0fUL;
  x *= 0x0101010101010101UL;
  return x >> 56;
#endif
}

// Return number of set bits, must have at most 3 set bits.
static inline uint64 bit_count_3(uint64 x) {
  x -= (x >> 1) & 0x5555555555555555UL;
  x *= 0x5555555555555555UL;
  return x >> 62;
}

// Return number of set bits, must have at most 15 set bits.
static inline uint64 bit_count_15(uint64 x) {
  x -= (x >> 1) & 0x5555555555555555UL;
  x = ((x >> 2) & 0x3333333333333333UL) + (x & 0x3333333333333333UL);
  x *= 0x1111111111111111UL;
  return x >> 60;
}

// Return number of bits set.
static inline uint64 bit_count_sparse(uint64 x) {
  uint64 n = 0;
  do {
    n += (x != 0);
    x &= (x - 1);
    n += (x != 0);
    x &= (x - 1);
    n += (x != 0);
    x &= (x - 1);
    n += (x != 0);
    x &= (x - 1);
  } while (x);
  return n;
}

// Return number of bits set.
static inline uint64 bit_count_dense(uint64 x) {
  return BITS_PER_WORD - bit_count_sparse(~x);
}

// Return number of bits in a special word form 00...0001...11
static inline uint64 bit_count_01(uint64 x) {
#if defined BITS_USE_ASM
  if (1 >= x) return x;
  x = asm_bsr(x);
  return x + 1;
#else
  uint64 ct = 0;
  uint64 a;

  a = (x & (1UL << 32)) >> (32 - 5);
  x >>= a;
  ct += a;

  a = (x & (1UL << 16)) >> (16 - 4);
  x >>= a;
  ct += a;

  a = (x & (1UL << 8)) >> (8 - 3);
  x >>= a;
  ct += a;

  a = (x & (1UL << 4)) >> (4 - 2);
  x >>= a;
  ct += a;

  a = (x & (1UL << 2)) >> (2 - 1);
  x >>= a;
  ct += a;

  a = (x & (1UL << 1)) >> (1 - 0);
  x >>= a;
  ct += a;

  ct += x & 1;
  return ct;
#endif
}

// Return a with bits at positions [k1] and [k2] swapped.
// k1==k2 is allowed (a is unchanged then)
static inline uint64 bit_swap(uint64 a, uint64 k1, uint64 k2) {
  uint64 x = ((a >> k1) ^ (a >> k2)) & 1;
  a ^= (x << k2);
  a ^= (x << k1);
  return a;
}

// Return a with bits at positions [k1] and [k2] swapped.
// Bits must have different values (!)
// (i.e. one is zero, the other one)
// k1==k2 is allowed (a is unchanged then)
static inline uint64 bit_swap_01(uint64 a, uint64 k1, uint64 k2) {
  return a ^ ((1UL << k1) ^ (1UL << k2));
}

// Return x with neighbor bits swapped.
static inline uint64 bit_swap_1(uint64 x) {
  uint64 m = 0x5555555555555555UL;
  return ((x & m) << 1) | ((x & (~m)) >> 1);
}

// Return x with groups of 2 bits swapped.
static inline uint64 bit_swap_2(uint64 x) {
  uint64 m = 0x3333333333333333UL;
  return ((x & m) << 2) | ((x & (~m)) >> 2);
}

// Return x with groups of 4 bits swapped.
static inline uint64 bit_swap_4(uint64 x) {
  uint64 m = 0x0f0f0f0f0f0f0f0fUL;
  return ((x & m) << 4) | ((x & (~m)) >> 4);
}

// Return x with groups of 8 bits swapped.
static inline uint64 bit_swap_8(uint64 x) {
  uint64 m = 0x00ff00ff00ff00ffUL;
  return ((x & m) << 8) | ((x & (~m)) >> 8);
}

// Return x with groups of 16 bits swapped.
static inline uint64 bit_swap_16(uint64 x) {
  uint64 m = 0x0000ffff0000ffffUL;
  return ((x & m) << 16) | ((x & (~m)) >> 16);
}

// Return x with groups of 32 bits swapped.
static inline uint64 bit_swap_32(uint64 x) { return (x << 32) | (x >> 32); }

// Return whether floor(log2(x))==floor(log2(y)).
static inline bool ld_eq(uint64 x, uint64 y) { return ((x ^ y) <= (x & y)); }

// Return whether floor(log2(x))!=floor(log2(y))
static inline bool ld_neq(uint64 x, uint64 y) { return ((x ^ y) > (x & y)); }

// Return index of lowest bit set.
// Bit index ranges from zero to BITS_PER_WORD-1.
// Examples:
//    ***1 --> 0
//    **10 --> 1
//    *100 --> 2
// Return 0 (also) if no bit is set.
static inline uint64 lowest_one_idx(uint64 x) {
#if defined BITS_USE_ASM
  return asm_bsf(x);
#else
  uint64 r = 0;
  x &= -x;
  r |= ((x & 0xffffffff00000000UL) != 0);
  r <<= 1;
  r |= ((x & 0xffff0000ffff0000UL) != 0);
  r <<= 1;
  r |= ((x & 0xff00ff00ff00ff00UL) != 0);
  r <<= 1;
  r |= ((x & 0xf0f0f0f0f0f0f0f0UL) != 0);
  r <<= 1;
  r |= ((x & 0xccccccccccccccccUL) != 0);
  r <<= 1;
  r |= ((x & 0xaaaaaaaaaaaaaaaaUL) != 0);
  return r;
#endif
}

// Return word where only the lowest set bit in x is set.
// Return 0 if no bit is set.
static inline uint64 lowest_one(uint64 x) { return x & -x; }

// Return word where only the lowest unset bit in x is set.
// Return 0 if all bits are set.
static inline uint64 lowest_zero(uint64 x) {
  x = ~x;
  return x & -x;
}

// Return word where the lowest bit set in x is cleared.
// Return 0 for input == 0.
static inline uint64 clear_lowest_one(uint64 x) { return x & (x - 1); }

// Return word where the lowest unset bit in x is set.
// Return ~0 for input == ~0.
static inline uint64 set_lowest_zero(uint64 x) { return x | (x + 1); }

// Return word where only the highest bit in x is set.
// Return 0 if no bit is set.
static inline uint64 highest_one(uint64 x) {
#if defined BITS_USE_ASM
  if (0 == x) return 0;
  x = asm_bsr(x);
  return 1UL << x;
#else
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return x ^ (x >> 1);
#endif
}

// Return word where only the highest unset bit in x is set.
// Return 0 if all bits are set.
static inline uint64 highest_zero(uint64 x) { return highest_one(~x); }

// Return word where the highest unset bit in x is set.
// Return ~0 for input == ~0.
static inline uint64 set_highest_zero(uint64 x) { return x | highest_one(~x); }

// Return index of highest bit set.
// Return 0 if no bit is set.
static inline uint64 highest_one_idx(uint64 x) {
#if defined BITS_USE_ASM
  return asm_bsr(x);
#else
  uint64 MU0 = 0x5555555555555555ULL;
  uint64 MU1 = 0x3333333333333333ULL;
  uint64 MU2 = 0x0f0f0f0f0f0f0f0fULL;
  uint64 MU3 = 0x00ff00ff00ff00ffULL;
  uint64 MU4 = 0x0000ffff0000ffffULL;
  uint64 MU5 = 0x00000000ffffffffULL;

  uint64 r =
      (uint64)ld_neq(x, x & MU0) + ((uint64)ld_neq(x, x & MU1) << 1) +
      ((uint64)ld_neq(x, x & MU2) << 2) + ((uint64)ld_neq(x, x & MU3) << 3) +
      ((uint64)ld_neq(x, x & MU4) << 4) + ((uint64)ld_neq(x, x & MU5) << 5);
  return r;
#endif
}

static inline uint64 highest_zero_idx(uint64 x) { return highest_one_idx(~x); }

// Return floor(log2(x)),
// i.e. return k so that 2^k <= x < 2^(k+1)
// If x==0, then 0 is returned (!)
static inline uint64 ld(uint64 x) {
#if defined BITS_USE_ASM
  return asm_bsr(x);
#else
  return highest_one_idx(x);
#endif
}

// Return whether x == 0(!) or x == 2**k
static inline bool is_pow_of_2(uint64 x) { return !(x & (x - 1)); }

// Return whether x \in {1,2,4,8,16,...}
static inline bool one_bit_q(uint64 x) {
  uint64 m = x - 1;
  return (((x ^ m) >> 1) == m);
}

// Return x if x=2**k
// else return 2**ceil(log_2(x))
// Exception: returns 0 for x==0
static inline uint64 next_pow_of_2(uint64 x) {
  if (is_pow_of_2(x)) return x;
#if defined BITS_USE_ASM
  uint64 n = 1UL << ld(x);
  return n << 1;
#else
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return x + 1;
#endif
}

// Return k if x=2**k else return k+1.
// Exception: returns 0 for x==0.
static inline uint64 next_exp_of_2(uint64 x) {
  if (x <= 1) return 0;
  return ld(x - 1) + 1;
}

// Return index of the i-th set bit of x
// where 0 <= i < bit_count(x).
static inline uint64 ith_one_idx(uint64 x, uint64 i) {
  uint64 x2 = x - ((x >> 1) & 0x5555555555555555ULL);
  uint64 x4 =
      ((x2 >> 2) & 0x3333333333333333ULL) + (x2 & 0x3333333333333333ULL);
  uint64 x8 = ((x4 >> 4) + x4) & 0x0f0f0f0f0f0f0f0fULL;
  uint64 ct = (x8 * 0x0101010101010101ULL) >> 56;

  ++i;
  if (ct < i) return ~0ULL;

  uint64 x16 =
      (0x00ff00ff00ff00ffULL & x8) + (0x00ff00ff00ff00ffULL & (x8 >> 8));
  uint64 x32 =
      (0x0000ffff0000ffffULL & x16) + (0x0000ffff0000ffffULL & (x16 >> 16));

  uint64 w, s = 0;
  w = x32 & 0xffffffffULL;
  if (w < i) {
    s += 32;
    i -= w;
  }

  x16 >>= s;
  w = x16 & 0xffff;
  if (w < i) {
    s += 16;
    i -= w;
  }

  x8 >>= s;
  w = x8 & 0xff;
  if (w < i) {
    s += 8;
    i -= w;
  }

  x4 >>= s;
  w = x4 & 0xf;
  if (w < i) {
    s += 4;
    i -= w;
  }

  x2 >>= s;
  w = x2 & 3;
  if (w < i) {
    s += 2;
    i -= w;
  }

  x >>= s;
  s += ((x & 1) != i);
  return s;
}
