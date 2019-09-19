
#pragma once

#include "types.h"

#if defined __GNUC__

#ifndef DISABLE_BIT_ASM
#ifdef __i386__
#define BITS_USE_ASM
#define BITS_USE_ASM_I386

// Bit Scan Forward: return index of lowest one.
static inline uint64 asm_bsf(uint64 x) {
  asm("bsfl %0, %0" : "=r"(x) : "0"(x));
  return x;
}

// Bit Scan Reverse: return index of highest one.
static inline uint64 asm_bsr(uint64 x) {
  asm("bsrl %0, %0" : "=r"(x) : "0"(x));
  return x;
}

// Byte swap
static inline uint64 asm_bswap(uint64 x) {
  asm("bswap %0" : "=r"(x) : "0"(x));
  return x;
}

// Rotate Left
static inline uint64 asm_rol(uint64 x, uint64 r) {
  asm("roll   %%cl, %0" : "=r"(x) : "0"(x), "c"(r));
  return x;
}

// Rotate Right
static inline uint64 asm_ror(uint64 x, uint64 r) {
  asm("rorl   %%cl, %0" : "=r"(x) : "0"(x), "c"(r));
  return x;
}

// Return the parity of x (which is the
// _complement_ of x86's parity flag).
// As parity is for the low byte only,
// therefore we need to prepend
// x ^= (x>>16);  x ^= (x>>8);
// in the code
static inline uint64 asm_parity(uint64 x) {
  x ^= (x >> 16);
  x ^= (x >> 8);
  asm("addl  $0, %0  \n"
      "setnp %%al    \n"
      "movzx %%al, %0"
      : "=r"(x)
      : "0"(x)
      : "eax");

  return x;
}

// Bit Test
static inline uint64 asm_bt(const uint64 *f, uint64 i) {
  uint64 ret;
  asm("btl  %2, %1 \n"  // carry = 0 or 1
      "sbbl %0, %0"     // ret = 0 or -1
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Test and Set
static inline uint64 asm_bts(uint64 *f, uint64 i) {
  uint64 ret;
  asm("btsl %2, %1 \n"
      "sbbl %0, %0"
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Set
static inline void asm_b_s(uint64 *f, uint64 i) {
  asm("btsl  %1, %0 \n"
      : /* void */
      : "m"(*f), "r"(i));
}

// Bit Test and Reset
static inline uint64 asm_btr(uint64 *f, uint64 i) {
  uint64 ret;
  asm("btrl  %2, %1 \n"
      "sbbl %0, %0"
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Reset
static inline void asm_b_r(uint64 *f, uint64 i) {
  asm("btrl  %1, %0 \n"
      : /* void */
      : "m"(*f), "r"(i));
}

// Bit Test and Complement
static inline uint64 asm_btc(uint64 *f, uint64 i) {
  uint64 ret;
  asm("btcl  %2, %1 \n"
      "sbbl %0, %0"
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Complement
static inline void asm_b_c(uint64 *f, uint64 i) {
  asm("btcl  %1, %0 \n"
      : /* void */
      : "m"(*f), "r"(i));
}

#endif

#ifdef __x86_64__
#define BITS_USE_ASM        // use asm code
#define BITS_USE_ASM_AMD64  // use AMD64 asm code

static inline uint64 asm_bit_count(uint64 x) {
  asm("popcntq %0, %0" : "=r"(x) : "0"(x));
  return x;
}

// AMD 24592:
// The BSF and BSR instructions search a source operand for the
// least-significant (BSF) or most-significant (BSR) bit that is set
// to 1. If a set bit is found, its bit index is loaded into the
// destination operand, and the zero flag (ZF) is set. If no set bit
// is found , the zero flag is cleared and the contents of the
// destination are undefined.

// Bit Scan Forward
static inline uint64 asm_bsf(uint64 x) {
  asm("bsfq %0, %0" : "=r"(x) : "0"(x));
  return x;
}

// Bit Scan Reverse
static inline uint64 asm_bsr(uint64 x) {
  asm("bsrq %0, %0" : "=r"(x) : "0"(x));
  return x;
}

// Byte swap
static inline uint64 asm_bswap(uint64 x) {
  asm("bswap %0" : "=r"(x) : "0"(x));
  return x;
}

// Rotate Left
static inline uint64 asm_rol(uint64 x, uint64 r) {
  asm("rolq   %%cl, %0" : "=r"(x) : "0"(x), "c"(r));
  return x;
}

// Rotate Right
static inline uint64 asm_ror(uint64 x, uint64 r) {
  asm("rorq   %%cl, %0" : "=r"(x) : "0"(x), "c"(r));
  return x;
}

// Return the parity of x (which is the
// _complement_ of AMD64's parity flag).
// As parity is for the low byte only,
// therefore we need to prepend
// x ^= (x>>32);  x ^= (x>>16);  x ^= (x>>8);
// in the code
static inline uint64 asm_parity(uint64 x) {
  x ^= (x >> 32);
  x ^= (x >> 16);
  x ^= (x >> 8);
  asm("addq  $0, %0  \n"
      "setnp %%al    \n"
      "movzx %%al, %0"
      : "=r"(x)
      : "0"(x)
      : "rax");

  return x;
}

// AMD 24592:
// The BTx instructions copy a specified bit in the first operand to
// the carry flag (CF) and leave the source bit unchanged (BT), or
// complement the source bit (BTC), or clear the source bit to 0
// (BTR), or set the source bit to 1 (BTS).

// Bit Test
static inline uint64 asm_bt(const uint64 *f, uint64 i) {
  uint64 ret;
  asm("btq  %2, %1 \n"  // carry = 0 or 1
      "sbbq %0, %0"     // ret = 0 or -1
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Test and Set
static inline uint64 asm_bts(uint64 *f, uint64 i) {
  uint64 ret;
  asm("btsq %2, %1 \n"
      "sbbq %0, %0"
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Set
static inline void asm_b_s(uint64 *f, uint64 i) {
  asm("btsq  %1, %0 \n"
      : /* void */
      : "m"(*f), "r"(i));
}

// Bit Test and Reset
static inline uint64 asm_btr(uint64 *f, uint64 i) {
  uint64 ret;
  asm("btrq  %2, %1 \n"
      "sbbq %0, %0"
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}

// Bit Reset
static inline void asm_b_r(uint64 *f, uint64 i) {
  asm("btrq  %1, %0 \n"
      : /* void */
      : "m"(*f), "r"(i));
}

// Bit Test and Complement
static inline uint64 asm_btc(uint64 *f, uint64 i) {
  uint64 ret;
  asm("btcq  %2, %1 \n"
      "sbbq %0, %0"
      : "=r"(ret)
      : "m"(*f), "r"(i));
  return ret;
}
// -------------------------

// Bit Complement
static inline void asm_b_c(uint64 *f, uint64 i) {
  asm("btcq  %1, %0 \n"
      : /* void */
      : "m"(*f), "r"(i));
}

#endif
#endif

#endif
