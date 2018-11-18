/**********************************************************************
 * Copyright (c) 2013, 2014 Pieter Wuille                             *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_FIELD_INNER5X52_IMPL_H_
#define _SECP256K1_FIELD_INNER5X52_IMPL_H_

#include <stdint.h>

#ifdef VERIFY
#define VERIFY_BITS(x, n) VERIFY_CHECK(((x) >> (n)) == 0)
#else
#define VERIFY_BITS(x, n) do { } while(0)
#endif

SECP256K1_INLINE void secp256k1_fe_mul_inner(uint64_t *r, const uint64_t *a, const uint64_t * SECP256K1_RESTRICT b);
SECP256K1_INLINE void secp256k1_fe_sqr_inner(uint64_t *r, const uint64_t *a);

#endif
