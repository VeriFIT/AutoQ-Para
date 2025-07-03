#pragma once

#include <cstdint>

#define DEBUG 1
#if DEBUG
#define do_on_debug(code) { code };
#else
#define do_on_debug(code) { };
#endif

#define INSERT_LEX_LT_CODE(left_field, right_field) \
    if (left_field < right_field) return true; \
    if (left_field > right_field) return false;


typedef int8_t   s8;
typedef int64_t  s16;
typedef int64_t  s32;
typedef int64_t  s64;

typedef uint8_t  u8;
typedef uint64_t u16;
typedef uint64_t u32;
typedef uint64_t u64;

