#pragma once

#define DEBUG 1
#if DEBUG
#define do_on_debug(code) { code };
#else
#define do_on_debug(code) { };
#endif

#define INSERT_LEX_LT_CODE(left_field, right_field) \
    if (left_field < right_field) return true; \
    if (left_field > right_field) return false;

