#pragma once
#include <stdint.h>
#include <stddef.h>

#ifdef _WIN32
#  define ULAM_API __declspec(dllexport)
#else
#  define ULAM_API __attribute__((visibility("default")))
#endif

typedef struct UlamState UlamState;

ULAM_API UlamState *ulam_init(int track_pairs);
ULAM_API void       ulam_free(UlamState *s);
ULAM_API int        ulam_compute_up_to(UlamState *s, uint64_t n);

ULAM_API int        ulam_is_ulam(UlamState *s, uint64_t n);

ULAM_API int        ulam_save(UlamState *s, const char *path);
ULAM_API UlamState *ulam_load(const char *path);

ULAM_API uint64_t   ulam_max_computed(UlamState *s);
ULAM_API uint64_t   ulam_count(UlamState *s);
ULAM_API uint64_t   ulam_get(UlamState *s, uint64_t index);
ULAM_API uint64_t  *ulam_get_ptr(UlamState *s);
ULAM_API uint64_t  *ulam_find_adjacent_pairs(UlamState *s, uint64_t *out_count);

ULAM_API void       ulam_free_ptr(void *p);
ULAM_API int        ulam_open_adj_stream(UlamState *s, const char *path);
ULAM_API void       ulam_close_adj_stream(UlamState *s);

ULAM_API int        ulam_track_pairs(UlamState *s);
ULAM_API uint64_t   ulam_pairs_count(UlamState *s);
ULAM_API void       ulam_get_pair(UlamState *s, uint64_t idx, uint64_t *a, uint64_t *b);
ULAM_API uint64_t  *ulam_pairs_ptr(UlamState *s);
