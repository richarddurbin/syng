#ifndef CSYNCMER_FAST_H
#define CSYNCMER_FAST_H

// Fast closed syncmer detection using ntHash32 - Pure C implementation
// Includes AVX2 SIMD acceleration when available

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
#include <atomic>
#define CSYNCMER_ATOMIC_INT std::atomic<int>
#define CSYNCMER_ATOMIC_LOAD(x) (x).load()
#define CSYNCMER_ATOMIC_STORE(x, v) (x).store(v)
#define CSYNCMER_THREAD_LOCAL thread_local
#define CSYNCMER_ALIGNAS(n) alignas(n)
#else
#include <stdatomic.h>
#include <stdalign.h>
#define CSYNCMER_ATOMIC_INT atomic_int
#define CSYNCMER_ATOMIC_LOAD(x) atomic_load(&(x))
#define CSYNCMER_ATOMIC_STORE(x, v) atomic_store(&(x), v)
#define CSYNCMER_THREAD_LOCAL _Thread_local
#define CSYNCMER_ALIGNAS(n) alignas(n)
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif

#ifdef __BMI2__
#include <x86intrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// ntHash32 Constants and Helpers
// ============================================================================

static const uint32_t CSYNCMER_NTHASH32_F[4] = {
    0x95c60474,  // A
    0x62a02b4c,  // C
    0x82572324,  // T
    0x4be24456   // G
};

static inline uint32_t csyncmer_rotl7_32(uint32_t x) {
    return (x << 7) | (x >> 25);
}

static inline void csyncmer_make_f_rot_32(size_t S, uint32_t f_rot[4]) {
    // ntHash uses (k-1)*R rotation, matching nthash.hpp::make_f_rot
    uint32_t rot = (uint32_t)(((S - 1) * 7) & 31);
    for (int i = 0; i < 4; i++) {
        uint32_t x = CSYNCMER_NTHASH32_F[i];
        f_rot[i] = rot ? ((x << rot) | (x >> (32 - rot))) : x;
    }
}

// ============================================================================
// ASCII Lookup Tables
// ============================================================================

static inline void csyncmer_init_ascii_hash_table(uint32_t table[256]) {
    memset(table, 0, 256 * sizeof(uint32_t));
    table['A'] = table['a'] = CSYNCMER_NTHASH32_F[0];
    table['C'] = table['c'] = CSYNCMER_NTHASH32_F[1];
    table['T'] = table['t'] = CSYNCMER_NTHASH32_F[2];
    table['G'] = table['g'] = CSYNCMER_NTHASH32_F[3];
}

static inline void csyncmer_init_ascii_to_idx(uint8_t table[256]) {
    memset(table, 0, 256);
    table['A'] = table['a'] = 0;
    table['C'] = table['c'] = 1;
    table['T'] = table['t'] = 2;
    table['G'] = table['g'] = 3;
}

// ============================================================================
// Scalar Implementation: Fused RESCAN with Branch-free Updates
// Macro-generated for forward (CANONICAL=0) and canonical (CANONICAL=1) variants.
// Count vs positions handled via runtime if(out_positions) — compiler eliminates
// the branch when inlined with NULL.
// ============================================================================

static inline uint32_t csyncmer_rotr7_32(uint32_t x) {
    return (x >> 7) | (x << 25);
}

#define CSYNCMER_DEFINE_RESCAN_32(FUNC_NAME, CANONICAL)                        \
static inline size_t FUNC_NAME(                                                \
    const char* sequence,                                                      \
    size_t length,                                                             \
    size_t K,                                                                  \
    size_t S,                                                                  \
    uint32_t* out_positions,                                                   \
    uint8_t* out_strands,                                                      \
    size_t max_positions                                                        \
) {                                                                            \
    if (!sequence || length < K || S == 0 || S >= K) {                         \
        return 0;                                                              \
    }                                                                          \
    if (out_positions && max_positions == 0) return 0;                          \
                                                                               \
    /* Initialize lookup tables */                                             \
    uint32_t F_ASCII[256];                                                     \
    uint8_t IDX_ASCII[256];                                                    \
    uint32_t f_rot[4];                                                         \
    csyncmer_init_ascii_hash_table(F_ASCII);                                   \
    csyncmer_init_ascii_to_idx(IDX_ASCII);                                     \
    csyncmer_make_f_rot_32(S, f_rot);                                          \
                                                                               \
    /* Canonical-only: RC rotation tables */                                   \
    uint32_t c_rot[4];                                                         \
    uint32_t RC32[4];                                                          \
    uint32_t RC32_ROTR7[4];                                                    \
    if (CANONICAL) {                                                           \
        uint32_t rot = (uint32_t)(((S - 1) * 7) & 31);                        \
        RC32[0] = 0x82572324; RC32[1] = 0x4be24456;                            \
        RC32[2] = 0x95c60474; RC32[3] = 0x62a02b4c;                            \
        RC32_ROTR7[0] = 0x4904ae46; RC32_ROTR7[1] = 0xac97c488;               \
        RC32_ROTR7[2] = 0xe92b8c08; RC32_ROTR7[3] = 0x98c54056;              \
        for (int ci = 0; ci < 4; ci++) {                                       \
            c_rot[ci] = rot ? ((RC32[ci] << rot) | (RC32[ci] >> (32 - rot))) : RC32[ci]; \
        }                                                                      \
    } else {                                                                   \
        (void)c_rot; (void)RC32; (void)RC32_ROTR7;                            \
    }                                                                          \
                                                                               \
    size_t window_size = K - S + 1;                                            \
    size_t num_smers = length - S + 1;                                         \
    const uint8_t* seq = (const uint8_t*)sequence;                             \
                                                                               \
    /* Power-of-2 circular buffer */                                           \
    size_t buf_size = 1;                                                       \
    while (buf_size < window_size) buf_size <<= 1;                             \
    size_t buf_mask = buf_size - 1;                                            \
    uint32_t* hash_buffer = (uint32_t*)malloc(buf_size * sizeof(uint32_t));     \
    if (!hash_buffer) return 0;                                                \
                                                                               \
    size_t syncmer_count = 0;                                                  \
                                                                               \
    /* First s-mer hash (forward) */                                           \
    uint32_t fw_hash = 0;                                                      \
    for (size_t j = 0; j < S; ++j) {                                           \
        fw_hash = csyncmer_rotl7_32(fw_hash) ^ F_ASCII[seq[j]];               \
    }                                                                          \
                                                                               \
    /* Canonical: also compute RC hash for first s-mer */                      \
    uint32_t rc_hash = 0;                                                      \
    if (CANONICAL) {                                                           \
        for (size_t j = 0; j < S; ++j) {                                       \
            rc_hash = csyncmer_rotl7_32(rc_hash) ^ RC32[IDX_ASCII[seq[S - 1 - j]]]; \
        }                                                                      \
    }                                                                          \
                                                                               \
    uint32_t effective_hash = fw_hash;                                         \
    if (CANONICAL) {                                                           \
        effective_hash = (fw_hash <= rc_hash) ? fw_hash : rc_hash;             \
    }                                                                          \
    hash_buffer[0] = effective_hash;                                           \
                                                                               \
    /* Rolling state for forward */                                            \
    uint32_t fw = fw_hash ^ f_rot[IDX_ASCII[seq[0]]];                         \
                                                                               \
    /* Fill initial window */                                                  \
    for (size_t i = 1; i < window_size; ++i) {                                 \
        uint8_t new_last = seq[i + S - 1];                                     \
        uint8_t new_first = seq[i];                                            \
                                                                               \
        /* Forward rolling */                                                  \
        fw_hash = csyncmer_rotl7_32(fw) ^ F_ASCII[new_last];                   \
        fw = fw_hash ^ f_rot[IDX_ASCII[new_first]];                           \
                                                                               \
        if (CANONICAL) {                                                       \
            /* RC rolling: rc = rotr7(rc) ^ rotr7(RC[old_first]) ^ c_rot[new_last] */ \
            uint8_t old_first = seq[i - 1];                                    \
            uint8_t idx_old = IDX_ASCII[old_first];                            \
            uint8_t idx_new = IDX_ASCII[new_last];                             \
            rc_hash = csyncmer_rotr7_32(rc_hash) ^ RC32_ROTR7[idx_old] ^ c_rot[idx_new]; \
            effective_hash = (fw_hash <= rc_hash) ? fw_hash : rc_hash;         \
        } else {                                                               \
            effective_hash = fw_hash;                                          \
        }                                                                      \
        hash_buffer[i & buf_mask] = effective_hash;                            \
    }                                                                          \
                                                                               \
    /* Find minimum in first window */                                         \
    uint32_t min_hash = UINT32_MAX;                                            \
    size_t min_pos = 0;                                                        \
    for (size_t i = 0; i < window_size; ++i) {                                 \
        if (hash_buffer[i] < min_hash) {                                       \
            min_hash = hash_buffer[i];                                         \
            min_pos = i;                                                       \
        }                                                                      \
    }                                                                          \
                                                                               \
    /* Check first k-mer */                                                    \
    if (min_pos == 0 || min_pos == window_size - 1) {                          \
        if (out_positions) {                                                   \
            out_positions[syncmer_count] = 0;                                  \
            if (CANONICAL && out_strands) {                                    \
                out_strands[syncmer_count] = 0; /* strand not tracked in rescan */ \
            }                                                                  \
        }                                                                      \
        syncmer_count++;                                                       \
        if (out_positions && syncmer_count >= max_positions) {                  \
            free(hash_buffer);                                                 \
            return syncmer_count;                                              \
        }                                                                      \
    }                                                                          \
                                                                               \
    /* Main loop with branch-free min update */                                \
    for (size_t kmer_idx = 1; kmer_idx < num_smers - window_size + 1; ++kmer_idx) { \
        size_t i = kmer_idx + window_size - 1;                                 \
        uint8_t new_last = seq[i + S - 1];                                     \
        uint8_t new_first = seq[i];                                            \
                                                                               \
        /* Forward rolling */                                                  \
        fw_hash = csyncmer_rotl7_32(fw) ^ F_ASCII[new_last];                   \
        fw = fw_hash ^ f_rot[IDX_ASCII[new_first]];                           \
                                                                               \
        if (CANONICAL) {                                                       \
            uint8_t old_first = seq[i - 1];                                    \
            uint8_t idx_old = IDX_ASCII[old_first];                            \
            uint8_t idx_new = IDX_ASCII[new_last];                             \
            rc_hash = csyncmer_rotr7_32(rc_hash) ^ RC32_ROTR7[idx_old] ^ c_rot[idx_new]; \
            effective_hash = (fw_hash <= rc_hash) ? fw_hash : rc_hash;         \
        } else {                                                               \
            effective_hash = fw_hash;                                          \
        }                                                                      \
        hash_buffer[i & buf_mask] = effective_hash;                            \
                                                                               \
        /* RESCAN: update minimum */                                           \
        if (min_pos < kmer_idx) {                                              \
            min_hash = UINT32_MAX;                                             \
            for (size_t j = kmer_idx; j <= i; ++j) {                           \
                uint32_t h = hash_buffer[j & buf_mask];                        \
                if (h < min_hash) {                                            \
                    min_hash = h;                                              \
                    min_pos = j;                                               \
                }                                                              \
            }                                                                  \
        } else {                                                               \
            if (effective_hash < min_hash) {                                   \
                min_hash = effective_hash;                                      \
                min_pos = i;                                                   \
            }                                                                  \
        }                                                                      \
                                                                               \
        /* Check closed syncmer condition */                                   \
        size_t min_offset = min_pos - kmer_idx;                                \
        if (min_offset == 0 || min_offset == window_size - 1) {                \
            if (out_positions) {                                               \
                out_positions[syncmer_count] = (uint32_t)kmer_idx;             \
                if (CANONICAL && out_strands) {                                \
                    out_strands[syncmer_count] = 0;                            \
                }                                                              \
            }                                                                  \
            syncmer_count++;                                                   \
            if (out_positions && syncmer_count >= max_positions) {              \
                free(hash_buffer);                                             \
                return syncmer_count;                                          \
            }                                                                  \
        }                                                                      \
    }                                                                          \
                                                                               \
    free(hash_buffer);                                                         \
    return syncmer_count;                                                      \
}

/* Generate forward and canonical RESCAN internal functions */
CSYNCMER_DEFINE_RESCAN_32(csyncmer_rescan_32_, 0)
CSYNCMER_DEFINE_RESCAN_32(csyncmer_canonical_rescan_32_, 1)

/* Public API wrappers preserving existing signatures */

static inline size_t csyncmer_rescan_32_count(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    size_t* num_syncmers
) {
    size_t c = csyncmer_rescan_32_(sequence, length, K, S, NULL, NULL, 0);
    if (num_syncmers) *num_syncmers = c;
    return c;
}

static inline size_t csyncmer_rescan_32_positions(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    uint32_t* out_positions,
    size_t max_positions
) {
    return csyncmer_rescan_32_(sequence, length, K, S, out_positions, NULL, max_positions);
}

static inline size_t csyncmer_canonical_rescan_32_count(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    return csyncmer_canonical_rescan_32_(sequence, length, K, S, NULL, NULL, 0);
}

static inline size_t csyncmer_canonical_rescan_32_positions(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    uint32_t* out_positions,
    uint8_t* out_strands,
    size_t max_positions
) {
    return csyncmer_canonical_rescan_32_(sequence, length, K, S, out_positions, out_strands, max_positions);
}

// ============================================================================
// ntHash64 Constants and Helpers
// ============================================================================

static const uint64_t CSYNCMER_NTHASH64_F[4] = {
    0x3c8bfbb395c60474ULL,  // A
    0x3193c18562a02b4cULL,  // C
    0x20323ed082572324ULL,  // T
    0x295549f54be24456ULL   // G
};

// Reverse complement constants: A↔T (0↔2), C↔G (1↔3)
// RC[base] = F[complement(base)]
// Index mapping: A=0, C=1, T=2, G=3
static const uint64_t CSYNCMER_NTHASH64_RC[4] = {
    0x20323ed082572324ULL,  // RC[A] (idx 0) = F[T] = F[2]
    0x295549f54be24456ULL,  // RC[C] (idx 1) = F[G] = F[3]
    0x3c8bfbb395c60474ULL,  // RC[T] (idx 2) = F[A] = F[0]
    0x3193c18562a02b4cULL   // RC[G] (idx 3) = F[C] = F[1]
};

// Pre-computed rotr7(RC[base]) for optimized canonical iterator
// Eliminates one rotation per iteration in RC rolling formula
static const uint64_t CSYNCMER_NTHASH64_RC_ROTR7[4] = {
    0x4840647da104ae46ULL,  // rotr7(RC[A]) = rotr7(0x20323ed082572324)
    0xac52aa93ea97c488ULL,  // rotr7(RC[C]) = rotr7(0x295549f54be24456)
    0xe87917f7672b8c08ULL,  // rotr7(RC[T]) = rotr7(0x3c8bfbb395c60474)
    0x986327830ac54056ULL   // rotr7(RC[G]) = rotr7(0x3193c18562a02b4c)
};

static inline uint64_t csyncmer_rotl7_64(uint64_t x) {
    return (x << 7) | (x >> 57);
}

static inline uint64_t csyncmer_rotr7_64(uint64_t x) {
    return (x >> 7) | (x << 57);
}

static inline void csyncmer_make_f_rot_64(size_t S, uint64_t f_rot[4]) {
    // ntHash uses (k-1)*R rotation
    uint32_t rot = (uint32_t)(((S - 1) * 7) & 63);
    for (int i = 0; i < 4; i++) {
        uint64_t x = CSYNCMER_NTHASH64_F[i];
        f_rot[i] = rot ? ((x << rot) | (x >> (64 - rot))) : x;
    }
}

// RC rotation table for canonical hashing
// Used for: rol(RC[new_last], (S-1)*7) in RC rolling formula
static inline void csyncmer_make_c_rot_64(size_t S, uint64_t c_rot[4]) {
    uint32_t rot = (uint32_t)(((S - 1) * 7) & 63);
    for (int i = 0; i < 4; i++) {
        uint64_t x = CSYNCMER_NTHASH64_RC[i];
        c_rot[i] = rot ? ((x << rot) | (x >> (64 - rot))) : x;  // left rotation
    }
}

static inline void csyncmer_init_ascii_hash_table_64(uint64_t table[256]) {
    memset(table, 0, 256 * sizeof(uint64_t));
    table['A'] = table['a'] = CSYNCMER_NTHASH64_F[0];
    table['C'] = table['c'] = CSYNCMER_NTHASH64_F[1];
    table['T'] = table['t'] = CSYNCMER_NTHASH64_F[2];
    table['G'] = table['g'] = CSYNCMER_NTHASH64_F[3];
}

static inline void csyncmer_init_ascii_rc_table_64(uint64_t table[256]) {
    memset(table, 0, 256 * sizeof(uint64_t));
    table['A'] = table['a'] = CSYNCMER_NTHASH64_RC[0];  // A→T
    table['C'] = table['c'] = CSYNCMER_NTHASH64_RC[1];  // C→G
    table['T'] = table['t'] = CSYNCMER_NTHASH64_RC[2];  // T→A
    table['G'] = table['g'] = CSYNCMER_NTHASH64_RC[3];  // G→C
}

// Pre-computed rotr7(RC[base]) lookup table for canonical iterator optimization
static inline void csyncmer_init_ascii_rc_rotr7_table_64(uint64_t table[256]) {
    memset(table, 0, 256 * sizeof(uint64_t));
    table['A'] = table['a'] = CSYNCMER_NTHASH64_RC_ROTR7[0];
    table['C'] = table['c'] = CSYNCMER_NTHASH64_RC_ROTR7[1];
    table['T'] = table['t'] = CSYNCMER_NTHASH64_RC_ROTR7[2];
    table['G'] = table['g'] = CSYNCMER_NTHASH64_RC_ROTR7[3];
}

// ============================================================================
// Shared Static Lookup Tables for 64-bit Iterator
// ============================================================================

static uint64_t CSYNCMER_F_ASCII_64[256];
static uint64_t CSYNCMER_RC_ASCII_64[256];
static uint64_t CSYNCMER_RC_ROTR7_ASCII_64[256];  // Pre-computed rotr7(RC[base])
static uint8_t CSYNCMER_IDX_ASCII[256];
static CSYNCMER_ATOMIC_INT CSYNCMER_TABLES_INITIALIZED = {0};

static inline void csyncmer_ensure_tables_initialized(void) {
    if (CSYNCMER_ATOMIC_LOAD(CSYNCMER_TABLES_INITIALIZED)) return;
    // Initialization is idempotent, so concurrent writes are harmless
    csyncmer_init_ascii_hash_table_64(CSYNCMER_F_ASCII_64);
    csyncmer_init_ascii_rc_table_64(CSYNCMER_RC_ASCII_64);
    csyncmer_init_ascii_rc_rotr7_table_64(CSYNCMER_RC_ROTR7_ASCII_64);
    csyncmer_init_ascii_to_idx(CSYNCMER_IDX_ASCII);
    CSYNCMER_ATOMIC_STORE(CSYNCMER_TABLES_INITIALIZED, 1);
}

// ============================================================================
// 64-bit Iterator API (scalar, portable, exact)
// ============================================================================

// Optimized iterator state - shared lookup tables, local variable caching
typedef struct CsyncmerIterator64 {
    const uint8_t* seq;
    uint64_t* hash_buffer;
    uint64_t hash;
    uint64_t fw;
    uint64_t min_hash;
    size_t min_pos;
    size_t kmer_idx;
    size_t num_kmers;
    size_t window_size;
    size_t S;
    size_t buf_mask;
    uint64_t f_rot[4];  // Small rotation table (fits in cache line)
} CsyncmerIterator64;

// Create iterator for closed syncmer detection
// K = total syncmer length, S = smer size (window_size = K - S + 1)
static inline CsyncmerIterator64* csyncmer_iterator_create_64(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    if (!sequence || length < K || S == 0 || S >= K) {
        return NULL;
    }

    // Ensure shared lookup tables are initialized
    csyncmer_ensure_tables_initialized();

    CsyncmerIterator64* iter = (CsyncmerIterator64*)malloc(sizeof(CsyncmerIterator64));
    if (!iter) return NULL;

    iter->seq = (const uint8_t*)sequence;
    iter->S = S;
    iter->window_size = K - S + 1;
    iter->num_kmers = length - K + 1;

    // Initialize S-dependent rotation table
    csyncmer_make_f_rot_64(S, iter->f_rot);

    // Allocate power-of-2 circular buffer
    size_t buf_size = 1;
    while (buf_size < iter->window_size) buf_size <<= 1;
    iter->buf_mask = buf_size - 1;
    iter->hash_buffer = (uint64_t*)malloc(buf_size * sizeof(uint64_t));
    if (!iter->hash_buffer) {
        free(iter);
        return NULL;
    }

    // Use shared static tables for initialization
    const uint64_t* F_ASCII = CSYNCMER_F_ASCII_64;
    const uint8_t* IDX_ASCII = CSYNCMER_IDX_ASCII;

    // Compute first s-mer hash
    uint64_t hash = 0;
    for (size_t j = 0; j < S; ++j) {
        hash = csyncmer_rotl7_64(hash) ^ F_ASCII[iter->seq[j]];
    }
    iter->hash = hash;
    iter->fw = hash ^ iter->f_rot[IDX_ASCII[iter->seq[0]]];
    iter->hash_buffer[0] = hash;

    // Fill initial window with remaining s-mer hashes
    for (size_t i = 1; i < iter->window_size; ++i) {
        hash = csyncmer_rotl7_64(iter->fw) ^ F_ASCII[iter->seq[i + S - 1]];
        iter->fw = hash ^ iter->f_rot[IDX_ASCII[iter->seq[i]]];
        iter->hash_buffer[i & iter->buf_mask] = hash;
    }
    iter->hash = hash;

    // Find minimum in first window
    uint64_t min_hash = UINT64_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < iter->window_size; ++i) {
        uint64_t h = iter->hash_buffer[i];
        if (h < min_hash) {
            min_hash = h;
            min_pos = i;
        }
    }
    iter->min_hash = min_hash;
    iter->min_pos = min_pos;
    iter->kmer_idx = 0;

    return iter;
}

// Get next syncmer position. Returns 1 if valid, 0 when exhausted.
static inline int csyncmer_iterator_next_64(
    CsyncmerIterator64* iter,
    size_t* pos
) {
    if (!iter) return 0;

    // Cache hot state in local variables for register allocation
    const uint8_t* seq = iter->seq;
    size_t kmer_idx = iter->kmer_idx;
    const size_t num_kmers = iter->num_kmers;
    const size_t window_size = iter->window_size;
    const size_t S = iter->S;
    uint64_t hash = iter->hash;
    uint64_t fw = iter->fw;
    uint64_t min_hash = iter->min_hash;
    size_t min_pos = iter->min_pos;
    const size_t buf_mask = iter->buf_mask;
    uint64_t* hash_buffer = iter->hash_buffer;
    const uint64_t* f_rot = iter->f_rot;

    // Use shared static tables
    const uint64_t* F_ASCII = CSYNCMER_F_ASCII_64;
    const uint8_t* IDX_ASCII = CSYNCMER_IDX_ASCII;

    // Check first k-mer on first call
    if (kmer_idx == 0) {
        size_t min_offset = min_pos;
        if (min_offset == 0 || min_offset == window_size - 1) {
            *pos = 0;
            iter->kmer_idx = 1;
            return 1;
        }
        kmer_idx = 1;
    }

    // Main loop
    while (kmer_idx < num_kmers) {
        size_t i = kmer_idx + window_size - 1;

        // Compute new hash for the new s-mer entering the window
        hash = csyncmer_rotl7_64(fw) ^ F_ASCII[seq[i + S - 1]];
        fw = hash ^ f_rot[IDX_ASCII[seq[i]]];
        hash_buffer[i & buf_mask] = hash;

        // RESCAN: update minimum
        if (min_pos < kmer_idx) {
            // Minimum fell out of window - rescan
            min_hash = UINT64_MAX;
            for (size_t j = kmer_idx; j <= i; ++j) {
                uint64_t h = hash_buffer[j & buf_mask];
                if (h < min_hash) {
                    min_hash = h;
                    min_pos = j;
                }
            }
        } else {
            // Check if new hash is smaller
            if (hash < min_hash) {
                min_hash = hash;
                min_pos = i;
            }
        }

        // Check closed syncmer condition
        size_t min_offset = min_pos - kmer_idx;
        size_t current_kmer = kmer_idx;
        kmer_idx++;

        if (min_offset == 0 || min_offset == window_size - 1) {
            // Write back modified state
            iter->hash = hash;
            iter->fw = fw;
            iter->min_hash = min_hash;
            iter->min_pos = min_pos;
            iter->kmer_idx = kmer_idx;
            *pos = current_kmer;
            return 1;
        }
    }

    // Write back state on exhaustion
    iter->hash = hash;
    iter->fw = fw;
    iter->min_hash = min_hash;
    iter->min_pos = min_pos;
    iter->kmer_idx = kmer_idx;
    return 0;  // Exhausted
}

// Free iterator
static inline void csyncmer_iterator_destroy_64(CsyncmerIterator64* iter) {
    if (iter) {
        if (iter->hash_buffer) {
            free(iter->hash_buffer);
        }
        free(iter);
    }
}

// ============================================================================
// 64-bit Canonical Iterator API (strand-independent, scalar, portable, exact)
// ============================================================================

// Canonical iterator state - computes min(forward, reverse_complement) hash
typedef struct CsyncmerIteratorCanonical64 {
    const uint8_t* seq;
    uint64_t* hash_buffer;     // Stores min(fw, rc) for each s-mer
    uint8_t* strand_buffer;    // Stores strand info: 0=forward, 1=RC
    uint64_t fw_hash;          // Current forward hash
    uint64_t rc_hash;          // Current reverse complement hash
    uint64_t min_hash;
    size_t min_pos;
    size_t kmer_idx;
    size_t num_kmers;
    size_t window_size;
    size_t S;
    size_t buf_mask;
    uint64_t f_rot[4];         // Forward rotation table: F[base] << ((S-1)*7)
    uint64_t c_rot[4];         // RC rotation table: RC[base] << ((S-1)*7)
} CsyncmerIteratorCanonical64;

// Create canonical iterator for strand-independent syncmer detection
// K = total syncmer length, S = smer size (window_size = K - S + 1)
// Returns NULL on invalid input or allocation failure
static inline CsyncmerIteratorCanonical64* csyncmer_iterator_create_canonical_64(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    if (!sequence || length < K || S == 0 || S >= K) {
        return NULL;
    }

    // Ensure shared lookup tables are initialized
    csyncmer_ensure_tables_initialized();

    CsyncmerIteratorCanonical64* iter = (CsyncmerIteratorCanonical64*)malloc(sizeof(CsyncmerIteratorCanonical64));
    if (!iter) return NULL;

    iter->seq = (const uint8_t*)sequence;
    iter->S = S;
    iter->window_size = K - S + 1;
    iter->num_kmers = length - K + 1;

    // Initialize rotation tables
    csyncmer_make_f_rot_64(S, iter->f_rot);
    csyncmer_make_c_rot_64(S, iter->c_rot);

    // Allocate power-of-2 circular buffer for canonical hashes
    size_t buf_size = 1;
    while (buf_size < iter->window_size) buf_size <<= 1;
    iter->buf_mask = buf_size - 1;
    iter->hash_buffer = (uint64_t*)malloc(buf_size * sizeof(uint64_t));
    iter->strand_buffer = (uint8_t*)malloc(buf_size * sizeof(uint8_t));
    if (!iter->hash_buffer || !iter->strand_buffer) {
        free(iter->hash_buffer);
        free(iter->strand_buffer);
        free(iter);
        return NULL;
    }

    // Use shared static tables
    const uint64_t* F_ASCII = CSYNCMER_F_ASCII_64;
    const uint64_t* RC_ASCII = CSYNCMER_RC_ASCII_64;
    const uint64_t* RC_ROTR7_ASCII = CSYNCMER_RC_ROTR7_ASCII_64;
    const uint8_t* IDX_ASCII = CSYNCMER_IDX_ASCII;

    // Compute first s-mer hash (both forward and RC)
    // Forward: process bases left to right with rotl7
    // RC: process bases right to left with rotl7 (same rotation, different order)
    uint64_t fw_hash = 0;
    uint64_t rc_hash = 0;
    for (size_t j = 0; j < S; ++j) {
        fw_hash = csyncmer_rotl7_64(fw_hash) ^ F_ASCII[iter->seq[j]];
        rc_hash = csyncmer_rotl7_64(rc_hash) ^ RC_ASCII[iter->seq[S - 1 - j]];
    }

    // Canonical hash = min(fw, rc), track strand
    // Fused comparison: single compare for both hash selection and strand
    int is_fw = (fw_hash <= rc_hash);
    uint64_t canon_hash = is_fw ? fw_hash : rc_hash;
    uint8_t strand = !is_fw;

    iter->hash_buffer[0] = canon_hash;
    iter->strand_buffer[0] = strand;
    iter->fw_hash = fw_hash;
    iter->rc_hash = rc_hash;

    // Fill initial window with remaining s-mer hashes using rolling
    // Optimized: pre-computed rotr7(RC[base]) table eliminates one rotation
    // ILP: interleave fw/rc computations for parallel execution
    for (size_t i = 1; i < iter->window_size; ++i) {
        uint8_t old_first = iter->seq[i - 1];
        uint8_t new_last = iter->seq[i + S - 1];

        // Cache index lookups (computed once, used twice)
        uint8_t idx_old = IDX_ASCII[old_first];
        uint8_t idx_new = IDX_ASCII[new_last];

        // Forward rolling (start fw_state computation)
        uint64_t fw_state = fw_hash ^ iter->f_rot[idx_old];

        // RC rotation (start early for ILP - CPU can execute in parallel with fw)
        uint64_t rc_rotr = csyncmer_rotr7_64(rc_hash);

        // Complete forward rolling
        fw_hash = csyncmer_rotl7_64(fw_state) ^ F_ASCII[new_last];

        // Complete RC rolling (use pre-computed rotr7 table - no rotation needed!)
        rc_hash = rc_rotr ^ RC_ROTR7_ASCII[old_first] ^ iter->c_rot[idx_new];

        // Fused comparison: single compare for both hash selection and strand
        is_fw = (fw_hash <= rc_hash);
        canon_hash = is_fw ? fw_hash : rc_hash;
        strand = !is_fw;

        iter->hash_buffer[i & iter->buf_mask] = canon_hash;
        iter->strand_buffer[i & iter->buf_mask] = strand;
    }

    iter->fw_hash = fw_hash;
    iter->rc_hash = rc_hash;

    // Find minimum in first window
    uint64_t min_hash = UINT64_MAX;
    size_t min_pos = 0;
    for (size_t i = 0; i < iter->window_size; ++i) {
        uint64_t h = iter->hash_buffer[i];
        if (h < min_hash) {
            min_hash = h;
            min_pos = i;
        }
    }
    iter->min_hash = min_hash;
    iter->min_pos = min_pos;
    iter->kmer_idx = 0;

    return iter;
}

// Get next syncmer position and strand. Returns 1 if valid, 0 when exhausted.
// pos: output for syncmer position in sequence
// strand: output for which strand had minimal s-mer (0=forward, 1=RC)
static inline int csyncmer_iterator_next_canonical_64(
    CsyncmerIteratorCanonical64* iter,
    size_t* pos,
    int* strand
) {
    if (!iter) return 0;

    // Cache hot state in local variables for register allocation
    const uint8_t* seq = iter->seq;
    size_t kmer_idx = iter->kmer_idx;
    const size_t num_kmers = iter->num_kmers;
    const size_t window_size = iter->window_size;
    const size_t S = iter->S;
    uint64_t fw_hash = iter->fw_hash;
    uint64_t rc_hash = iter->rc_hash;
    uint64_t min_hash = iter->min_hash;
    size_t min_pos = iter->min_pos;
    const size_t buf_mask = iter->buf_mask;
    uint64_t* hash_buffer = iter->hash_buffer;
    uint8_t* strand_buffer = iter->strand_buffer;
    const uint64_t* f_rot = iter->f_rot;
    const uint64_t* c_rot = iter->c_rot;

    // Use shared static tables
    const uint64_t* F_ASCII = CSYNCMER_F_ASCII_64;
    const uint64_t* RC_ROTR7_ASCII = CSYNCMER_RC_ROTR7_ASCII_64;
    const uint8_t* IDX_ASCII = CSYNCMER_IDX_ASCII;

    // Check first k-mer on first call
    if (kmer_idx == 0) {
        size_t min_offset = min_pos;
        if (min_offset == 0 || min_offset == window_size - 1) {
            *pos = 0;
            *strand = strand_buffer[min_pos & buf_mask];
            iter->kmer_idx = 1;
            return 1;
        }
        kmer_idx = 1;
    }

    // Main loop - optimized with pre-computed rotr7 and ILP restructuring
    while (kmer_idx < num_kmers) {
        size_t i = kmer_idx + window_size - 1;  // s-mer position entering window

        // Rolling hash computation
        uint8_t old_first = seq[i - 1];
        uint8_t new_last = seq[i + S - 1];

        // Cache index lookups (computed once, used twice)
        uint8_t idx_old = IDX_ASCII[old_first];
        uint8_t idx_new = IDX_ASCII[new_last];

        // Forward rolling (start fw_state computation)
        uint64_t fw_state = fw_hash ^ f_rot[idx_old];

        // RC rotation (start early for ILP - CPU can execute in parallel with fw)
        uint64_t rc_rotr = csyncmer_rotr7_64(rc_hash);

        // Complete forward rolling
        fw_hash = csyncmer_rotl7_64(fw_state) ^ F_ASCII[new_last];

        // Complete RC rolling (use pre-computed rotr7 table - no rotation needed!)
        rc_hash = rc_rotr ^ RC_ROTR7_ASCII[old_first] ^ c_rot[idx_new];

        // Fused comparison: single compare for both hash selection and strand
        int is_fw = (fw_hash <= rc_hash);
        uint64_t canon_hash = is_fw ? fw_hash : rc_hash;
        uint8_t cur_strand = !is_fw;

        hash_buffer[i & buf_mask] = canon_hash;
        strand_buffer[i & buf_mask] = cur_strand;

        // RESCAN: update minimum
        if (min_pos < kmer_idx) {
            // Minimum fell out of window - rescan
            min_hash = UINT64_MAX;
            for (size_t j = kmer_idx; j <= i; ++j) {
                uint64_t h = hash_buffer[j & buf_mask];
                if (h < min_hash) {
                    min_hash = h;
                    min_pos = j;
                }
            }
        } else {
            // Check if new hash is smaller
            if (canon_hash < min_hash) {
                min_hash = canon_hash;
                min_pos = i;
            }
        }

        // Check closed syncmer condition
        size_t min_offset = min_pos - kmer_idx;
        size_t current_kmer = kmer_idx;
        kmer_idx++;

        if (min_offset == 0 || min_offset == window_size - 1) {
            // Write back modified state
            iter->fw_hash = fw_hash;
            iter->rc_hash = rc_hash;
            iter->min_hash = min_hash;
            iter->min_pos = min_pos;
            iter->kmer_idx = kmer_idx;
            *pos = current_kmer;
            *strand = strand_buffer[min_pos & buf_mask];
            return 1;
        }
    }

    // Write back state on exhaustion
    iter->fw_hash = fw_hash;
    iter->rc_hash = rc_hash;
    iter->min_hash = min_hash;
    iter->min_pos = min_pos;
    iter->kmer_idx = kmer_idx;
    return 0;  // Exhausted
}

// Free canonical iterator
static inline void csyncmer_iterator_destroy_canonical_64(CsyncmerIteratorCanonical64* iter) {
    if (iter) {
        free(iter->hash_buffer);
        free(iter->strand_buffer);
        free(iter);
    }
}

// ============================================================================
// Optimized Two-Stack with SIMD 2-bit Packing and Parallel Hash
#ifdef __AVX2__
// Uses techniques from simd-minimizers:
// - 2-bit DNA packing (4 bases/byte)
// - SIMD table lookup for parallel hash computation
// - SIMD gather for loading all 8 chunks simultaneously
// - Split 32-bit hash/position buffers for exact correctness
// ============================================================================

// 2-bit pack table: A=0, C=1, T=2, G=3
// Note: Uses the same layout as base_to_bits (A=65, a=97, etc.)
static inline uint8_t csyncmer_pack_base(char c) {
    // Simple branchless lookup: A/a=0, C/c=1, T/t=2, G/g=3
    static const uint8_t table[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 0-15
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 16-31
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 32-47
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 48-63
        0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,0, // 64-79:  A=65->0, C=67->1, G=71->3
        0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0, // 80-95:  T=84->2
        0,0,0,1,0,0,0,3,0,0,0,0,0,0,0,0, // 96-111: a=97->0, c=99->1, g=103->3
        0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0, // 112-127: t=116->2
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 128-143
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 144-159
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 160-175
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 176-191
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 192-207
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 208-223
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 224-239
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  // 240-255
    };
    return table[(uint8_t)c];
}

#ifdef __BMI2__
// Pack 8 ASCII chars to 2 bytes using BMI2 PEXT
// Extracts bits 1-2 from each byte: A=0, C=1, T=2, G=3
static inline uint16_t csyncmer_pack_8_pext(const uint8_t* ascii) {
    uint64_t chars;
    memcpy(&chars, ascii, 8);
    return (uint16_t)_pext_u64(chars, 0x0606060606060606ULL);
}
#endif

// SIMD hash table (8 elements for permutevar8x32, duplicated)
// Using lower 32 bits of 64-bit ntHash constants
alignas(32) static const uint32_t CSYNCMER_SIMD_F32[8] = {
    0x95c60474,  // A (from 0x3c8bfbb395c60474)
    0x62a02b4c,  // C (from 0x3193c18562a02b4c)
    0x82572324,  // T (from 0x2032ed0082572324)
    0x4be24456,  // G (from 0x295549f54be24456)
    0x95c60474,  // A (duplicate for permutevar)
    0x62a02b4c,  // C
    0x82572324,  // T
    0x4be24456   // G
};

// 32-bit RC (reverse complement) constants: RC[base] = F[complement(base)]
// A↔T (0↔2), C↔G (1↔3)
alignas(32) static const uint32_t CSYNCMER_SIMD_RC32[8] = {
    0x82572324,  // RC[A] (idx 0) = F[T] = F[2]
    0x4be24456,  // RC[C] (idx 1) = F[G] = F[3]
    0x95c60474,  // RC[T] (idx 2) = F[A] = F[0]
    0x62a02b4c,  // RC[G] (idx 3) = F[C] = F[1]
    0x82572324,  // RC[A] (duplicate for permutevar)
    0x4be24456,  // RC[C]
    0x95c60474,  // RC[T]
    0x62a02b4c   // RC[G]
};

// Pre-computed rotr7(RC[base]) for optimized canonical SIMD computation
// Eliminates one rotation per iteration in RC rolling formula
// rotr7(x) = (x >> 7) | (x << 25) for 32-bit
alignas(32) static const uint32_t CSYNCMER_SIMD_RC32_ROTR7[8] = {
    0x4904ae46,  // rotr7(RC[A]) = rotr7(0x82572324)
    0xac97c488,  // rotr7(RC[C]) = rotr7(0x4be24456)
    0xe92b8c08,  // rotr7(RC[T]) = rotr7(0x95c60474)
    0x98c54056,  // rotr7(RC[G]) = rotr7(0x62a02b4c)
    0x4904ae46,  // duplicate for permutevar
    0xac97c488,
    0xe92b8c08,
    0x98c54056
};

// SIMD table lookup using permutevar8x32
static inline __m256i csyncmer_simd_lookup(__m256i table, __m256i indices) {
    return _mm256_permutevar8x32_epi32(table, indices);
}

// SIMD rotate left by 7 bits (32-bit elements)
static inline __m256i csyncmer_simd_rotl7(__m256i x) {
    return _mm256_or_si256(_mm256_slli_epi32(x, 7), _mm256_srli_epi32(x, 25));
}

// SIMD rotate right by 7 bits (32-bit elements) - for canonical RC rolling
static inline __m256i csyncmer_simd_rotr7(__m256i x) {
    return _mm256_or_si256(_mm256_srli_epi32(x, 7), _mm256_slli_epi32(x, 25));
}

// Pack sequence to 2-bit representation (4 bases per byte)
static inline void csyncmer_pack_seq_2bit(
    const char* seq,
    size_t len,
    uint8_t* out
) {
#ifdef __BMI2__
    // Fast path: pack 8 chars (2 output bytes) at a time using PEXT
    size_t i = 0;
    size_t byte_idx = 0;
    const uint8_t* useq = (const uint8_t*)seq;

    for (; i + 8 <= len; i += 8, byte_idx += 2) {
        uint16_t packed = csyncmer_pack_8_pext(useq + i);
        memcpy(out + byte_idx, &packed, 2);
    }

    // Remainder (0-7 chars)
    for (; i < len; i++) {
        size_t bidx = i / 4;
        size_t bit_offset = (i % 4) * 2;
        out[bidx] |= (csyncmer_pack_base(seq[i]) << bit_offset);
    }
#else
    // Scalar fallback
    for (size_t i = 0; i < len; i++) {
        size_t byte_idx = i / 4;
        size_t bit_offset = (i % 4) * 2;
        out[byte_idx] |= (csyncmer_pack_base(seq[i]) << bit_offset);
    }
#endif
}

// 8x8 matrix transpose for batch output collection (non-canonical positions)
static inline void csyncmer_transpose_8x8(__m256i m[8], __m256i t[8]) {
    __m256 r0 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[0]), _mm256_castsi256_ps(m[1]));
    __m256 r1 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[0]), _mm256_castsi256_ps(m[1]));
    __m256 r2 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[2]), _mm256_castsi256_ps(m[3]));
    __m256 r3 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[2]), _mm256_castsi256_ps(m[3]));
    __m256 r4 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[4]), _mm256_castsi256_ps(m[5]));
    __m256 r5 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[4]), _mm256_castsi256_ps(m[5]));
    __m256 r6 = _mm256_unpacklo_ps(_mm256_castsi256_ps(m[6]), _mm256_castsi256_ps(m[7]));
    __m256 r7 = _mm256_unpackhi_ps(_mm256_castsi256_ps(m[6]), _mm256_castsi256_ps(m[7]));

    __m256 s0 = _mm256_shuffle_ps(r0, r2, 0x44);
    __m256 s1 = _mm256_shuffle_ps(r0, r2, 0xEE);
    __m256 s2 = _mm256_shuffle_ps(r1, r3, 0x44);
    __m256 s3 = _mm256_shuffle_ps(r1, r3, 0xEE);
    __m256 s4 = _mm256_shuffle_ps(r4, r6, 0x44);
    __m256 s5 = _mm256_shuffle_ps(r4, r6, 0xEE);
    __m256 s6 = _mm256_shuffle_ps(r5, r7, 0x44);
    __m256 s7 = _mm256_shuffle_ps(r5, r7, 0xEE);

    t[0] = _mm256_castps_si256(_mm256_permute2f128_ps(s0, s4, 0x20));
    t[1] = _mm256_castps_si256(_mm256_permute2f128_ps(s1, s5, 0x20));
    t[2] = _mm256_castps_si256(_mm256_permute2f128_ps(s2, s6, 0x20));
    t[3] = _mm256_castps_si256(_mm256_permute2f128_ps(s3, s7, 0x20));
    t[4] = _mm256_castps_si256(_mm256_permute2f128_ps(s0, s4, 0x31));
    t[5] = _mm256_castps_si256(_mm256_permute2f128_ps(s1, s5, 0x31));
    t[6] = _mm256_castps_si256(_mm256_permute2f128_ps(s2, s6, 0x31));
    t[7] = _mm256_castps_si256(_mm256_permute2f128_ps(s3, s7, 0x31));
}

// Pre-computed dedup shuffle table (256 entries for all 8-bit masks)
// Entry i contains indices to pack non-masked elements to the left
// Using Lemire's algorithm for branchless filtering
alignas(32) static const uint32_t CSYNCMER_UNIQSHUF[256][8] = {
    {0,1,2,3,4,5,6,7}, {1,2,3,4,5,6,7,0}, {0,2,3,4,5,6,7,1}, {2,3,4,5,6,7,0,1},
    {0,1,3,4,5,6,7,2}, {1,3,4,5,6,7,0,2}, {0,3,4,5,6,7,1,2}, {3,4,5,6,7,0,1,2},
    {0,1,2,4,5,6,7,3}, {1,2,4,5,6,7,0,3}, {0,2,4,5,6,7,1,3}, {2,4,5,6,7,0,1,3},
    {0,1,4,5,6,7,2,3}, {1,4,5,6,7,0,2,3}, {0,4,5,6,7,1,2,3}, {4,5,6,7,0,1,2,3},
    {0,1,2,3,5,6,7,4}, {1,2,3,5,6,7,0,4}, {0,2,3,5,6,7,1,4}, {2,3,5,6,7,0,1,4},
    {0,1,3,5,6,7,2,4}, {1,3,5,6,7,0,2,4}, {0,3,5,6,7,1,2,4}, {3,5,6,7,0,1,2,4},
    {0,1,2,5,6,7,3,4}, {1,2,5,6,7,0,3,4}, {0,2,5,6,7,1,3,4}, {2,5,6,7,0,1,3,4},
    {0,1,5,6,7,2,3,4}, {1,5,6,7,0,2,3,4}, {0,5,6,7,1,2,3,4}, {5,6,7,0,1,2,3,4},
    {0,1,2,3,4,6,7,5}, {1,2,3,4,6,7,0,5}, {0,2,3,4,6,7,1,5}, {2,3,4,6,7,0,1,5},
    {0,1,3,4,6,7,2,5}, {1,3,4,6,7,0,2,5}, {0,3,4,6,7,1,2,5}, {3,4,6,7,0,1,2,5},
    {0,1,2,4,6,7,3,5}, {1,2,4,6,7,0,3,5}, {0,2,4,6,7,1,3,5}, {2,4,6,7,0,1,3,5},
    {0,1,4,6,7,2,3,5}, {1,4,6,7,0,2,3,5}, {0,4,6,7,1,2,3,5}, {4,6,7,0,1,2,3,5},
    {0,1,2,3,6,7,4,5}, {1,2,3,6,7,0,4,5}, {0,2,3,6,7,1,4,5}, {2,3,6,7,0,1,4,5},
    {0,1,3,6,7,2,4,5}, {1,3,6,7,0,2,4,5}, {0,3,6,7,1,2,4,5}, {3,6,7,0,1,2,4,5},
    {0,1,2,6,7,3,4,5}, {1,2,6,7,0,3,4,5}, {0,2,6,7,1,3,4,5}, {2,6,7,0,1,3,4,5},
    {0,1,6,7,2,3,4,5}, {1,6,7,0,2,3,4,5}, {0,6,7,1,2,3,4,5}, {6,7,0,1,2,3,4,5},
    {0,1,2,3,4,5,7,6}, {1,2,3,4,5,7,0,6}, {0,2,3,4,5,7,1,6}, {2,3,4,5,7,0,1,6},
    {0,1,3,4,5,7,2,6}, {1,3,4,5,7,0,2,6}, {0,3,4,5,7,1,2,6}, {3,4,5,7,0,1,2,6},
    {0,1,2,4,5,7,3,6}, {1,2,4,5,7,0,3,6}, {0,2,4,5,7,1,3,6}, {2,4,5,7,0,1,3,6},
    {0,1,4,5,7,2,3,6}, {1,4,5,7,0,2,3,6}, {0,4,5,7,1,2,3,6}, {4,5,7,0,1,2,3,6},
    {0,1,2,3,5,7,4,6}, {1,2,3,5,7,0,4,6}, {0,2,3,5,7,1,4,6}, {2,3,5,7,0,1,4,6},
    {0,1,3,5,7,2,4,6}, {1,3,5,7,0,2,4,6}, {0,3,5,7,1,2,4,6}, {3,5,7,0,1,2,4,6},
    {0,1,2,5,7,3,4,6}, {1,2,5,7,0,3,4,6}, {0,2,5,7,1,3,4,6}, {2,5,7,0,1,3,4,6},
    {0,1,5,7,2,3,4,6}, {1,5,7,0,2,3,4,6}, {0,5,7,1,2,3,4,6}, {5,7,0,1,2,3,4,6},
    {0,1,2,3,4,7,5,6}, {1,2,3,4,7,0,5,6}, {0,2,3,4,7,1,5,6}, {2,3,4,7,0,1,5,6},
    {0,1,3,4,7,2,5,6}, {1,3,4,7,0,2,5,6}, {0,3,4,7,1,2,5,6}, {3,4,7,0,1,2,5,6},
    {0,1,2,4,7,3,5,6}, {1,2,4,7,0,3,5,6}, {0,2,4,7,1,3,5,6}, {2,4,7,0,1,3,5,6},
    {0,1,4,7,2,3,5,6}, {1,4,7,0,2,3,5,6}, {0,4,7,1,2,3,5,6}, {4,7,0,1,2,3,5,6},
    {0,1,2,3,7,4,5,6}, {1,2,3,7,0,4,5,6}, {0,2,3,7,1,4,5,6}, {2,3,7,0,1,4,5,6},
    {0,1,3,7,2,4,5,6}, {1,3,7,0,2,4,5,6}, {0,3,7,1,2,4,5,6}, {3,7,0,1,2,4,5,6},
    {0,1,2,7,3,4,5,6}, {1,2,7,0,3,4,5,6}, {0,2,7,1,3,4,5,6}, {2,7,0,1,3,4,5,6},
    {0,1,7,2,3,4,5,6}, {1,7,0,2,3,4,5,6}, {0,7,1,2,3,4,5,6}, {7,0,1,2,3,4,5,6},
    {0,1,2,3,4,5,6,7}, {1,2,3,4,5,6,0,7}, {0,2,3,4,5,6,1,7}, {2,3,4,5,6,0,1,7},
    {0,1,3,4,5,6,2,7}, {1,3,4,5,6,0,2,7}, {0,3,4,5,6,1,2,7}, {3,4,5,6,0,1,2,7},
    {0,1,2,4,5,6,3,7}, {1,2,4,5,6,0,3,7}, {0,2,4,5,6,1,3,7}, {2,4,5,6,0,1,3,7},
    {0,1,4,5,6,2,3,7}, {1,4,5,6,0,2,3,7}, {0,4,5,6,1,2,3,7}, {4,5,6,0,1,2,3,7},
    {0,1,2,3,5,6,4,7}, {1,2,3,5,6,0,4,7}, {0,2,3,5,6,1,4,7}, {2,3,5,6,0,1,4,7},
    {0,1,3,5,6,2,4,7}, {1,3,5,6,0,2,4,7}, {0,3,5,6,1,2,4,7}, {3,5,6,0,1,2,4,7},
    {0,1,2,5,6,3,4,7}, {1,2,5,6,0,3,4,7}, {0,2,5,6,1,3,4,7}, {2,5,6,0,1,3,4,7},
    {0,1,5,6,2,3,4,7}, {1,5,6,0,2,3,4,7}, {0,5,6,1,2,3,4,7}, {5,6,0,1,2,3,4,7},
    {0,1,2,3,4,6,5,7}, {1,2,3,4,6,0,5,7}, {0,2,3,4,6,1,5,7}, {2,3,4,6,0,1,5,7},
    {0,1,3,4,6,2,5,7}, {1,3,4,6,0,2,5,7}, {0,3,4,6,1,2,5,7}, {3,4,6,0,1,2,5,7},
    {0,1,2,4,6,3,5,7}, {1,2,4,6,0,3,5,7}, {0,2,4,6,1,3,5,7}, {2,4,6,0,1,3,5,7},
    {0,1,4,6,2,3,5,7}, {1,4,6,0,2,3,5,7}, {0,4,6,1,2,3,5,7}, {4,6,0,1,2,3,5,7},
    {0,1,2,3,6,4,5,7}, {1,2,3,6,0,4,5,7}, {0,2,3,6,1,4,5,7}, {2,3,6,0,1,4,5,7},
    {0,1,3,6,2,4,5,7}, {1,3,6,0,2,4,5,7}, {0,3,6,1,2,4,5,7}, {3,6,0,1,2,4,5,7},
    {0,1,2,6,3,4,5,7}, {1,2,6,0,3,4,5,7}, {0,2,6,1,3,4,5,7}, {2,6,0,1,3,4,5,7},
    {0,1,6,2,3,4,5,7}, {1,6,0,2,3,4,5,7}, {0,6,1,2,3,4,5,7}, {6,0,1,2,3,4,5,7},
    {0,1,2,3,4,5,6,7}, {1,2,3,4,5,0,6,7}, {0,2,3,4,5,1,6,7}, {2,3,4,5,0,1,6,7},
    {0,1,3,4,5,2,6,7}, {1,3,4,5,0,2,6,7}, {0,3,4,5,1,2,6,7}, {3,4,5,0,1,2,6,7},
    {0,1,2,4,5,3,6,7}, {1,2,4,5,0,3,6,7}, {0,2,4,5,1,3,6,7}, {2,4,5,0,1,3,6,7},
    {0,1,4,5,2,3,6,7}, {1,4,5,0,2,3,6,7}, {0,4,5,1,2,3,6,7}, {4,5,0,1,2,3,6,7},
    {0,1,2,3,5,4,6,7}, {1,2,3,5,0,4,6,7}, {0,2,3,5,1,4,6,7}, {2,3,5,0,1,4,6,7},
    {0,1,3,5,2,4,6,7}, {1,3,5,0,2,4,6,7}, {0,3,5,1,2,4,6,7}, {3,5,0,1,2,4,6,7},
    {0,1,2,5,3,4,6,7}, {1,2,5,0,3,4,6,7}, {0,2,5,1,3,4,6,7}, {2,5,0,1,3,4,6,7},
    {0,1,5,2,3,4,6,7}, {1,5,0,2,3,4,6,7}, {0,5,1,2,3,4,6,7}, {5,0,1,2,3,4,6,7},
    {0,1,2,3,4,5,6,7}, {1,2,3,4,0,5,6,7}, {0,2,3,4,1,5,6,7}, {2,3,4,0,1,5,6,7},
    {0,1,3,4,2,5,6,7}, {1,3,4,0,2,5,6,7}, {0,3,4,1,2,5,6,7}, {3,4,0,1,2,5,6,7},
    {0,1,2,4,3,5,6,7}, {1,2,4,0,3,5,6,7}, {0,2,4,1,3,5,6,7}, {2,4,0,1,3,5,6,7},
    {0,1,4,2,3,5,6,7}, {1,4,0,2,3,5,6,7}, {0,4,1,2,3,5,6,7}, {4,0,1,2,3,5,6,7},
    {0,1,2,3,4,5,6,7}, {1,2,3,0,4,5,6,7}, {0,2,3,1,4,5,6,7}, {2,3,0,1,4,5,6,7},
    {0,1,3,2,4,5,6,7}, {1,3,0,2,4,5,6,7}, {0,3,1,2,4,5,6,7}, {3,0,1,2,4,5,6,7},
    {0,1,2,3,4,5,6,7}, {1,2,0,3,4,5,6,7}, {0,2,1,3,4,5,6,7}, {2,0,1,3,4,5,6,7},
    {0,1,2,3,4,5,6,7}, {1,0,2,3,4,5,6,7}, {0,1,2,3,4,5,6,7}, {0,1,2,3,4,5,6,7},
};

// Filtered append using dedup shuffle table
// skip_mask: bit i=1 means skip lane i
// Returns new write index
static inline size_t csyncmer_append_filtered(
    __m256i vals, int skip_mask, uint32_t* out, size_t write_idx
) {
    int num = 8 - __builtin_popcount(skip_mask);
    if (num == 0) return write_idx;
    __m256i shuf = _mm256_load_si256((__m256i*)CSYNCMER_UNIQSHUF[skip_mask]);
    __m256i packed = _mm256_permutevar8x32_epi32(vals, shuf);
    _mm256_storeu_si256((__m256i*)(out + write_idx), packed);
    return write_idx + num;
}

// Optimized two-stack with SIMD hash computation (32-bit hash, 8 parallel lanes)
// Macro-generated for forward/canonical and count/positions variants.
//
// LIMITATIONS:
// - Uses 16-bit hash packing, causing ~0.00004% discrepancy vs reference when
//   two s-mer hashes have identical upper 16 bits.
// - Falls back to scalar RESCAN when window_size > 64 (i.e., K - S + 1 > 64)
//   or num_kmers < 64. The fallback uses full 32-bit hash resolution, so results
//   may differ from the SIMD path in the rare case of 16-bit hash collisions.
//   This is inconsequential for the num_kmers < 64 case (too few s-mers for
//   collisions), but could in theory matter for the window_size > 64 case on
//   long sequences.
// For exact results, use csyncmer_rescan_32_count() or the 64-bit
// iterator API (csyncmer_iterator_*_64).

#define CSYNCMER_DEFINE_TWOSTACK_SIMD_32(FUNC_NAME, CANONICAL, COLLECT_POSITIONS) \
static inline size_t FUNC_NAME(                                                \
    const char* sequence,                                                      \
    size_t length,                                                             \
    size_t K,                                                                  \
    size_t S,                                                                  \
    uint32_t* out_positions,                                                   \
    uint8_t* out_strands,                                                      \
    size_t max_positions                                                        \
) {                                                                            \
    if (!sequence || length < K || S == 0 || S >= K) return 0;                 \
    if (COLLECT_POSITIONS && max_positions == 0) return 0;                      \
                                                                               \
    size_t window_size = K - S + 1;                                            \
    size_t num_smers = length - S + 1;                                         \
    size_t num_kmers = num_smers - window_size + 1;                            \
                                                                               \
    /* Fall back to scalar for small inputs or large windows */                \
    if (num_kmers < 64 || window_size > 64) {                                 \
        if (CANONICAL) {                                                       \
            if (COLLECT_POSITIONS) {                                           \
                return csyncmer_canonical_rescan_32_(                           \
                    sequence, length, K, S, out_positions, out_strands,         \
                    max_positions);                                             \
            } else {                                                           \
                return csyncmer_canonical_rescan_32_(                           \
                    sequence, length, K, S, NULL, NULL, 0);                    \
            }                                                                  \
        } else {                                                               \
            if (COLLECT_POSITIONS) {                                           \
                return csyncmer_rescan_32_(                                     \
                    sequence, length, K, S, out_positions, NULL, max_positions);\
            } else {                                                           \
                return csyncmer_rescan_32_(                                     \
                    sequence, length, K, S, NULL, NULL, 0);                    \
            }                                                                  \
        }                                                                      \
    }                                                                          \
                                                                               \
    /* Phase 1: Pack sequence to 2-bit representation */                       \
    size_t packed_size = (length + 3) / 4 + 32;                                \
    uint8_t* packed = (uint8_t*)aligned_alloc(32, packed_size);                \
    if (!packed) return 0;                                                     \
    memset(packed, 0, packed_size);                                            \
    csyncmer_pack_seq_2bit(sequence, length, packed);                          \
                                                                               \
    /* Phase 2: Calculate chunk parameters */                                  \
    size_t chunk_len = (num_kmers + 7) / 8;                                    \
    chunk_len = ((chunk_len + 3) / 4) * 4;                                     \
    size_t bytes_per_chunk = chunk_len / 4;                                    \
    size_t max_iter = chunk_len + window_size - 1 + S - 1;                     \
                                                                               \
    alignas(32) uint32_t chunk_starts[8];                                      \
    alignas(32) uint32_t chunk_ends[8];                                        \
    for (int i = 0; i < 8; i++) {                                              \
        chunk_starts[i] = (uint32_t)(i * chunk_len);                           \
        chunk_ends[i] = (uint32_t)((i + 1) * chunk_len);                      \
        if (chunk_ends[i] > num_kmers) chunk_ends[i] = (uint32_t)num_kmers;   \
    }                                                                          \
                                                                               \
    /* Thread-local lane buffers (unique per expansion via ##) */              \
    static CSYNCMER_THREAD_LOCAL uint32_t* ts_pos_bufs_##FUNC_NAME[8];        \
    static CSYNCMER_THREAD_LOCAL uint8_t* ts_str_bufs_##FUNC_NAME[8];         \
    static CSYNCMER_THREAD_LOCAL size_t ts_buf_cap_##FUNC_NAME[8];            \
    static CSYNCMER_THREAD_LOCAL int ts_init_##FUNC_NAME;                     \
    size_t lane_counts[8] = {0};                                               \
    size_t max_per_lane = chunk_len + 64;                                      \
                                                                               \
    if (COLLECT_POSITIONS) {                                                   \
        if (!ts_init_##FUNC_NAME) {                                            \
            memset(ts_pos_bufs_##FUNC_NAME, 0, sizeof(ts_pos_bufs_##FUNC_NAME)); \
            memset(ts_str_bufs_##FUNC_NAME, 0, sizeof(ts_str_bufs_##FUNC_NAME)); \
            memset(ts_buf_cap_##FUNC_NAME, 0, sizeof(ts_buf_cap_##FUNC_NAME));\
            ts_init_##FUNC_NAME = 1;                                           \
        }                                                                      \
        for (int i = 0; i < 8; i++) {                                          \
            if (ts_buf_cap_##FUNC_NAME[i] < max_per_lane) {                   \
                free(ts_pos_bufs_##FUNC_NAME[i]);                              \
                ts_pos_bufs_##FUNC_NAME[i] = (uint32_t*)aligned_alloc(         \
                    32, max_per_lane * sizeof(uint32_t));                       \
                if (!ts_pos_bufs_##FUNC_NAME[i]) {                             \
                    ts_buf_cap_##FUNC_NAME[i] = 0;                             \
                    free(packed); return 0;                                     \
                }                                                              \
                if (CANONICAL) {                                               \
                    free(ts_str_bufs_##FUNC_NAME[i]);                          \
                    ts_str_bufs_##FUNC_NAME[i] = (uint8_t*)aligned_alloc(      \
                        32, max_per_lane);                                     \
                    if (!ts_str_bufs_##FUNC_NAME[i]) {                         \
                        free(ts_pos_bufs_##FUNC_NAME[i]);                      \
                        ts_pos_bufs_##FUNC_NAME[i] = NULL;                     \
                        ts_buf_cap_##FUNC_NAME[i] = 0;                         \
                        free(packed); return 0;                                \
                    }                                                          \
                }                                                              \
                ts_buf_cap_##FUNC_NAME[i] = max_per_lane;                      \
            }                                                                  \
        }                                                                      \
    }                                                                          \
                                                                               \
    size_t syncmer_count = 0;                                                  \
                                                                               \
    /* Phase 3: Initialize SIMD state */                                       \
    __m256i f_table = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_F32);    \
                                                                               \
    uint32_t rot = ((S - 1) * 7) & 31;                                        \
    alignas(32) uint32_t f_rot_arr[8];                                         \
    for (int i = 0; i < 8; i++) {                                              \
        uint32_t x = CSYNCMER_SIMD_F32[i];                                    \
        f_rot_arr[i] = (x << rot) | (x >> (32 - rot));                        \
    }                                                                          \
    __m256i f_rot_table = _mm256_load_si256((const __m256i*)f_rot_arr);        \
                                                                               \
    /* Canonical: RC tables */                                                 \
    __m256i rc_table_v = _mm256_setzero_si256();                               \
    __m256i rc_rotr7_table = _mm256_setzero_si256();                           \
    __m256i c_rot_table = _mm256_setzero_si256();                              \
    if (CANONICAL) {                                                           \
        rc_table_v = _mm256_load_si256((const __m256i*)CSYNCMER_SIMD_RC32);    \
        rc_rotr7_table = _mm256_load_si256(                                    \
            (const __m256i*)CSYNCMER_SIMD_RC32_ROTR7);                        \
        alignas(32) uint32_t c_rot_arr[8];                                     \
        for (int i = 0; i < 8; i++) {                                          \
            uint32_t cx = CSYNCMER_SIMD_RC32[i];                               \
            c_rot_arr[i] = (cx << rot) | (cx >> (32 - rot));                   \
        }                                                                      \
        c_rot_table = _mm256_load_si256((const __m256i*)c_rot_arr);            \
    }                                                                          \
                                                                               \
    __m256i fw = _mm256_setzero_si256();                                       \
    __m256i rc = _mm256_setzero_si256();                                       \
    __m256i prev_remove_base = _mm256_setzero_si256();                         \
                                                                               \
    /* Delay buffer */                                                         \
    size_t delay_size = 1;                                                     \
    while (delay_size < S) delay_size *= 2;                                    \
    size_t delay_mask = delay_size - 1;                                        \
    __m256i* delay_buf = (__m256i*)aligned_alloc(32,                           \
        delay_size * sizeof(__m256i));                                         \
    if (!delay_buf) { free(packed); return 0; }                                \
    for (size_t i = 0; i < delay_size; i++)                                    \
        delay_buf[i] = _mm256_setzero_si256();                                \
    size_t wr_idx = 0, rd_idx = 0;                                            \
                                                                               \
    /* Two-stack: [hash upper 16 bits][position lower 16 bits] */              \
    alignas(32) __m256i ring_buf[64];                                          \
    for (size_t i = 0; i < window_size; i++)                                   \
        ring_buf[i] = _mm256_set1_epi32((int)UINT32_MAX);                     \
    __m256i prefix_min = _mm256_set1_epi32((int)UINT32_MAX);                   \
    __m256i val_mask = _mm256_set1_epi32((int)0xFFFF0000);                     \
    __m256i pos_mask = _mm256_set1_epi32(0x0000FFFF);                          \
                                                                               \
    /* Canonical + positions: strand tracking */                               \
    alignas(32) uint8_t strand_ring[64];                                       \
    uint8_t prefix_strand = 0;                                                 \
    if (CANONICAL && COLLECT_POSITIONS) {                                       \
        for (size_t i = 0; i < window_size; i++) strand_ring[i] = 0;          \
    }                                                                          \
                                                                               \
    size_t ring_idx = 0;                                                       \
    uint32_t pos = 0;                                                          \
    uint32_t pos_offset = 0;                                                   \
    __m256i pos_offset_vec = _mm256_setzero_si256();                           \
    const uint32_t max_pos_val = 0xFFFF;                                       \
                                                                               \
                                                                               \
    /* Batch buffer for non-canonical positions (transpose approach) */         \
    alignas(32) __m256i batch_pos[8];                                          \
    size_t batch_count = 0;                                                    \
    size_t batch_base_kmer = 0;                                                \
    alignas(32) uint32_t idx_arr[8] = {0, 1, 2, 3, 4, 5, 6, 7};              \
    __m256i idx_offsets = _mm256_load_si256((const __m256i*)idx_arr);           \
    __m256i w_minus_1_vec = _mm256_set1_epi32((int)(window_size - 1));         \
                                                                               \
    __m256i mask_2bit = _mm256_set1_epi32(0x03);                               \
    alignas(32) int32_t gather_base[8];                                        \
    for (int i = 0; i < 8; i++)                                                \
        gather_base[i] = (int32_t)(i * bytes_per_chunk);                       \
    __m256i gather_base_vec = _mm256_load_si256((const __m256i*)gather_base);   \
                                                                               \
    /* Validity threshold: min lane length across all lanes.                   \
     * Must handle chunk_ends[i] <= chunk_starts[i] for short sequences        \
     * where later lanes have no valid k-mers. */                              \
    size_t last_lane_limit = chunk_len;                                        \
    for (int i = 0; i < 8; i++) {                                              \
        if (chunk_ends[i] > chunk_starts[i]) {                                 \
            size_t ll = chunk_ends[i] - chunk_starts[i];                       \
            if (ll < last_lane_limit) last_lane_limit = ll;                    \
        } else { last_lane_limit = 0; }                                        \
    }                                                                          \
                                                                               \
    __m256i cur_data = _mm256_setzero_si256();                                 \
    size_t buf_pos = 16;                                                       \
    size_t seq_pos = 0;                                                        \
                                                                               \
    /* === MAIN LOOP === */                                                    \
    for (size_t iter = 0; iter < max_iter; iter++) {                           \
        if (buf_pos >= 16) {                                                   \
            size_t byte_offset = seq_pos / 4;                                  \
            __m256i byte_off_vec = _mm256_set1_epi32((int32_t)byte_offset);    \
            __m256i bv_idx = _mm256_add_epi32(gather_base_vec, byte_off_vec);  \
            cur_data = _mm256_i32gather_epi32((const int*)packed, bv_idx, 1);  \
            buf_pos = seq_pos % 16;                                            \
        }                                                                      \
        __m256i add_base = _mm256_and_si256(                                   \
            _mm256_srli_epi32(cur_data, (int)(buf_pos * 2)), mask_2bit);       \
        buf_pos++;                                                             \
        seq_pos++;                                                             \
                                                                               \
        __m256i remove_base = delay_buf[rd_idx];                               \
        delay_buf[wr_idx] = add_base;                                          \
        wr_idx = (wr_idx + 1) & delay_mask;                                    \
        if (iter >= S - 1) rd_idx = (rd_idx + 1) & delay_mask;                \
                                                                               \
        /* Forward hash */                                                     \
        __m256i fw_rotated = csyncmer_simd_rotl7(fw);                          \
        __m256i add_hash_fw = csyncmer_simd_lookup(f_table, add_base);         \
        __m256i fw_hash = _mm256_xor_si256(fw_rotated, add_hash_fw);           \
                                                                               \
        __m256i hash_out;                                                      \
        uint8_t strand_mask = 0;                                               \
                                                                               \
        if (CANONICAL) {                                                       \
            __m256i rc_hash;                                                   \
            if (iter < S - 1) {                                                \
                fw = fw_hash;                                                  \
                prev_remove_base = remove_base;                                \
                continue;                                                      \
            } else if (iter == S - 1) {                                        \
                rc = _mm256_setzero_si256();                                   \
                for (size_t j = 0; j < S; j++) {                               \
                    size_t bi = (S - 1 - j) & delay_mask;                     \
                    __m256i bj = delay_buf[bi];                                \
                    __m256i rv = csyncmer_simd_lookup(rc_table_v, bj);         \
                    rc = _mm256_xor_si256(csyncmer_simd_rotl7(rc), rv);       \
                }                                                              \
                rc_hash = rc;                                                  \
                __m256i fw_rm = csyncmer_simd_lookup(f_rot_table, remove_base);\
                fw = _mm256_xor_si256(fw_hash, fw_rm);                        \
                prev_remove_base = remove_base;                                \
            } else {                                                           \
                __m256i rc_rotr = csyncmer_simd_rotr7(rc);                     \
                __m256i rc_rm = csyncmer_simd_lookup(                          \
                    rc_rotr7_table, prev_remove_base);                         \
                __m256i rc_ad = csyncmer_simd_lookup(c_rot_table, add_base);   \
                rc_hash = _mm256_xor_si256(                                    \
                    _mm256_xor_si256(rc_rotr, rc_rm), rc_ad);                 \
                __m256i fw_rm = csyncmer_simd_lookup(f_rot_table, remove_base);\
                fw = _mm256_xor_si256(fw_hash, fw_rm);                        \
                rc = rc_hash;                                                  \
                prev_remove_base = remove_base;                                \
            }                                                                  \
            hash_out = _mm256_min_epu32(fw_hash, rc_hash);                     \
            if (COLLECT_POSITIONS) {                                           \
                __m256i cmp = _mm256_cmpgt_epi32(                              \
                    _mm256_xor_si256(fw_hash,                                  \
                        _mm256_set1_epi32((int)0x80000000)),                   \
                    _mm256_xor_si256(rc_hash,                                  \
                        _mm256_set1_epi32((int)0x80000000)));                  \
                strand_mask = (uint8_t)_mm256_movemask_ps(                    \
                    _mm256_castsi256_ps(cmp));                                \
            }                                                                  \
        } else {                                                               \
            if (iter >= S - 1) {                                               \
                __m256i rm = csyncmer_simd_lookup(f_rot_table, remove_base);   \
                fw = _mm256_xor_si256(fw_hash, rm);                           \
            } else {                                                           \
                fw = fw_hash;                                                  \
                continue;                                                      \
            }                                                                  \
            hash_out = fw_hash;                                                \
        }                                                                      \
                                                                               \
        /* Pack hash + position */                                             \
        __m256i pos_vec = _mm256_set1_epi32((int)pos);                         \
        __m256i elem = _mm256_or_si256(                                        \
            _mm256_and_si256(hash_out, val_mask),                              \
            _mm256_and_si256(pos_vec, pos_mask));                              \
        ring_buf[ring_idx] = elem;                                             \
                                                                               \
        /* Update prefix minimum */                                            \
        if (CANONICAL && COLLECT_POSITIONS) {                                  \
            strand_ring[ring_idx] = strand_mask;                               \
            __m256i pcmp = _mm256_cmpgt_epi32(                                \
                _mm256_xor_si256(prefix_min,                                   \
                    _mm256_set1_epi32((int)0x80000000)),                       \
                _mm256_xor_si256(elem,                                         \
                    _mm256_set1_epi32((int)0x80000000)));                      \
            uint8_t umask = (uint8_t)_mm256_movemask_ps(                      \
                _mm256_castsi256_ps(pcmp));                                    \
            prefix_min = _mm256_min_epu32(prefix_min, elem);                   \
            prefix_strand = (prefix_strand & ~umask) |                        \
                            (strand_mask & umask);                             \
        } else {                                                               \
            prefix_min = _mm256_min_epu32(prefix_min, elem);                   \
        }                                                                      \
                                                                               \
        /* Handle 16-bit position overflow */                                  \
        if (pos == max_pos_val) {                                              \
            uint32_t delta = (1 << 16) - 2 - 2 * window_size;                 \
            pos -= delta;                                                      \
            pos_offset += delta;                                               \
            pos_offset_vec = _mm256_set1_epi32((int)pos_offset);              \
            __m256i dv = _mm256_set1_epi32((int)delta);                       \
            prefix_min = _mm256_sub_epi32(prefix_min, dv);                    \
            for (size_t j = 0; j < window_size; j++)                           \
                ring_buf[j] = _mm256_sub_epi32(ring_buf[j], dv);             \
        }                                                                      \
                                                                               \
        pos++;                                                                 \
        ring_idx++;                                                            \
                                                                               \
        /* Suffix recomputation when ring wraps */                             \
        if (ring_idx == window_size) {                                         \
            ring_idx = 0;                                                      \
            __m256i smin = ring_buf[window_size - 1];                          \
            if (CANONICAL && COLLECT_POSITIONS) {                              \
                uint8_t sstr = strand_ring[window_size - 1];                   \
                for (size_t j = window_size - 1; j > 0; j--) {                \
                    __m256i prev = ring_buf[j - 1];                            \
                    __m256i sc = _mm256_cmpgt_epi32(                           \
                        _mm256_xor_si256(smin,                                 \
                            _mm256_set1_epi32((int)0x80000000)),               \
                        _mm256_xor_si256(prev,                                 \
                            _mm256_set1_epi32((int)0x80000000)));              \
                    uint8_t sm = (uint8_t)_mm256_movemask_ps(                 \
                        _mm256_castsi256_ps(sc));                              \
                    smin = _mm256_min_epu32(smin, prev);                       \
                    sstr = (sstr & ~sm) | (strand_ring[j - 1] & sm);          \
                    ring_buf[j - 1] = smin;                                    \
                    strand_ring[j - 1] = sstr;                                 \
                }                                                              \
            } else {                                                           \
                for (size_t j = window_size - 1; j > 0; j--) {                \
                    smin = _mm256_min_epu32(smin, ring_buf[j - 1]);            \
                    ring_buf[j - 1] = smin;                                    \
                }                                                              \
            }                                                                  \
            prefix_min = _mm256_set1_epi32((int)UINT32_MAX);                   \
            if (CANONICAL && COLLECT_POSITIONS) prefix_strand = 0;             \
        }                                                                      \
                                                                               \
        /* === SYNCMER CHECK === */                                            \
        if (iter < S - 1 + window_size - 1) continue;                         \
                                                                               \
        __m256i suf_min = ring_buf[ring_idx];                                  \
        __m256i min_elem;                                                      \
        uint8_t overall_strand = 0;                                            \
                                                                               \
        if (CANONICAL && COLLECT_POSITIONS) {                                  \
            uint8_t suf_str = strand_ring[ring_idx];                           \
            __m256i oc = _mm256_cmpgt_epi32(                                  \
                _mm256_xor_si256(prefix_min,                                   \
                    _mm256_set1_epi32((int)0x80000000)),                       \
                _mm256_xor_si256(suf_min,                                      \
                    _mm256_set1_epi32((int)0x80000000)));                      \
            uint8_t om = (uint8_t)_mm256_movemask_ps(                         \
                _mm256_castsi256_ps(oc));                                      \
            min_elem = _mm256_min_epu32(prefix_min, suf_min);                  \
            overall_strand = (prefix_strand & ~om) | (suf_str & om);          \
        } else {                                                               \
            min_elem = _mm256_min_epu32(prefix_min, suf_min);                  \
        }                                                                      \
                                                                               \
        __m256i min_pos_vec = _mm256_add_epi32(                                \
            _mm256_and_si256(min_elem, pos_mask), pos_offset_vec);             \
        size_t kmer_idx = iter - S + 1 - window_size + 1;                     \
                                                                               \
        if (COLLECT_POSITIONS && !CANONICAL) {                                  \
            /* Non-canonical positions: batch + transpose approach.             \
             * Keeps the syncmer check off the hot loop critical path. */       \
            if (batch_count % 8 == 0) batch_base_kmer = kmer_idx;             \
            batch_pos[batch_count % 8] = min_pos_vec;                          \
            batch_count++;                                                     \
                                                                               \
            if (batch_count % 8 == 0) {                                        \
                __m256i tp[8];                                                 \
                csyncmer_transpose_8x8(batch_pos, tp);                         \
                __m256i bb = _mm256_set1_epi32((int)batch_base_kmer);          \
                __m256i bfs = _mm256_add_epi32(bb, idx_offsets);               \
                __m256i bls = _mm256_add_epi32(bfs, w_minus_1_vec);            \
                for (int lane = 0; lane < 8; lane++) {                         \
                    __m256i lo = _mm256_set1_epi32(                            \
                        (int)chunk_starts[lane]);                              \
                    __m256i ap = _mm256_add_epi32(lo, bfs);                    \
                    __m256i bisy = _mm256_or_si256(                            \
                        _mm256_cmpeq_epi32(tp[lane], bfs),                     \
                        _mm256_cmpeq_epi32(tp[lane], bls));                    \
                    if (batch_base_kmer + 7 >= last_lane_limit) {              \
                        __m256i lim = _mm256_set1_epi32(                       \
                            (int)(chunk_ends[lane] - chunk_starts[lane]));     \
                        bisy = _mm256_and_si256(bisy,                          \
                            _mm256_cmpgt_epi32(lim, _mm256_add_epi32(         \
                                bb, idx_offsets)));                            \
                    }                                                          \
                    int bkm = _mm256_movemask_ps(                              \
                        _mm256_castsi256_ps(bisy));                            \
                    if (bkm && lane_counts[lane] + 8 <= max_per_lane) {        \
                        int skm = (~bkm) & 0xFF;                               \
                        lane_counts[lane] = csyncmer_append_filtered(          \
                            ap, skm,                                           \
                            ts_pos_bufs_##FUNC_NAME[lane],                     \
                            lane_counts[lane]);                                \
                    }                                                          \
                }                                                              \
            }                                                                  \
        } else if (COLLECT_POSITIONS && CANONICAL) {                           \
            /* Canonical positions: per-iteration scatter.                      \
             * Avoids old code's 64 scalar strand bit-extractions. */          \
            __m256i fs = _mm256_set1_epi32((int)kmer_idx);                    \
            __m256i ls = _mm256_set1_epi32(                                   \
                (int)(kmer_idx + window_size - 1));                            \
            __m256i isy = _mm256_or_si256(                                    \
                _mm256_cmpeq_epi32(min_pos_vec, fs),                          \
                _mm256_cmpeq_epi32(min_pos_vec, ls));                         \
            if (kmer_idx >= last_lane_limit) {                                 \
                __m256i lims = _mm256_sub_epi32(                              \
                    _mm256_load_si256((const __m256i*)chunk_ends),             \
                    _mm256_load_si256((const __m256i*)chunk_starts));          \
                isy = _mm256_and_si256(isy,                                   \
                    _mm256_cmpgt_epi32(lims, fs));                            \
            }                                                                  \
            int km = _mm256_movemask_ps(_mm256_castsi256_ps(isy));            \
            while (km) {                                                       \
                int lane = __builtin_ctz(km);                                  \
                size_t lc = lane_counts[lane];                                 \
                ts_pos_bufs_##FUNC_NAME[lane][lc] =                            \
                    (uint32_t)(chunk_starts[lane] + kmer_idx);                 \
                ts_str_bufs_##FUNC_NAME[lane][lc] =                            \
                    (overall_strand >> lane) & 1;                              \
                lane_counts[lane] = lc + 1;                                    \
                km &= km - 1;                                                  \
            }                                                                  \
        } else {                                                               \
            /* Count-only: syncmer check + popcount */                         \
            __m256i fs = _mm256_set1_epi32((int)kmer_idx);                    \
            __m256i ls = _mm256_set1_epi32(                                   \
                (int)(kmer_idx + window_size - 1));                            \
            __m256i isy = _mm256_or_si256(                                    \
                _mm256_cmpeq_epi32(min_pos_vec, fs),                          \
                _mm256_cmpeq_epi32(min_pos_vec, ls));                         \
            if (kmer_idx >= last_lane_limit) {                                 \
                __m256i lims = _mm256_sub_epi32(                              \
                    _mm256_load_si256((const __m256i*)chunk_ends),             \
                    _mm256_load_si256((const __m256i*)chunk_starts));          \
                isy = _mm256_and_si256(isy,                                   \
                    _mm256_cmpgt_epi32(lims, fs));                            \
            }                                                                  \
            syncmer_count += __builtin_popcount(                               \
                _mm256_movemask_ps(_mm256_castsi256_ps(isy)));                \
        }                                                                      \
    } /* end main loop */                                                      \
                                                                               \
    /* Process remaining partial batch / concatenate (positions only) */        \
    if (COLLECT_POSITIONS && !CANONICAL) {                                     \
        size_t partial = batch_count % 8;                                      \
        if (partial > 0) {                                                     \
            for (size_t pi = partial; pi < 8; pi++)                            \
                batch_pos[pi] = _mm256_setzero_si256();                        \
            __m256i tp[8];                                                     \
            csyncmer_transpose_8x8(batch_pos, tp);                            \
            __m256i bb = _mm256_set1_epi32((int)batch_base_kmer);             \
            for (int lane = 0; lane < 8; lane++) {                             \
                alignas(32) uint32_t absp[8];                                  \
                for (int j = 0; j < 8; j++)                                    \
                    absp[j] = (uint32_t)(chunk_starts[lane] +                  \
                                         batch_base_kmer + j);                \
                __m256i ap = _mm256_load_si256((__m256i*)absp);                \
                __m256i bfs = ap;                                              \
                __m256i bls = _mm256_add_epi32(ap,                             \
                    _mm256_set1_epi32((int)(window_size - 1)));                \
                __m256i lo = _mm256_set1_epi32((int)chunk_starts[lane]);       \
                __m256i amp = _mm256_add_epi32(tp[lane], lo);                  \
                __m256i bisy = _mm256_or_si256(                                \
                    _mm256_cmpeq_epi32(amp, bfs),                              \
                    _mm256_cmpeq_epi32(amp, bls));                             \
                alignas(32) uint32_t va[8];                                    \
                for (int j = 0; j < 8; j++) {                                 \
                    size_t ak = chunk_starts[lane] + batch_base_kmer + j;     \
                    va[j] = (j < (int)partial && ak < chunk_ends[lane])       \
                        ? 0xFFFFFFFF : 0;                                      \
                }                                                              \
                bisy = _mm256_and_si256(bisy,                                  \
                    _mm256_load_si256((__m256i*)va));                           \
                int bkm = _mm256_movemask_ps(                                  \
                    _mm256_castsi256_ps(bisy));                                \
                if (bkm && lane_counts[lane] + 8 <= max_per_lane) {            \
                    int skm = (~bkm) & 0xFF;                                   \
                    lane_counts[lane] = csyncmer_append_filtered(              \
                        ap, skm,                                               \
                        ts_pos_bufs_##FUNC_NAME[lane],                         \
                        lane_counts[lane]);                                    \
                }                                                              \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    if (COLLECT_POSITIONS) {                                                   \
        /* Concatenate lane buffers */                                         \
        size_t total = 0;                                                      \
        for (int i = 0; i < 8; i++) {                                          \
            if (lane_counts[i] > 0 &&                                          \
                total + lane_counts[i] <= max_positions) {                     \
                memcpy(out_positions + total,                                   \
                       ts_pos_bufs_##FUNC_NAME[i],                             \
                       lane_counts[i] * sizeof(uint32_t));                     \
                if (CANONICAL && out_strands) {                                \
                    memcpy(out_strands + total,                                 \
                           ts_str_bufs_##FUNC_NAME[i],                         \
                           lane_counts[i]);                                    \
                }                                                              \
                total += lane_counts[i];                                       \
            }                                                                  \
        }                                                                      \
        syncmer_count = total;                                                 \
    }                                                                          \
                                                                               \
    free(delay_buf);                                                           \
    free(packed);                                                              \
    return syncmer_count;                                                      \
}

/* Generate all 4 TWOSTACK variants */
CSYNCMER_DEFINE_TWOSTACK_SIMD_32(csyncmer_twostack_simd_32_count_,              0, 0)
CSYNCMER_DEFINE_TWOSTACK_SIMD_32(csyncmer_twostack_simd_32_positions_,          0, 1)
CSYNCMER_DEFINE_TWOSTACK_SIMD_32(csyncmer_twostack_simd_32_canonical_count_,    1, 0)
CSYNCMER_DEFINE_TWOSTACK_SIMD_32(csyncmer_twostack_simd_32_canonical_positions_,1, 1)

/* Public API wrappers preserving existing signatures */

static inline size_t csyncmer_twostack_simd_32_count(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    return csyncmer_twostack_simd_32_count_(sequence, length, K, S,
                                             NULL, NULL, 0);
}

static inline size_t csyncmer_twostack_simd_32_positions(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    uint32_t* out_positions,
    size_t max_positions
) {
    return csyncmer_twostack_simd_32_positions_(sequence, length, K, S,
                                                 out_positions, NULL,
                                                 max_positions);
}

static inline size_t csyncmer_twostack_simd_32_canonical_count(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S
) {
    return csyncmer_twostack_simd_32_canonical_count_(sequence, length, K, S,
                                                       NULL, NULL, 0);
}

static inline size_t csyncmer_twostack_simd_32_canonical_positions(
    const char* sequence,
    size_t length,
    size_t K,
    size_t S,
    uint32_t* out_positions,
    uint8_t* out_strands,
    size_t max_positions
) {
    return csyncmer_twostack_simd_32_canonical_positions_(
        sequence, length, K, S, out_positions, out_strands, max_positions);
}

// Note on TG-count optimization (from simd-minimizers)
// ============================================================================
//
// simd-minimizers uses TG-count for strand SELECTION in minimizers:
// - "Canonical" in minimizers means: pick ONE strand consistently based on TG-count
// - TG-count > 0 means more TG bases, so pick RC strand (which has more AC)
//
// For canonical SYNCMERS, we use a different definition:
// - "Canonical" means: use min(forward_hash, rc_hash) for each s-mer
// - We need the ACTUAL minimum hash value to find the minimum s-mer position
// - TG-count can't skip hash computation because hash values are position-dependent
//
// The TG-count trick doesn't apply to canonical syncmer hash computation.

#endif  // __AVX2__

#ifdef __cplusplus
}
#endif

#endif // CSYNCMER_FAST_H
