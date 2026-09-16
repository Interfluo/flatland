// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#include "fl_deflate.hpp"

#include <algorithm>
#include <cstring>

namespace flatland {

namespace {

/* ------------------------------------------------------------------------ */
/* Bit output                                                               */
/* ------------------------------------------------------------------------ */
/*
 * DEFLATE packs bits into bytes starting at the LEAST significant bit. Huffman
 * codes are the exception: their bits are emitted starting from the MOST
 * significant, so they need reversing on the way out. Everything else — the
 * block header, the extra bits that follow a length or distance code — is
 * written low bit first. Getting these two the wrong way round produces a
 * stream that looks plausible and inflates to garbage, so they are separate
 * functions rather than one with a flag.
 */
class BitWriter {
public:
    explicit BitWriter(std::vector<uint8_t>& out) : out_(out) {}

    void bits(uint32_t value, int count) {         // low bit first
        while (count-- > 0) {
            acc_ |= (value & 1u) << nbits_;
            value >>= 1;
            if (++nbits_ == 8) { out_.push_back((uint8_t)acc_); acc_ = 0; nbits_ = 0; }
        }
    }

    void code(uint32_t huff, int len) {            // high bit first
        while (len-- > 0) bits((huff >> len) & 1u, 1);
    }

    void align() {                                  // pad to a byte boundary
        if (nbits_) { out_.push_back((uint8_t)acc_); acc_ = 0; nbits_ = 0; }
    }

private:
    std::vector<uint8_t>& out_;
    uint32_t acc_ = 0;
    int nbits_ = 0;
};

/* ------------------------------------------------------------------------ */
/* Fixed Huffman tables (RFC 1951, 3.2.6)                                    */
/* ------------------------------------------------------------------------ */

// Literal/length symbol -> (code, bit length).
inline void fixed_litlen(unsigned sym, uint32_t& code, int& len) {
    if (sym < 144)      { code = 0x30  + sym;         len = 8; }
    else if (sym < 256) { code = 0x190 + (sym - 144); len = 9; }
    else if (sym < 280) { code = 0x00  + (sym - 256); len = 7; }
    else                { code = 0xC0  + (sym - 280); len = 8; }
}

// Length 3..258 -> symbol 257..285, plus its extra bits.
const uint16_t LEN_BASE[29] = {
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31,
    35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
};
const uint8_t LEN_EXTRA[29] = {
    0,0,0,0,0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, 0
};

// Distance 1..32768 -> symbol 0..29, plus its extra bits.
const uint16_t DIST_BASE[30] = {
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193,
    257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145,
    8193, 12289, 16385, 24577
};
const uint8_t DIST_EXTRA[30] = {
    0,0,0,0, 1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,8, 9,9,
    10,10, 11,11, 12,12, 13,13
};

int length_symbol(int len) {
    for (int i = 28; i >= 0; --i) if (len >= LEN_BASE[i]) return i;
    return 0;
}

int distance_symbol(int dist) {
    for (int i = 29; i >= 0; --i) if (dist >= DIST_BASE[i]) return i;
    return 0;
}

/* ------------------------------------------------------------------------ */
/* LZ77 match finding                                                        */
/* ------------------------------------------------------------------------ */

const int MIN_MATCH = 3;
const int MAX_MATCH = 258;
const int WINDOW    = 32768;          // DEFLATE's maximum distance
const int HASH_BITS = 15;
const int HASH_SIZE = 1 << HASH_BITS;

inline uint32_t hash3(const uint8_t* p) {
    return (uint32_t)(((uint32_t)p[0] << 10) ^ ((uint32_t)p[1] << 5) ^ p[2]) & (HASH_SIZE - 1);
}

struct Matcher {
    const uint8_t* data;
    size_t size;
    std::vector<int32_t> head;        // hash -> most recent position
    std::vector<int32_t> prev;        // position -> previous position, same hash
    int max_chain;

    Matcher(const uint8_t* d, size_t n, int chain)
        : data(d), size(n), head(HASH_SIZE, -1), prev(n, -1), max_chain(chain) {}

    void insert(size_t pos) {
        if (pos + MIN_MATCH > size) return;
        const uint32_t h = hash3(data + pos);
        prev[pos] = head[h];
        head[h] = (int32_t)pos;
    }

    // Longest match for the bytes at `pos`. Returns its length (0 if none) and
    // sets `dist`. Only looks backwards, never past the 32 KiB window.
    int find(size_t pos, int& dist) const {
        if (pos + MIN_MATCH > size) return 0;
        const size_t max_len = std::min<size_t>(MAX_MATCH, size - pos);
        if (max_len < MIN_MATCH) return 0;

        int best_len = 0, best_dist = 0;
        int chain = max_chain;
        int32_t cand = head[hash3(data + pos)];

        while (cand >= 0 && chain-- > 0) {
            const size_t d = pos - (size_t)cand;
            if (d == 0 || d > (size_t)WINDOW) break;

            // Cheap rejection before the full comparison: if the byte that
            // would extend the current best does not match, this candidate
            // cannot beat it.
            if (best_len > 0 && data[cand + best_len] != data[pos + best_len]) {
                cand = prev[cand];
                continue;
            }
            size_t len = 0;
            while (len < max_len && data[cand + len] == data[pos + len]) ++len;

            if ((int)len > best_len) {
                best_len = (int)len;
                best_dist = (int)d;
                if (len >= max_len) break;          // cannot do better
            }
            cand = prev[cand];
        }
        if (best_len >= MIN_MATCH) { dist = best_dist; return best_len; }
        return 0;
    }
};

void emit_literal(BitWriter& bw, uint8_t byte) {
    uint32_t code; int len;
    fixed_litlen(byte, code, len);
    bw.code(code, len);
}

void emit_match(BitWriter& bw, int length, int dist) {
    const int ls = length_symbol(length);
    uint32_t code; int len;
    fixed_litlen(257 + ls, code, len);
    bw.code(code, len);
    if (LEN_EXTRA[ls]) bw.bits((uint32_t)(length - LEN_BASE[ls]), LEN_EXTRA[ls]);

    const int ds = distance_symbol(dist);
    bw.code((uint32_t)ds, 5);                        // fixed distance codes are 5 bits
    if (DIST_EXTRA[ds]) bw.bits((uint32_t)(dist - DIST_BASE[ds]), DIST_EXTRA[ds]);
}

// A single stored (uncompressed) DEFLATE block chain, for level 0 and for the
// degenerate empty input.
void store_blocks(BitWriter& bw, const uint8_t* data, size_t size, std::vector<uint8_t>& out) {
    size_t pos = 0;
    do {
        const uint16_t chunk = (uint16_t)std::min<size_t>(65535, size - pos);
        const bool last = (pos + chunk >= size);
        bw.bits(last ? 1 : 0, 1);
        bw.bits(0, 2);                               // BTYPE 00, stored
        bw.align();
        out.push_back((uint8_t)(chunk & 0xFF));
        out.push_back((uint8_t)(chunk >> 8));
        out.push_back((uint8_t)(~chunk & 0xFF));
        out.push_back((uint8_t)((~chunk >> 8) & 0xFF));
        out.insert(out.end(), data + pos, data + pos + chunk);
        pos += chunk;
    } while (pos < size);
}

/* The CRC-32 table (IEEE 802.3 polynomial, reflected: 0xEDB88320).
 *
 * Built at compile time. It began life as a function-local `static uint32_t
 * table[256]` behind a `static bool ready` flag, which is a data race: every
 * worker writing a PNG raced to fill the same table and set the same flag.
 * Each thread wrote identical bytes, so it never produced a wrong checksum in
 * practice, but it is undefined behaviour and ThreadSanitizer rightly failed
 * on it -- intermittently, which is worse than failing every time.
 *
 * constexpr removes the problem rather than synchronising it: there is no
 * initialisation left to race over, and no atomic to check on each call. */
struct Crc32Table { uint32_t v[256]; };

constexpr Crc32Table make_crc32_table() {
    Crc32Table t{};
    for (uint32_t i = 0; i < 256; ++i) {
        uint32_t c = i;
        for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
        t.v[i] = c;
    }
    return t;
}

constexpr Crc32Table CRC32_TABLE = make_crc32_table();

} // namespace

/* ------------------------------------------------------------------------ */
/* Checksums                                                                 */
/* ------------------------------------------------------------------------ */

uint32_t crc32_bytes(const uint8_t* data, size_t size, uint32_t seed) {
    uint32_t c = seed ^ 0xFFFFFFFFu;
    for (size_t i = 0; i < size; ++i) c = CRC32_TABLE.v[(c ^ data[i]) & 0xFF] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

uint32_t adler32_bytes(const uint8_t* data, size_t size) {
    uint32_t a = 1, b = 0;
    const uint32_t MOD = 65521;
    // Chunked so the accumulators cannot overflow before the modulo.
    while (size) {
        const size_t n = std::min<size_t>(size, 5552);
        for (size_t i = 0; i < n; ++i) { a += data[i]; b += a; }
        a %= MOD; b %= MOD;
        data += n; size -= n;
    }
    return (b << 16) | a;
}

/* ------------------------------------------------------------------------ */
/* Entry point                                                               */
/* ------------------------------------------------------------------------ */

std::vector<uint8_t> zlib_compress(const uint8_t* data, size_t size, int level) {
    std::vector<uint8_t> out;
    out.reserve(size / 3 + 64);

    // zlib header (RFC 1950): deflate, 32 KiB window, no preset dictionary.
    // The check bits make the 16-bit header a multiple of 31.
    // These are the canonical FLG bytes for CMF 0x78; each already makes
    // (CMF << 8 | FLG) a multiple of 31, which is what the check requires.
    // FLEVEL is advisory only and nothing reads it.
    const uint8_t cmf = 0x78;
    const uint8_t flg = (level <= 1) ? 0x01 : (level >= 9 ? 0xDA : 0x9C);
    out.push_back(cmf);
    out.push_back(flg);

    {
        BitWriter bw(out);
        if (size == 0 || level <= 0) {
            if (size == 0) {
                bw.bits(1, 1); bw.bits(0, 2); bw.align();
                out.push_back(0); out.push_back(0); out.push_back(0xFF); out.push_back(0xFF);
            } else {
                store_blocks(bw, data, size, out);
            }
        } else {
            const size_t body_start = out.size();
            const int chain = (level <= 3) ? 16 : (level <= 6 ? 64 : 256);
            Matcher m(data, size, chain);

            bw.bits(1, 1);                           // final block
            bw.bits(1, 2);                           // BTYPE 01, fixed Huffman

            size_t pos = 0;
            while (pos < size) {
                int dist = 0;
                const int len = m.find(pos, dist);
                if (len >= MIN_MATCH) {
                    emit_match(bw, len, dist);
                    // Every position inside the match still has to enter the
                    // hash chain, or later matches lose candidates.
                    for (int k = 0; k < len; ++k) m.insert(pos + k);
                    pos += len;
                } else {
                    emit_literal(bw, data[pos]);
                    m.insert(pos);
                    ++pos;
                }
            }
            uint32_t code; int clen;
            fixed_litlen(256, code, clen);           // end-of-block
            bw.code(code, clen);
            bw.align();

            // Incompressible input (noise, already-compressed data) makes fixed
            // Huffman coding slightly EXPAND it — 8-bit codes for bytes that
            // carry 8 bits of entropy, plus block overhead. Falling back to
            // stored blocks means the output is never larger than the input by
            // more than the framing, which matters because a raster of pure
            // noise is a perfectly legitimate thing to render.
            const size_t stored_overhead = 5 * (size / 65535 + 1);
            if (out.size() - body_start > size + stored_overhead) {
                out.resize(body_start);
                BitWriter sw(out);
                store_blocks(sw, data, size, out);
            }
        }
    }

    const uint32_t adler = adler32_bytes(data, size);
    out.push_back((uint8_t)(adler >> 24));
    out.push_back((uint8_t)(adler >> 16));
    out.push_back((uint8_t)(adler >> 8));
    out.push_back((uint8_t)(adler));
    return out;
}

} // namespace flatland
