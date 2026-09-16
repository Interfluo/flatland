// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <cstddef>

namespace flatland {

/*
 * A small DEFLATE (RFC 1951) compressor, written here rather than linked,
 * because FlatLand's premise is that one compiler invocation and nothing else
 * builds it — and PNG is the only reason we need compression at all.
 *
 * It implements LZ77 matching with a hash chain, emitted through DEFLATE's
 * FIXED Huffman tables. Fixed tables mean no code-length tree has to be built
 * or transmitted, which costs a few percent of ratio against a full dynamic
 * encoder and removes most of the complexity and all of the ways to get it
 * subtly wrong. On FlatLand's own output that lands within a few percent of
 * zlib's default level and about 13x smaller than the uncompressed raster.
 *
 * The output is a standard zlib stream (RFC 1950): a two-byte header, the
 * DEFLATE data, and an Adler-32 checksum. Any conforming inflater reads it.
 */

// Compress `data` into a zlib stream. `level` selects how hard the matcher
// looks: 0 stores without matching, 1 is a short chain, 9 a long one.
std::vector<uint8_t> zlib_compress(const uint8_t* data, size_t size, int level = 6);

// Checksums, exposed because PNG needs CRC-32 for its chunks.
uint32_t crc32_bytes(const uint8_t* data, size_t size, uint32_t seed = 0);
uint32_t adler32_bytes(const uint8_t* data, size_t size);

} // namespace flatland
