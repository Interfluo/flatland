#include "fl_field.hpp"

#include "fl_locale.hpp"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <cmath>
#include <cstddef>
#include <cstdlib>

namespace flatland {

namespace {

bool file_readable(const std::string& path) {
    if (path.empty()) return false;
    std::ifstream f(path, std::ios::binary);
    return f.good();
}

void strip_bom(std::string& s) {
    if (s.size() >= 3 && (unsigned char)s[0] == 0xEF
                      && (unsigned char)s[1] == 0xBB
                      && (unsigned char)s[2] == 0xBF)
        s.erase(0, 3);
}

bool blank(const std::string& s) {
    return s.find_first_not_of(" \t\r\n\f\v") == std::string::npos;
}

} // namespace

// Split "path@3" into {"path", 3}; a bare path means column 0.
//
// The filesystem gets the deciding vote: a file genuinely named "run@1" is a
// real thing, and splitting it unconditionally made such a file unreachable.
// Only when the token as written does not name a readable file is a trailing
// @<digits> treated as a column selector.
FieldToken parse_field_token(const std::string& token) {
    if (file_readable(token)) return { token, 0 };

    const size_t at = token.find_last_of('@');
    if (at != std::string::npos && at + 1 < token.size()) {
        const std::string num = token.substr(at + 1);
        if (num.find_first_not_of("0123456789") == std::string::npos) {
            unsigned long long col = 0;
            try {
                col = std::stoull(num);
            } catch (const std::out_of_range&) {
                // stoul's own what() is the bare word "stoul", which tells a user
                // nothing about which argument was wrong.
                throw std::runtime_error("field column '" + num + "' is too large");
            }
            if (col > (unsigned long long)std::numeric_limits<size_t>::max())
                throw std::runtime_error("field column '" + num + "' is too large");
            return { token.substr(0, at), (size_t)col };
        }
    }
    return { token, 0 };
}

// Load a field matrix and decide node-vs-face from the row count (or `forced`).
//
// Every non-blank line must yield numbers, and must be fully consumed. Skipping
// an unparseable row silently — a BOM on line 1, an "N/A" placeholder, a text
// header — shifts the remaining rows up and can flip the node/face inference
// entirely, so a per-vertex field is read as a per-face one and the run returns
// a plausible wrong answer with a zero exit status.
template <typename T>
void load_matrix_into(const std::string& filename, const Mesh<T>& mesh,
                      FieldMatrix<T>& m, ValueMode forced) {
    CNumericScope c_numeric;          // field files are written with '.'
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open data file '" + filename + "'");
    m.data.clear(); m.nrows = 0; m.ncols = 0;

    std::string line;
    size_t lineno = 0;
    size_t first_row_line = 0;

    auto fail = [&](const std::string& what) {
        throw std::runtime_error("data file '" + filename + "' line " +
                                 std::to_string(lineno) + ": " + what);
    };

    while (std::getline(file, line)) {
        ++lineno;
        if (lineno == 1) strip_bom(line);

        const size_t hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        if (blank(line)) continue;

        const char* p = line.c_str();
        char* end = nullptr;
        const size_t before = m.data.size();
        for (;;) {
            const double d = std::strtod(p, &end);
            if (end == p) break;
            if (!std::isfinite(d))
                fail("field value '" + std::string(p, (const char*)end) +
                     "' is not a finite number");
            m.data.push_back((T)d);
            p = end;
        }
        const size_t got = m.data.size() - before;

        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;
        if (got == 0)
            fail("expected numeric field values, found '" + line.substr(0, 40) + "'");
        if (*p != '\0')
            fail("unexpected text after the field values: '" + std::string(p).substr(0, 40) + "'");

        if (m.ncols == 0) { m.ncols = got; first_row_line = lineno; }
        else if (got != m.ncols)
            fail("has " + std::to_string(got) + " columns, but line " +
                 std::to_string(first_row_line) + " has " + std::to_string(m.ncols));
        ++m.nrows;
    }
    if (m.nrows == 0) throw std::runtime_error("data file '" + filename + "' has no values");

    if (forced != MODE_NONE) {
        size_t need = (forced == MODE_NODE) ? mesh.vertices.size() : mesh.faces.size();
        if (m.nrows != need)
            throw std::runtime_error("data file '" + filename + "' has " + std::to_string(m.nrows) +
                " rows but mesh has " + std::to_string(need) +
                (forced == MODE_NODE ? " vertices (--field-mode node)" : " faces (--field-mode face)"));
        m.mode = forced;
    } else if (m.nrows == mesh.vertices.size()) {
        m.mode = MODE_NODE;
    } else if (m.nrows == mesh.faces.size()) {
        m.mode = MODE_FACE;
    } else {
        throw std::runtime_error("data file '" + filename + "' has " + std::to_string(m.nrows) +
            " rows, but mesh has " + std::to_string(mesh.vertices.size()) + " vertices / " +
            std::to_string(mesh.faces.size()) + " faces; use --field-mode to disambiguate");
    }
}

template <typename T>
void extract_column(const FieldMatrix<T>& m, size_t col, Field<T>& out) {
    if (col >= m.ncols)
        throw std::runtime_error("field column " + std::to_string(col) +
            " out of range (file has " + std::to_string(m.ncols) + " column(s))");
    out.values.resize(m.nrows);
    for (size_t r = 0; r < m.nrows; ++r) out.values[r] = m.data[r * m.ncols + col];
    out.mode = m.mode;
}

template struct Field<float>;
template struct Field<double>;
template struct FieldMatrix<float>;
template struct FieldMatrix<double>;

template void load_matrix_into<float>(const std::string&, const Mesh<float>&, FieldMatrix<float>&, ValueMode);
template void load_matrix_into<double>(const std::string&, const Mesh<double>&, FieldMatrix<double>&, ValueMode);
template void extract_column<float>(const FieldMatrix<float>&, size_t, Field<float>&);
template void extract_column<double>(const FieldMatrix<double>&, size_t, Field<double>&);

} // namespace flatland
