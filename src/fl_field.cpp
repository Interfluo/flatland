#include "fl_field.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <cstddef>
#include <cstdlib>

namespace flatland {

// Split "path@3" into {"path", 3}; a bare path means column 0.
FieldToken parse_field_token(const std::string& token) {
    size_t at = token.find_last_of('@');
    if (at != std::string::npos && at + 1 < token.size()) {
        const std::string num = token.substr(at + 1);
        if (num.find_first_not_of("0123456789") == std::string::npos)
            return { token.substr(0, at), (size_t)std::stoul(num) };
    }
    return { token, 0 };
}

// Load a field matrix and decide node-vs-face from the row count (or `forced`).
template <typename T>
void load_matrix_into(const std::string& filename, const Mesh<T>& mesh,
                      FieldMatrix<T>& m, ValueMode forced) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) throw std::runtime_error("cannot open data file '" + filename + "'");
    m.data.clear(); m.nrows = 0; m.ncols = 0;

    std::string line;
    while (std::getline(file, line)) {
        size_t c = line.find('#');
        if (c != std::string::npos) line = line.substr(0, c);
        const char* p = line.c_str(); char* end;
        size_t before = m.data.size();
        double d = std::strtod(p, &end);
        while (end != p) { m.data.push_back((T)d); p = end; d = std::strtod(p, &end); }
        size_t got = m.data.size() - before;
        if (got == 0) continue;                 // blank / comment-only line
        if (m.ncols == 0) m.ncols = got;
        else if (got != m.ncols)
            throw std::runtime_error("data file '" + filename + "' row " + std::to_string(m.nrows+1) +
                " has " + std::to_string(got) + " columns, expected " + std::to_string(m.ncols));
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
