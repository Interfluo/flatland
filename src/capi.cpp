/*
 * capi.cpp — implementation of the public C ABI declared in include/flatland.h.
 *
 * Everything here is a thin, defensive layer over the templated engine. Its
 * jobs are: validate what crosses the boundary, translate C++ exceptions into
 * status codes plus a thread-local message, and keep the caller's arrays and
 * the engine's storage clearly separated. No exception may escape into C.
 */

/*
 * Which side of the Windows import/export boundary this translation unit is on
 * is a property of the LINK, not of the source, so the build system says it:
 *   shared library  ->  -DFLATLAND_BUILD_SHARED  (FL_API = __declspec(dllexport))
 *   static library  ->  -DFLATLAND_STATIC        (FL_API = empty)
 * Defining it here unconditionally, as this file used to, stamped dllexport on
 * the STATIC archive too, which makes every consumer of the .lib re-export the
 * C ABI from its own binary. Elsewhere (ELF/Mach-O) neither macro has any
 * effect: FL_API is visibility("default") either way.
 *
 * The fallback keeps a hand-rolled compile working: with neither macro set the
 * header would resolve FL_API to dllimport, and a definition marked dllimport
 * is ill-formed, so default to the export side rather than fail obscurely.
 */
#if !defined(FLATLAND_BUILD_SHARED) && !defined(FLATLAND_STATIC)
#  define FLATLAND_BUILD_SHARED 1
#endif
#include "flatland.h"

#include "fl_batch.hpp"
#include "fl_field.hpp"
#include "fl_image.hpp"
#include "fl_io.hpp"
#include "fl_mesh.hpp"
#include "fl_parallel.hpp"
#include "fl_project.hpp"
#include "fl_raster.hpp"
#include "fl_vec.hpp"

#include <algorithm>
#include <memory>
#include <mutex>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>

using namespace flatland;

/* ------------------------------------------------------------------------- */
/* Error reporting                                                            */
/* ------------------------------------------------------------------------- */

namespace {

// Thread-local so concurrent calls cannot overwrite each other's diagnostics.
thread_local std::string g_error;

fl_status ok() { g_error.clear(); return FL_OK; }

fl_status fail(fl_status s, const std::string& msg) {
    g_error = msg;
    return s;
}

// Run `f`, mapping any escaping exception onto a status code. `on_error` is the
// category appropriate to what `f` was doing — the call site knows whether it
// was reading a file or projecting geometry, which the exception itself does not
// record, so classification happens here rather than by inspecting messages.
template <typename F>
fl_status guard(fl_status on_error, F&& f) {
    try {
        f();
        return ok();
    } catch (const std::bad_alloc&) {
        return fail(FL_ERR_OUT_OF_MEMORY, "out of memory");
    } catch (const std::exception& e) {
        return fail(on_error, e.what());
    } catch (...) {
        return fail(FL_ERR_INTERNAL, "unknown error");
    }
}

} // namespace

/* ------------------------------------------------------------------------- */
/* Handles                                                                    */
/* ------------------------------------------------------------------------- */

// A mesh is held canonically in double. The float copy is built on first use
// rather than eagerly, so a caller who only ever runs double precision does not
// pay for it; std::once_flag keeps that lazy build safe for the concurrent
// sharing the header promises.
struct fl_mesh {
    Mesh<double> md;
    std::once_flag f_once;
    std::unique_ptr<Mesh<float>> mf;

    const Mesh<float>& as_float() {
        std::call_once(f_once, [&]() {
            mf.reset(new Mesh<float>(narrow_mesh<float, double>(md)));
        });
        return *mf;
    }
};

struct fl_image {
    int32_t w = 0, h = 0;
    std::vector<uint8_t> mask;
    std::vector<double> values;      // empty when the view carried no field
    double min_v = 0, max_v = 0;
    bool has_values = false;
};

struct fl_field_matrix {
    FieldMatrix<double> m;
};

/* ------------------------------------------------------------------------- */
/* Shared helpers                                                             */
/* ------------------------------------------------------------------------- */

namespace {

const fl_options& effective(const fl_options* opts, fl_options& storage) {
    if (opts) return *opts;
    fl_options_init(&storage);
    return storage;
}

// Reading an enum-typed field whose stored value is outside the enumeration is
// undefined behaviour, and validating exactly that is this function's job: a
// foreign caller can put any bit pattern in the struct. An optimizer permitted
// to assume fl_precision only ever holds 0 or 1 could delete the check. So the
// raw bytes are copied out as int32 and compared as integers.
int32_t enum_bits(const void* p) {
    int32_t v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

fl_status check_options(const fl_options& o) {
    if (!std::isfinite(o.resolution) || o.resolution <= 0)
        return fail(FL_ERR_INVALID_ARGUMENT, "options.resolution must be a positive, finite number");
    const int32_t prec = enum_bits(&o.precision);
    const int32_t fmode = enum_bits(&o.field_mode);
    if (prec != FL_PRECISION_FLOAT && prec != FL_PRECISION_DOUBLE)
        return fail(FL_ERR_INVALID_ARGUMENT, "options.precision is not a valid fl_precision value");
    if (fmode != FL_FIELD_AUTO && fmode != FL_FIELD_NODE && fmode != FL_FIELD_FACE)
        return fail(FL_ERR_INVALID_ARGUMENT, "options.field_mode is not a valid fl_field_mode value");
    if (o.threads < 0)
        return fail(FL_ERR_INVALID_ARGUMENT, "options.threads cannot be negative (0 means one per core)");
    return FL_OK;
}

// Decide whether a field of `len` values is per-vertex or per-face.
fl_status resolve_field_mode(const fl_mesh* mesh, size_t len, fl_field_mode requested,
                             ValueMode& out) {
    const size_t nv = mesh->md.vertices.size();
    const size_t nf = mesh->md.faces.size();
    switch (requested) {
        case FL_FIELD_NODE:
            if (len != nv)
                return fail(FL_ERR_DIMENSION, "field has " + std::to_string(len) +
                            " values but the mesh has " + std::to_string(nv) + " vertices");
            out = MODE_NODE; return FL_OK;
        case FL_FIELD_FACE:
            if (len != nf)
                return fail(FL_ERR_DIMENSION, "field has " + std::to_string(len) +
                            " values but the mesh has " + std::to_string(nf) + " faces");
            out = MODE_FACE; return FL_OK;
        default: break;
    }
    if (nv == nf && len == nv)
        return fail(FL_ERR_DIMENSION,
                    "this mesh has equally many vertices and faces (" + std::to_string(nv) +
                    "), so a field of that length is ambiguous; set options.field_mode explicitly");
    if (len == nv) { out = MODE_NODE; return FL_OK; }
    if (len == nf) { out = MODE_FACE; return FL_OK; }
    return fail(FL_ERR_DIMENSION, "field has " + std::to_string(len) +
                " values, but the mesh has " + std::to_string(nv) + " vertices / " +
                std::to_string(nf) + " faces");
}

template <typename T>
void fill_field(Field<T>& dst, const double* src, size_t len, ValueMode mode) {
    dst.values.resize(len);
    for (size_t i = 0; i < len; ++i) dst.values[i] = (T)src[i];
    dst.mode = mode;
}

template <typename T>
void to_result(const ViewResult<T>& r, fl_result& out) {
    std::memset(&out, 0, sizeof(out));
    out.area           = (double)r.area;
    out.covered_pixels = (int64_t)r.covered_pixels;
    out.width          = r.image_width;
    out.height         = r.image_height;
    out.has_field      = r.has_field ? 1 : 0;
    out.has_stats      = r.has_stats ? 1 : 0;
    if (r.has_stats) {
        out.average  = (double)r.average_value;
        out.integral = (double)r.integral;
        out.min      = (double)r.min_val;
        out.max      = (double)r.max_val;
    }
}

bool dir_ok(const double d[3]) {
    return std::isfinite(d[0]) && std::isfinite(d[1]) && std::isfinite(d[2]) &&
           !(d[0] == 0 && d[1] == 0 && d[2] == 0);
}

template <typename T>
const Mesh<T>& typed_mesh(fl_mesh* m);
template <> const Mesh<double>& typed_mesh<double>(fl_mesh* m) { return m->md; }
template <> const Mesh<float>&  typed_mesh<float>(fl_mesh* m)  { return m->as_float(); }

// One view, in the requested working precision.
template <typename T>
void project_one(fl_mesh* mesh, const double view[3], const double* field, size_t field_len,
                 ValueMode mode, const fl_options& o, fl_result* out, Renderer<T>& renderer) {
    Field<T> f;
    if (field) fill_field(f, field, field_len, mode);
    const ViewResult<T> r = process_view<T>({(T)view[0], (T)view[1], (T)view[2]},
                                            typed_mesh<T>(mesh), field ? &f : nullptr,
                                            (T)o.resolution, o.cull != 0, renderer);
    to_result(r, *out);
}

template <typename T>
void batch_impl(fl_mesh* mesh, const fl_batch_desc& d, const fl_options& o,
                ValueMode mode, fl_result* out) {
    const unsigned workers = choose_workers(d.view_count, (unsigned)o.threads);
    std::vector<Renderer<T>> renderers(workers);
    std::vector<Field<T>> fields(workers);

    const Mesh<T>& m = typed_mesh<T>(mesh);

    try {
        parallel_for(d.view_count, workers, [&](size_t k, unsigned w) {
            const double* vd = d.views + k*3;
            if (!dir_ok(vd))
                throw std::runtime_error("view direction must be a non-zero, finite vector");

            const Field<T>* fp = nullptr;
            if (d.field_matrix) {
                const size_t col = d.field_columns ? (size_t)d.field_columns[k] : 0;
                if (col >= d.field_cols)
                    throw std::runtime_error("field column " + std::to_string(col) +
                                             " is out of range (" + std::to_string(d.field_cols) +
                                             " column(s))");
                Field<T>& f = fields[w];
                f.values.resize(d.field_rows);
                for (size_t r = 0; r < d.field_rows; ++r)
                    f.values[r] = (T)d.field_matrix[r * d.field_cols + col];
                f.mode = mode;
                fp = &f;
            }

            double res = o.resolution;
            if (d.resolutions) {
                res = d.resolutions[k];
                if (!std::isfinite(res) || res <= 0)
                    throw std::runtime_error("per-view resolution must be a positive, finite number");
            }

            const ViewResult<T> r = process_view<T>({(T)vd[0], (T)vd[1], (T)vd[2]}, m, fp,
                                                    (T)res, o.cull != 0, renderers[w]);
            to_result(r, out[k]);
        });
    } catch (const ParallelError& e) {
        throw std::runtime_error("view " + std::to_string(e.index) + ": " + e.what());
    }
}

template <typename T>
void image_impl(fl_mesh* mesh, const double view[3], const double* field, size_t field_len,
                ValueMode mode, const fl_options& o, fl_image* img, fl_result* out) {
    Renderer<T> renderer;
    Field<T> f;
    if (field) fill_field(f, field, field_len, mode);
    const ViewResult<T> r = process_view<T>({(T)view[0], (T)view[1], (T)view[2]},
                                            typed_mesh<T>(mesh), field ? &f : nullptr,
                                            (T)o.resolution, o.cull != 0, renderer);
    if (out) to_result(r, *out);

    img->w = r.image_width;
    img->h = r.image_height;
    img->has_values = r.has_field && r.image_width > 0;
    img->min_v = (double)r.min_val;
    img->max_v = (double)r.max_val;

    const size_t n = (size_t)std::max(0, img->w) * (size_t)std::max(0, img->h);
    img->mask.assign(renderer.mask.begin(), renderer.mask.begin() + n);
    if (img->has_values) {
        img->values.resize(n);
        for (size_t i = 0; i < n; ++i) img->values[i] = (double)renderer.val_buffer[i];
    } else {
        img->values.clear();
    }
}

template <typename Src>
fl_status make_mesh(const Src* vertices, size_t nv, const int32_t* faces, size_t nf,
                    fl_mesh** out_mesh) {
    if (!vertices || !faces || !out_mesh)
        return fail(FL_ERR_INVALID_ARGUMENT, "vertices, faces and out_mesh must all be non-NULL");
    if (nv == 0) return fail(FL_ERR_DIMENSION, "a mesh needs at least one vertex");
    if (nf == 0) return fail(FL_ERR_DIMENSION, "a mesh needs at least one face");
    *out_mesh = nullptr;

    return guard(FL_ERR_DIMENSION, [&]() {
        RawMesh raw;
        raw.vertices.resize(nv);
        for (size_t i = 0; i < nv; ++i)
            raw.vertices[i] = { (double)vertices[i*3+0],
                                (double)vertices[i*3+1],
                                (double)vertices[i*3+2] };
        raw.faces.resize(nf);
        for (size_t i = 0; i < nf; ++i)
            raw.faces[i] = { faces[i*3+0], faces[i*3+1], faces[i*3+2] };

        std::unique_ptr<fl_mesh> m(new fl_mesh());
        m->md = build_mesh<double>(std::move(raw), "<arrays>");
        *out_mesh = m.release();
    });
}

} // namespace

/* ------------------------------------------------------------------------- */
/* Version and errors                                                         */
/* ------------------------------------------------------------------------- */

extern "C" {

const char* fl_status_string(fl_status s) {
    switch (s) {
        case FL_OK:                   return "ok";
        case FL_ERR_INVALID_ARGUMENT: return "invalid argument";
        case FL_ERR_OUT_OF_MEMORY:    return "out of memory";
        case FL_ERR_IO:               return "I/O error";
        case FL_ERR_PARSE:            return "parse error";
        case FL_ERR_DIMENSION:        return "dimension mismatch";
        case FL_ERR_NUMERIC:          return "numeric error";
        case FL_ERR_UNSUPPORTED:      return "unsupported";
        case FL_ERR_INTERNAL:         return "internal error";
    }
    return "unrecognized status";
}

const char* fl_last_error(void) { return g_error.c_str(); }

void fl_version(int* major, int* minor, int* patch) {
    if (major) *major = FL_VERSION_MAJOR;
    if (minor) *minor = FL_VERSION_MINOR;
    if (patch) *patch = FL_VERSION_PATCH;
}

const char* fl_version_string(void) {
    static const char* v = "3.0.0";
    return v;
}

void fl_angle_to_dir(double az_deg, double el_deg, double out_dir[3]) {
    if (!out_dir) return;
    const Vec3<double> d = angle_to_dir(az_deg, el_deg);
    out_dir[0] = d.x; out_dir[1] = d.y; out_dir[2] = d.z;
}

/* ------------------------------------------------------------------------- */
/* Options and descriptors                                                    */
/* ------------------------------------------------------------------------- */

void fl_options_init(fl_options* o) {
    if (!o) return;
    std::memset(o, 0, sizeof(*o));
    o->resolution = 0.001;
    o->cull       = 1;
    o->precision  = FL_PRECISION_FLOAT;
    o->field_mode = FL_FIELD_AUTO;
    o->threads    = 0;
}

void fl_batch_desc_init(fl_batch_desc* d) {
    if (!d) return;
    std::memset(d, 0, sizeof(*d));
}

/* ------------------------------------------------------------------------- */
/* Mesh                                                                       */
/* ------------------------------------------------------------------------- */

fl_status fl_mesh_create(const double* vertices, size_t nv,
                         const int32_t* faces, size_t nf, fl_mesh** out_mesh) {
    return make_mesh<double>(vertices, nv, faces, nf, out_mesh);
}

fl_status fl_mesh_create_f32(const float* vertices, size_t nv,
                             const int32_t* faces, size_t nf, fl_mesh** out_mesh) {
    return make_mesh<float>(vertices, nv, faces, nf, out_mesh);
}

fl_status fl_mesh_load(const char* path, fl_mesh** out_mesh) {
    if (!path || !out_mesh)
        return fail(FL_ERR_INVALID_ARGUMENT, "path and out_mesh must both be non-NULL");
    *out_mesh = nullptr;
    return guard(FL_ERR_PARSE, [&]() {
        std::unique_ptr<fl_mesh> m(new fl_mesh());
        m->md = load_mesh<double>(path);
        *out_mesh = m.release();
    });
}

void   fl_mesh_destroy(fl_mesh* m)            { delete m; }
size_t fl_mesh_vertex_count(const fl_mesh* m) { return m ? m->md.vertices.size() : 0; }
size_t fl_mesh_face_count(const fl_mesh* m)   { return m ? m->md.faces.size() : 0; }

fl_status fl_mesh_copy_vertices(const fl_mesh* m, double* out, size_t out_len) {
    if (!m || !out) return fail(FL_ERR_INVALID_ARGUMENT, "mesh and out must both be non-NULL");
    const size_t need = m->md.vertices.size() * 3;
    if (out_len < need)
        return fail(FL_ERR_DIMENSION, "out needs room for " + std::to_string(need) +
                    " doubles, got " + std::to_string(out_len));
    // Undo the internal recentering so the caller gets back the frame it gave us.
    for (size_t i = 0; i < m->md.vertices.size(); ++i) {
        out[i*3+0] = m->md.vertices[i].x + m->md.origin.x;
        out[i*3+1] = m->md.vertices[i].y + m->md.origin.y;
        out[i*3+2] = m->md.vertices[i].z + m->md.origin.z;
    }
    return ok();
}

fl_status fl_mesh_copy_faces(const fl_mesh* m, int32_t* out, size_t out_len) {
    if (!m || !out) return fail(FL_ERR_INVALID_ARGUMENT, "mesh and out must both be non-NULL");
    const size_t need = m->md.faces.size() * 3;
    if (out_len < need)
        return fail(FL_ERR_DIMENSION, "out needs room for " + std::to_string(need) +
                    " int32 values, got " + std::to_string(out_len));
    for (size_t i = 0; i < m->md.faces.size(); ++i) {
        out[i*3+0] = m->md.faces[i].v0_idx;
        out[i*3+1] = m->md.faces[i].v1_idx;
        out[i*3+2] = m->md.faces[i].v2_idx;
    }
    return ok();
}

/* ------------------------------------------------------------------------- */
/* Projection                                                                 */
/* ------------------------------------------------------------------------- */

fl_status fl_project(const fl_mesh* mesh_c, const double view_dir[3],
                     const double* field, size_t field_len,
                     const fl_options* opts, fl_result* out_result) {
    if (!mesh_c || !view_dir || !out_result)
        return fail(FL_ERR_INVALID_ARGUMENT, "mesh, view_dir and out_result must all be non-NULL");
    if (field && field_len == 0)
        return fail(FL_ERR_INVALID_ARGUMENT, "field is non-NULL but field_len is 0");
    if (!dir_ok(view_dir))
        return fail(FL_ERR_INVALID_ARGUMENT, "view_dir must be a non-zero, finite vector");

    fl_options storage;
    const fl_options& o = effective(opts, storage);
    if (fl_status s = check_options(o)) return s;

    fl_mesh* mesh = const_cast<fl_mesh*>(mesh_c);   // only the lazy float copy mutates
    ValueMode mode = MODE_NONE;
    if (field) {
        if (fl_status s = resolve_field_mode(mesh, field_len, o.field_mode, mode)) return s;
    }

    return guard(FL_ERR_NUMERIC, [&]() {
        if (enum_bits(&o.precision) == FL_PRECISION_DOUBLE) {
            Renderer<double> r;
            project_one<double>(mesh, view_dir, field, field_len, mode, o, out_result, r);
        } else {
            Renderer<float> r;
            project_one<float>(mesh, view_dir, field, field_len, mode, o, out_result, r);
        }
    });
}

fl_status fl_project_batch(const fl_mesh* mesh_c, const fl_batch_desc* desc,
                           const fl_options* opts, fl_result* out_results) {
    if (!mesh_c || !desc || !out_results)
        return fail(FL_ERR_INVALID_ARGUMENT, "mesh, desc and out_results must all be non-NULL");
    if (!desc->views || desc->view_count == 0)
        return fail(FL_ERR_INVALID_ARGUMENT, "desc->views must be non-NULL with view_count > 0");
    if (desc->field_matrix && (desc->field_rows == 0 || desc->field_cols == 0))
        return fail(FL_ERR_INVALID_ARGUMENT,
                    "desc->field_matrix is set but field_rows/field_cols are zero");

    fl_options storage;
    const fl_options& o = effective(opts, storage);
    if (fl_status s = check_options(o)) return s;

    fl_mesh* mesh = const_cast<fl_mesh*>(mesh_c);
    ValueMode mode = MODE_NONE;
    if (desc->field_matrix) {
        if (fl_status s = resolve_field_mode(mesh, desc->field_rows, o.field_mode, mode)) return s;
        if (desc->field_columns) {
            for (size_t k = 0; k < desc->view_count; ++k) {
                const int32_t c = desc->field_columns[k];
                if (c < 0 || (size_t)c >= desc->field_cols)
                    return fail(FL_ERR_DIMENSION, "field_columns[" + std::to_string(k) +
                                "] = " + std::to_string(c) + " is out of range (" +
                                std::to_string(desc->field_cols) + " column(s))");
            }
        }
    }

    return guard(FL_ERR_NUMERIC, [&]() {
        if (enum_bits(&o.precision) == FL_PRECISION_DOUBLE) batch_impl<double>(mesh, *desc, o, mode, out_results);
        else                                    batch_impl<float>(mesh, *desc, o, mode, out_results);
    });
}

/* ------------------------------------------------------------------------- */
/* Rasters                                                                    */
/* ------------------------------------------------------------------------- */

fl_status fl_project_image(const fl_mesh* mesh_c, const double view_dir[3],
                           const double* field, size_t field_len,
                           const fl_options* opts, fl_image** out_image,
                           fl_result* out_result) {
    if (!mesh_c || !view_dir || !out_image)
        return fail(FL_ERR_INVALID_ARGUMENT, "mesh, view_dir and out_image must all be non-NULL");
    if (!dir_ok(view_dir))
        return fail(FL_ERR_INVALID_ARGUMENT, "view_dir must be a non-zero, finite vector");
    *out_image = nullptr;

    fl_options storage;
    const fl_options& o = effective(opts, storage);
    if (fl_status s = check_options(o)) return s;

    fl_mesh* mesh = const_cast<fl_mesh*>(mesh_c);
    ValueMode mode = MODE_NONE;
    if (field) {
        if (field_len == 0) return fail(FL_ERR_INVALID_ARGUMENT, "field is non-NULL but field_len is 0");
        if (fl_status s = resolve_field_mode(mesh, field_len, o.field_mode, mode)) return s;
    }

    return guard(FL_ERR_NUMERIC, [&]() {
        std::unique_ptr<fl_image> img(new fl_image());
        if (enum_bits(&o.precision) == FL_PRECISION_DOUBLE)
            image_impl<double>(mesh, view_dir, field, field_len, mode, o, img.get(), out_result);
        else
            image_impl<float>(mesh, view_dir, field, field_len, mode, o, img.get(), out_result);
        *out_image = img.release();
    });
}

int32_t fl_image_width (const fl_image* i) { return i ? i->w : 0; }
int32_t fl_image_height(const fl_image* i) { return i ? i->h : 0; }

const uint8_t* fl_image_mask(const fl_image* i) {
    return (i && !i->mask.empty()) ? i->mask.data() : nullptr;
}

const double* fl_image_values(const fl_image* i) {
    return (i && i->has_values && !i->values.empty()) ? i->values.data() : nullptr;
}

fl_status fl_image_write_ppm(const fl_image* i, const char* path) {
    if (!i || !path) return fail(FL_ERR_INVALID_ARGUMENT, "image and path must both be non-NULL");
    if (i->w <= 0 || i->h <= 0)
        return fail(FL_ERR_DIMENSION, "this view covered no pixels, so there is no image to write");
    return guard(FL_ERR_IO, [&]() {
        write_ppm<double>(path, i->w, i->h, i->mask.data(),
                          i->has_values ? i->values.data() : nullptr, i->min_v, i->max_v);
    });
}

void fl_image_destroy(fl_image* i) { delete i; }

/* ------------------------------------------------------------------------- */
/* Field files                                                                */
/* ------------------------------------------------------------------------- */

fl_status fl_field_matrix_load(const char* path, const fl_mesh* mesh, fl_field_mode mode,
                               fl_field_matrix** out_matrix) {
    if (!path || !mesh || !out_matrix)
        return fail(FL_ERR_INVALID_ARGUMENT, "path, mesh and out_matrix must all be non-NULL");
    *out_matrix = nullptr;

    ValueMode forced = MODE_NONE;
    if (mode == FL_FIELD_NODE)      forced = MODE_NODE;
    else if (mode == FL_FIELD_FACE) forced = MODE_FACE;
    else if (mode != FL_FIELD_AUTO)
        return fail(FL_ERR_INVALID_ARGUMENT, "mode is not a valid fl_field_mode value");

    return guard(FL_ERR_PARSE, [&]() {
        std::unique_ptr<fl_field_matrix> fm(new fl_field_matrix());
        load_matrix_into(path, mesh->md, fm->m, forced);
        *out_matrix = fm.release();
    });
}

size_t fl_field_matrix_rows(const fl_field_matrix* m) { return m ? m->m.nrows : 0; }
size_t fl_field_matrix_cols(const fl_field_matrix* m) { return m ? m->m.ncols : 0; }

fl_field_mode fl_field_matrix_mode(const fl_field_matrix* m) {
    if (!m) return FL_FIELD_AUTO;
    switch (m->m.mode) {
        case MODE_NODE: return FL_FIELD_NODE;
        case MODE_FACE: return FL_FIELD_FACE;
        default:        return FL_FIELD_AUTO;
    }
}

const double* fl_field_matrix_data(const fl_field_matrix* m) {
    return (m && !m->m.data.empty()) ? m->m.data.data() : nullptr;
}

fl_status fl_field_matrix_copy_column(const fl_field_matrix* m, size_t col,
                                      double* out, size_t out_len) {
    if (!m || !out) return fail(FL_ERR_INVALID_ARGUMENT, "matrix and out must both be non-NULL");
    if (col >= m->m.ncols)
        return fail(FL_ERR_DIMENSION, "column " + std::to_string(col) + " is out of range (" +
                    std::to_string(m->m.ncols) + " column(s))");
    if (out_len < m->m.nrows)
        return fail(FL_ERR_DIMENSION, "out needs room for " + std::to_string(m->m.nrows) +
                    " doubles, got " + std::to_string(out_len));
    for (size_t r = 0; r < m->m.nrows; ++r) out[r] = m->m.data[r * m->m.ncols + col];
    return ok();
}

void fl_field_matrix_destroy(fl_field_matrix* m) { delete m; }

} // extern "C"
