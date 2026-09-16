/*
 * flatland.h — public C ABI for FlatLand.
 *
 * FlatLand computes the projected (visible) surface area of a triangle mesh from
 * arbitrary view directions by orthographic projection and software rasterization.
 * Given a scalar field defined per vertex or per face it also reports, over the
 * visible projection: the mean, the extremes, and the area integral
 *
 *     I = integral of f dA = sum over covered pixels of (value * pixel_area)
 *
 * The integral is deliberately generic. Supply a radiance field and I is a radiant
 * intensity; supply a pressure field and I is a force. FlatLand itself attaches no
 * physical meaning to the numbers it is given.
 *
 * This header is the stable interface. It is plain C (no C++ types cross the
 * boundary), so it is consumable from C, C++, Python via ctypes, MATLAB via
 * loadlibrary, Julia via ccall, and anything else that speaks the C ABI.
 *
 * -------------------------------------------------------------------------
 * CONVENTIONS
 *
 * Coordinates   Vertices are (x, y, z) triples of double, interleaved.
 * Faces         Triangles are triples of ZERO-BASED int32 vertex indices,
 *               interleaved. (Note that OBJ files on disk are 1-based; the
 *               loader converts. This API is 0-based throughout.)
 * View          A view direction is the direction the camera LOOKS ALONG, so
 *               the visible surface is the one facing back toward -direction.
 *               Magnitude is irrelevant; the vector is normalized internally.
 * Precision     Every value crossing this boundary is double, regardless of the
 *               precision the engine is asked to compute in. fl_precision selects
 *               the internal working type only.
 * Fields        A field has one value per vertex (node mode) or per face (face
 *               mode). Node fields are interpolated across each triangle; face
 *               fields are constant over it.
 * Ownership     Every array the caller passes in is borrowed for the duration of
 *               the call only; FlatLand copies whatever it needs to retain.
 *               Every handle FlatLand returns must be released with its
 *               corresponding _destroy function.
 * Threads       An fl_mesh is immutable once created and may be shared freely
 *               across threads. fl_project_batch parallelizes internally.
 *               fl_last_error is thread-local.
 *
 * -------------------------------------------------------------------------
 * MINIMAL EXAMPLE
 *
 *     double verts[] = { 0,0,0,  1,0,0,  0,1,0 };
 *     int32_t faces[] = { 0, 1, 2 };
 *     fl_mesh* mesh = NULL;
 *     if (fl_mesh_create(verts, 3, faces, 1, &mesh) != FL_OK)
 *         fprintf(stderr, "%s\n", fl_last_error());
 *
 *     fl_options opts;
 *     fl_options_init(&opts);
 *     opts.resolution = 0.01;
 *
 *     double view[3] = { 0, 0, 1 };
 *     fl_result r;
 *     fl_project(mesh, view, NULL, 0, &opts, &r);
 *     printf("area = %g\n", r.area);
 *
 *     fl_mesh_destroy(mesh);
 */

#ifndef FLATLAND_H
#define FLATLAND_H

#include <stddef.h>
#include <stdint.h>

/* ------------------------------------------------------------------------- */
/* Version                                                                    */
/* ------------------------------------------------------------------------- */

#define FL_VERSION_MAJOR 3
#define FL_VERSION_MINOR 0
#define FL_VERSION_PATCH 0

/* ------------------------------------------------------------------------- */
/* Linkage                                                                    */
/* ------------------------------------------------------------------------- */

/*
 * Define FLATLAND_STATIC when linking the static library. The shared library
 * defines FLATLAND_BUILD_SHARED while building itself; consumers need define
 * nothing on ELF/Mach-O targets, and nothing on Windows either, since the
 * default below is the import side.
 */
#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(FLATLAND_STATIC)
#    define FL_API
#  elif defined(FLATLAND_BUILD_SHARED)
#    define FL_API __declspec(dllexport)
#  else
#    define FL_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) && __GNUC__ >= 4
#  define FL_API __attribute__((visibility("default")))
#else
#  define FL_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Status codes                                                               */
/* ------------------------------------------------------------------------- */

typedef enum fl_status {
    FL_OK                   = 0,  /* success                                        */
    FL_ERR_INVALID_ARGUMENT = 1,  /* NULL where a pointer is required, bad enum, ... */
    FL_ERR_OUT_OF_MEMORY    = 2,  /* allocation failed                              */
    FL_ERR_IO               = 3,  /* file could not be opened or written            */
    FL_ERR_PARSE            = 4,  /* file opened but its contents are malformed     */
    FL_ERR_DIMENSION        = 5,  /* size mismatch: field rows vs mesh, column range */
    FL_ERR_NUMERIC          = 6,  /* non-finite input, or a raster too large to build */
    FL_ERR_UNSUPPORTED      = 7,  /* valid request this build cannot serve          */
    FL_ERR_INTERNAL         = 8   /* a bug in FlatLand; please report               */
} fl_status;

/* Static, human-readable name of a status code. Never NULL, never freed. */
FL_API const char* fl_status_string(fl_status status);

/*
 * Detail for the most recent failing call ON THE CALLING THREAD. Returns a
 * pointer valid until the next FlatLand call on that same thread, or an empty
 * string if the last call succeeded. Never NULL.
 */
FL_API const char* fl_last_error(void);

/* Compile-time version this library was built from. */
FL_API void        fl_version(int* major, int* minor, int* patch);
FL_API const char* fl_version_string(void);

/* ------------------------------------------------------------------------- */
/* Enumerations                                                               */
/* ------------------------------------------------------------------------- */

typedef enum fl_precision {
    FL_PRECISION_FLOAT  = 0,  /* faster; adequate for coverage and smooth fields */
    FL_PRECISION_DOUBLE = 1   /* slower; use for wide dynamic range or fine rasters */
} fl_precision;

typedef enum fl_field_mode {
    FL_FIELD_AUTO = 0,  /* infer from row count; error if vertex and face counts tie */
    FL_FIELD_NODE = 1,  /* one value per vertex, interpolated across each triangle  */
    FL_FIELD_FACE = 2   /* one value per face, constant across it                   */
} fl_field_mode;

/* ------------------------------------------------------------------------- */
/* Mesh                                                                       */
/* ------------------------------------------------------------------------- */

typedef struct fl_mesh fl_mesh;

/*
 * Build a mesh from caller-owned arrays.
 *
 *   vertices     vertex_count * 3 doubles, xyz interleaved
 *   faces        face_count * 3 int32, zero-based vertex indices, interleaved
 *
 * Both arrays are copied. Indices are validated against vertex_count and
 * coordinates are checked for finiteness; either failing is FL_ERR_DIMENSION or
 * FL_ERR_NUMERIC respectively. On success *out_mesh receives a handle the caller
 * must release with fl_mesh_destroy.
 */
FL_API fl_status fl_mesh_create(const double* vertices, size_t vertex_count,
                                const int32_t* faces,   size_t face_count,
                                fl_mesh** out_mesh);

/* As fl_mesh_create, for callers holding single-precision geometry. */
FL_API fl_status fl_mesh_create_f32(const float* vertices, size_t vertex_count,
                                    const int32_t* faces,  size_t face_count,
                                    fl_mesh** out_mesh);

/*
 * Load a mesh from an OBJ or STL file (binary or ASCII STL, auto-detected).
 * Polygonal OBJ faces are fan-triangulated. STL carries no shared vertices, so
 * a mesh loaded from STL has three vertices per triangle and its natural field
 * mode is face.
 */
FL_API fl_status fl_mesh_load(const char* path, fl_mesh** out_mesh);

FL_API void   fl_mesh_destroy(fl_mesh* mesh);
FL_API size_t fl_mesh_vertex_count(const fl_mesh* mesh);
FL_API size_t fl_mesh_face_count(const fl_mesh* mesh);

/*
 * Read geometry back out — chiefly so a caller who used fl_mesh_load can hand the
 * arrays to their own code. out_len is in ELEMENTS and must be at least
 * vertex_count*3 / face_count*3 respectively.
 *
 * Vertices come back in the original input frame. (FlatLand recenters geometry
 * internally for numerical conditioning; that is not visible here.)
 */
FL_API fl_status fl_mesh_copy_vertices(const fl_mesh* mesh, double*  out, size_t out_len);
FL_API fl_status fl_mesh_copy_faces   (const fl_mesh* mesh, int32_t* out, size_t out_len);

/* ------------------------------------------------------------------------- */
/* Options                                                                    */
/* ------------------------------------------------------------------------- */

typedef struct fl_options {
    double        resolution;  /* pixel edge length in mesh units; default 0.001   */
    int32_t       cull;        /* 1 = backface cull (default), 0 = render all faces */
    fl_precision  precision;   /* default FL_PRECISION_FLOAT                        */
    fl_field_mode field_mode;  /* default FL_FIELD_AUTO                             */
    int32_t       threads;     /* batch worker threads; 0 = one per core (default)  */
    uint32_t      reserved[8]; /* must be zero; reserved for future fields          */
} fl_options;

/*
 * Fill opts with defaults. ALWAYS call this before setting fields — it zeroes
 * the reserved words, which future versions rely on to stay ABI-compatible.
 */
FL_API void fl_options_init(fl_options* opts);

/* ------------------------------------------------------------------------- */
/* Results                                                                    */
/* ------------------------------------------------------------------------- */

typedef struct fl_result {
    double  area;            /* covered_pixels * resolution^2                     */
    double  average;         /* mean field value over covered pixels              */
    double  integral;        /* sum(value * pixel_area) — the area integral       */
    double  min;             /* least field value over covered pixels             */
    double  max;             /* greatest field value over covered pixels          */
    int64_t covered_pixels;  /* pixels whose center fell on a visible triangle    */
    int32_t width;           /* raster width actually used                        */
    int32_t height;          /* raster height actually used                       */
    int32_t has_field;       /* 1 if a field was supplied for this view           */
    /*
     * 1 if covered_pixels > 0 AND has_field. When 0, the four field statistics
     * above are NOT measurements — they are left at zero and must be ignored.
     * A view can legitimately cover no pixels: geometry thinner than one pixel,
     * or a mesh seen exactly edge-on. Check this before consuming a batch.
     */
    int32_t  has_stats;
    uint32_t reserved[8];
} fl_result;

/* ------------------------------------------------------------------------- */
/* Projection                                                                 */
/* ------------------------------------------------------------------------- */

/*
 * Project one view.
 *
 *   view_dir    3 doubles, the direction the camera looks along; need not be a
 *               unit vector, but must be non-zero and finite
 *   field       field_len doubles, or NULL for a geometry-only run
 *   field_len   vertex_count for a node field, face_count for a face field;
 *               with opts->field_mode == FL_FIELD_AUTO the length decides, and
 *               a mesh whose vertex and face counts are equal is ambiguous and
 *               must be disambiguated by setting field_mode explicitly
 *   opts        may be NULL, in which case defaults apply
 */
FL_API fl_status fl_project(const fl_mesh* mesh,
                            const double view_dir[3],
                            const double* field, size_t field_len,
                            const fl_options* opts,
                            fl_result* out_result);

/*
 * A time series: many views over one mesh, evaluated in parallel.
 *
 * This is the throughput path. The mesh is projected once per view but parsed,
 * validated and recentered only once overall, and the field matrix is read in
 * place rather than re-extracted per worker.
 */
typedef struct fl_batch_desc {
    const double* views;       /* view_count * 3 doubles, directions interleaved  */
    size_t        view_count;

    /*
     * Optional field matrix, row-major, field_rows x field_cols.
     * One ROW per mesh entity, one COLUMN per timestep — so an entire time
     * series lives in a single array. field_rows must equal the mesh vertex
     * count (node) or face count (face); see opts->field_mode.
     * Set field_matrix to NULL for a geometry-only batch.
     */
    const double* field_matrix;
    size_t        field_rows;
    size_t        field_cols;

    /*
     * Which column each view reads, view_count entries. NULL means every view
     * reads column 0. Entries must be < field_cols.
     */
    const int32_t* field_columns;

    /*
     * Per-view pixel size, view_count entries. NULL means every view uses
     * opts->resolution. Entries must be positive and finite.
     */
    const double* resolutions;

    uint32_t reserved[8];      /* must be zero */
} fl_batch_desc;

/* Fill desc with zeroes/defaults. Call before populating, as with fl_options_init. */
FL_API void fl_batch_desc_init(fl_batch_desc* desc);

/*
 * Run a batch. out_results must have room for desc->view_count entries and is
 * filled in view order regardless of the order workers happen to finish in.
 *
 * If any view fails the whole call fails and reports the LOWEST-indexed failure,
 * so the error is reproducible run to run. out_results is then unspecified.
 */
FL_API fl_status fl_project_batch(const fl_mesh* mesh,
                                  const fl_batch_desc* desc,
                                  const fl_options* opts,
                                  fl_result* out_results);

/* Convert azimuth/elevation in degrees to a unit direction. Azimuth sweeps around
 * +Z measured from +X; elevation rises from the XY plane toward +Z. */
FL_API void fl_angle_to_dir(double azimuth_deg, double elevation_deg, double out_dir[3]);

/* ------------------------------------------------------------------------- */
/* Rasters                                                                    */
/* ------------------------------------------------------------------------- */

/*
 * A rendered view: the coverage mask and, where a field was supplied, the
 * interpolated per-pixel values. Row 0 is the BOTTOM row in mesh space; index a
 * pixel as [y * width + x].
 */
typedef struct fl_image fl_image;

/*
 * Project one view and keep the raster. Arguments match fl_project; *out_image
 * receives a handle to release with fl_image_destroy. out_result may be NULL if
 * only the raster is wanted.
 *
 * A view covering no pixels yields a valid image of width 0 and height 0 rather
 * than a stale or partially-filled buffer.
 */
FL_API fl_status fl_project_image(const fl_mesh* mesh,
                                  const double view_dir[3],
                                  const double* field, size_t field_len,
                                  const fl_options* opts,
                                  fl_image** out_image,
                                  fl_result* out_result);

FL_API int32_t fl_image_width (const fl_image* image);
FL_API int32_t fl_image_height(const fl_image* image);

/*
 * Borrowed views of the raster, valid until fl_image_destroy. Each spans
 * width*height elements. fl_image_values returns NULL when the view carried no
 * field. Neither buffer is copied, so a binding can wrap them in a zero-copy
 * array as long as it keeps the fl_image alive.
 */
FL_API const uint8_t* fl_image_mask  (const fl_image* image);  /* 1 = covered */
FL_API const double*  fl_image_values(const fl_image* image);

/* Write a false-color PPM (P6). Covered pixels are ramped across the field range;
 * a fieldless image is written as a white silhouette. */
FL_API fl_status fl_image_write_ppm(const fl_image* image, const char* path);

FL_API void fl_image_destroy(fl_image* image);

/* ------------------------------------------------------------------------- */
/* Field files                                                                */
/* ------------------------------------------------------------------------- */

/*
 * Field matrices can also be read from disk, in the same whitespace-separated
 * text format the CLI accepts: one row per mesh entity, one column per timestep,
 * '#' starting a comment. This exists so a caller can reuse an existing on-disk
 * time series without reimplementing the parser, then feed the column straight
 * into fl_project or hand the whole matrix to fl_project_batch.
 */
typedef struct fl_field_matrix fl_field_matrix;

FL_API fl_status fl_field_matrix_load(const char* path, const fl_mesh* mesh,
                                      fl_field_mode mode, fl_field_matrix** out_matrix);

FL_API size_t        fl_field_matrix_rows(const fl_field_matrix* m);
FL_API size_t        fl_field_matrix_cols(const fl_field_matrix* m);
FL_API fl_field_mode fl_field_matrix_mode(const fl_field_matrix* m);

/* Borrowed row-major view, rows*cols doubles, valid until destroy. Suitable for
 * passing directly as fl_batch_desc::field_matrix. */
FL_API const double* fl_field_matrix_data(const fl_field_matrix* m);

/* Copy one column out. out_len must be at least fl_field_matrix_rows(m). */
FL_API fl_status fl_field_matrix_copy_column(const fl_field_matrix* m, size_t col,
                                             double* out, size_t out_len);

FL_API void fl_field_matrix_destroy(fl_field_matrix* m);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* FLATLAND_H */
