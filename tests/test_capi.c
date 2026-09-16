// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

/*
 * test_capi.c — exercises the public C ABI as a C consumer would.
 *
 * Deliberately compiled as C, not C++: the header's whole purpose is to be
 * consumable from C, ctypes and MATLAB, and building it with a C++ compiler
 * would hide exactly the mistakes that break those callers.
 *
 * Usage: test_capi [mesh.obj [field.txt]]
 * Exits non-zero if any check fails.
 */

#include "flatland.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int passed = 0, failed = 0;

static void ok_(const char* what) { printf("  PASS %s\n", what); ++passed; }
static void bad_(const char* what, const char* detail) {
    printf("  FAIL %s (%s)\n", what, detail); ++failed;
}

static void check(int cond, const char* what) {
    if (cond) ok_(what); else bad_(what, "condition false");
}

static void near_(double a, double e, double tol, const char* what) {
    char buf[160];
    if (fabs(a - e) <= tol) {
        snprintf(buf, sizeof buf, "%s (%.9g ~= %.9g)", what, a, e);
        ok_(buf);
    } else {
        snprintf(buf, sizeof buf, "got %.9g, want %.9g", a, e);
        bad_(what, buf);
    }
}

static void expect_err(fl_status s, const char* what) {
    char buf[220];
    if (s == FL_OK) { bad_(what, "expected a failure, got FL_OK"); return; }
    /* Every failure must leave a usable message behind, not an empty string. */
    if (fl_last_error()[0] == '\0') { bad_(what, "failed but left no error message"); return; }
    snprintf(buf, sizeof buf, "%s [%s: %.90s]", what, fl_status_string(s), fl_last_error());
    ok_(buf);
}

/* A closed unit box spanning [0,1]^3, CCW outward normals.
 * 8 vertices and 12 faces, so node and face fields are never ambiguous. */
static const double BOX_V[24] = {
    0,0,0,  1,0,0,  1,1,0,  0,1,0,
    0,0,1,  1,0,1,  1,1,1,  0,1,1
};
static const int32_t BOX_F[36] = {
    0,3,2,  0,2,1,   4,5,6,  4,6,7,
    0,1,5,  0,5,4,   1,2,6,  1,6,5,
    2,3,7,  2,7,6,   3,0,4,  3,4,7
};

int main(int argc, char** argv) {
    fl_mesh* box = NULL;
    fl_options o;
    fl_result r;
    fl_status s;

    printf("\n== version and defaults ==\n");
    {
        int maj = -1, min = -1, pat = -1;
        fl_version(&maj, &min, &pat);
        check(maj == FL_VERSION_MAJOR && min == FL_VERSION_MINOR && pat == FL_VERSION_PATCH,
              "fl_version matches the header macros");
        check(fl_version_string() != NULL && fl_version_string()[0] != '\0',
              "fl_version_string is non-empty");
        check(strcmp(fl_status_string(FL_OK), "ok") == 0, "fl_status_string(FL_OK)");
    }

    fl_options_init(&o);
    check(o.resolution == 0.001 && o.cull == 1 &&
          o.precision == FL_PRECISION_FLOAT && o.field_mode == FL_FIELD_AUTO && o.threads == 0,
          "fl_options_init sets the documented defaults");
    {
        int i, zeroed = 1;
        for (i = 0; i < 8; ++i) if (o.reserved[i] != 0) zeroed = 0;
        check(zeroed, "fl_options_init zeroes the reserved words");
    }

    printf("\n== mesh from caller arrays ==\n");
    s = fl_mesh_create(BOX_V, 8, BOX_F, 12, &box);
    check(s == FL_OK && box != NULL, "fl_mesh_create succeeds");
    check(fl_mesh_vertex_count(box) == 8, "vertex count round-trips");
    check(fl_mesh_face_count(box) == 12, "face count round-trips");

    {   /* Recentering is an internal detail; the caller's frame comes back. */
        double back[24];
        int i, same = 1;
        s = fl_mesh_copy_vertices(box, back, 24);
        for (i = 0; i < 24; ++i) if (fabs(back[i] - BOX_V[i]) > 1e-12) same = 0;
        check(s == FL_OK && same, "fl_mesh_copy_vertices returns the original frame");
    }
    {
        int32_t back[36];
        int i, same = 1;
        s = fl_mesh_copy_faces(box, back, 36);
        for (i = 0; i < 36; ++i) if (back[i] != BOX_F[i]) same = 0;
        check(s == FL_OK && same, "fl_mesh_copy_faces round-trips");
    }

    printf("\n== projection ==\n");
    o.resolution = 0.002;
    o.precision  = FL_PRECISION_DOUBLE;
    {
        const double view[3] = {1, 0, 0};
        s = fl_project(box, view, NULL, 0, &o, &r);
        check(s == FL_OK, "fl_project with no field");
        near_(r.area, 1.0, 0.01, "unit box projects to area 1");
        check(r.has_field == 0 && r.has_stats == 0, "no field means no statistics");
        check(r.covered_pixels > 0 && r.width > 0 && r.height > 0, "raster dimensions are reported");
    }

    {   /* Field = vertex x. Looking along +X the visible surface is x = 0, so
           the mean must be 0 — this is the assertion that catches an inverted
           cull, and it needs a field that VARIES across the mesh. */
        double fx[8];
        const double view[3] = {1, 0, 0};
        int i;
        for (i = 0; i < 8; ++i) fx[i] = BOX_V[i*3];

        s = fl_project(box, view, fx, 8, &o, &r);
        check(s == FL_OK, "fl_project with a node field");
        check(r.has_field == 1 && r.has_stats == 1, "a covered field view has statistics");
        near_(r.average, 0.0, 1e-6, "culled view resolves to the NEAR surface");
        near_(r.min, 0.0, 1e-6, "field min over the near surface");
        near_(r.max, 0.0, 1e-6, "field max over the near surface");

        o.cull = 0;
        s = fl_project(box, view, fx, 8, &o, &r);
        near_(r.average, 0.0, 1e-6, "culling agrees with no culling");
        o.cull = 1;
    }

    {   /* A constant face field: the area integral is the constant times the area. */
        double ff[12];
        const double view[3] = {0, 0, 1};
        int i;
        for (i = 0; i < 12; ++i) ff[i] = 4.0;
        s = fl_project(box, view, ff, 12, &o, &r);
        check(s == FL_OK, "fl_project with a face field");
        near_(r.average, 4.0, 1e-9, "constant face field average");
        near_(r.integral, 4.0 * r.area, 1e-9, "integral == constant * area");
    }

    {   /* float and double must agree to well within pixel discretization. */
        double a_f, a_d;
        const double view[3] = {1, 1, 1};
        o.precision = FL_PRECISION_FLOAT;
        fl_project(box, view, NULL, 0, &o, &r); a_f = r.area;
        o.precision = FL_PRECISION_DOUBLE;
        fl_project(box, view, NULL, 0, &o, &r); a_d = r.area;
        near_(a_f, a_d, 1e-3, "float and double agree on area");
    }

    printf("\n== batch ==\n");
    {
        /* 3 views, one field matrix of 8 rows x 3 columns holding 10 / 20 / 30,
           each view reading its own column. */
        double views[9]   = { 1,0,0,  0,1,0,  0,0,1 };
        double matrix[24];
        int32_t cols[3]   = {0, 1, 2};
        double res[3]     = {0.004, 0.004, 0.004};
        fl_result out[3];
        fl_batch_desc d;
        int i;
        for (i = 0; i < 8; ++i) { matrix[i*3+0]=10; matrix[i*3+1]=20; matrix[i*3+2]=30; }

        fl_batch_desc_init(&d);
        d.views = views; d.view_count = 3;
        d.field_matrix = matrix; d.field_rows = 8; d.field_cols = 3;
        d.field_columns = cols;
        d.resolutions = res;

        s = fl_project_batch(box, &d, &o, out);
        check(s == FL_OK, "fl_project_batch succeeds");
        near_(out[0].average, 10.0, 1e-9, "view 0 reads column 0");
        near_(out[1].average, 20.0, 1e-9, "view 1 reads column 1");
        near_(out[2].average, 30.0, 1e-9, "view 2 reads column 2");
        near_(out[0].area, 1.0, 0.02, "batch view 0 area");

        /* Results must not depend on how many workers ran them. */
        {
            fl_result single[3];
            fl_options o1 = o;
            o1.threads = 1;
            fl_project_batch(box, &d, &o1, single);
            check(single[0].average == out[0].average &&
                  single[1].average == out[1].average &&
                  single[2].average == out[2].average,
                  "batch results are independent of thread count");
        }

        d.field_columns = NULL;
        s = fl_project_batch(box, &d, &o, out);
        check(s == FL_OK, "NULL field_columns means every view reads column 0");
        near_(out[2].average, 10.0, 1e-9, "...confirmed on the last view");
    }

    printf("\n== rasters ==\n");
    {
        fl_image* img = NULL;
        const double view[3] = {1, 0, 0};
        double fx[8];
        int i;
        for (i = 0; i < 8; ++i) fx[i] = BOX_V[i*3+1];   /* field = vertex y */

        s = fl_project_image(box, view, fx, 8, &o, &img, &r);
        check(s == FL_OK && img != NULL, "fl_project_image succeeds");
        check(fl_image_width(img) == r.width && fl_image_height(img) == r.height,
              "image dimensions match the result");
        check(fl_image_mask(img) != NULL, "mask is exposed");
        check(fl_image_values(img) != NULL, "values are exposed for a field view");
        {
            const unsigned char* mask = fl_image_mask(img);
            long n = (long)fl_image_width(img) * fl_image_height(img), covered = 0, k;
            for (k = 0; k < n; ++k) if (mask[k]) ++covered;
            check(covered == r.covered_pixels, "mask agrees with covered_pixels");
        }
        {
            const char* path = "capi_test_out.ppm";
            FILE* fp;
            s = fl_image_write_ppm(img, path);
            check(s == FL_OK, "fl_image_write_ppm succeeds");
            fp = fopen(path, "rb");
            if (fp) {
                char hdr[3] = {0,0,0};
                size_t got = fread(hdr, 1, 2, fp);
                fclose(fp);
                check(got == 2 && hdr[0] == 'P' && hdr[1] == '6', "the PPM has a P6 header");
                remove(path);
            } else {
                bad_("the PPM has a P6 header", "could not reopen the file");
            }
        }
        fl_image_destroy(img);

        /* A view that covers nothing must still yield a valid, empty image
           rather than a stale raster from a previous call. */
        {
            fl_mesh* tri = NULL;
            const double tv[9] = {0,0,0,  1,0,0,  1,1,0};
            const int32_t tf[3] = {0,1,2};
            const double edge_on[3] = {1, 0, 0};
            fl_image* empty = NULL;
            fl_mesh_create(tv, 3, tf, 1, &tri);
            s = fl_project_image(tri, edge_on, NULL, 0, &o, &empty, &r);
            check(s == FL_OK, "an edge-on view is not an error");
            check(fl_image_width(empty) == 0 && fl_image_height(empty) == 0,
                  "an empty view yields a 0x0 image, not a stale one");
            check(r.covered_pixels == 0 && r.has_stats == 0,
                  "an empty view reports no coverage and no statistics");
            expect_err(fl_image_write_ppm(empty, "should_not_exist.ppm"),
                       "writing a PPM for an empty view is refused");
            fl_image_destroy(empty);
            fl_mesh_destroy(tri);
        }
    }

    printf("\n== zero-coverage statistics ==\n");
    {
        /* A sliver thinner than one pixel: covered_pixels is 0, so the field
           statistics are not measurements and has_stats must say so. */
        const double sv[9] = {0,0,0,  1,0,0,  1,0.000001,0};
        const int32_t sf[3] = {0,1,2};
        const double field[3] = {-5, -3, -1};
        const double view[3] = {0,0,-1};
        fl_mesh* sliver = NULL;
        fl_options so = o;
        so.resolution = 0.1;
        so.cull = 0;
        fl_mesh_create(sv, 3, sf, 1, &sliver);
        s = fl_project(sliver, view, field, 3, &so, &r);
        check(s == FL_OK, "a zero-coverage view is not an error");
        check(r.covered_pixels == 0, "the sliver covers no pixels");
        check(r.has_field == 1 && r.has_stats == 0,
              "has_field is set but has_stats is not");
        check(r.average == 0 && r.min == 0 && r.max == 0 && r.integral == 0,
              "the untrustworthy statistics are left at zero");
        fl_mesh_destroy(sliver);
    }

    printf("\n== error handling ==\n");
    {
        const double view[3] = {1,0,0};
        const double zero[3] = {0,0,0};
        double bad_field[7] = {0,0,0,0,0,0,0};
        fl_mesh* m2 = NULL;
        fl_options bo;

        expect_err(fl_mesh_create(NULL, 8, BOX_F, 12, &m2), "NULL vertices rejected");
        expect_err(fl_mesh_create(BOX_V, 0, BOX_F, 12, &m2), "zero vertex count rejected");
        expect_err(fl_mesh_create(BOX_V, 8, BOX_F, 0, &m2),  "zero face count rejected");
        expect_err(fl_project(NULL, view, NULL, 0, &o, &r),  "NULL mesh rejected");
        expect_err(fl_project(box, view, NULL, 0, &o, NULL), "NULL out_result rejected");
        expect_err(fl_project(box, zero, NULL, 0, &o, &r),   "zero view direction rejected");
        expect_err(fl_project(box, view, bad_field, 7, &o, &r), "wrong-length field rejected");
        expect_err(fl_mesh_load("no_such_file_here.obj", &m2), "missing mesh file rejected");

        {   /* An index past the end of the vertex array must be caught, not
               dereferenced. */
            const int32_t oob[3] = {0, 1, 99};
            expect_err(fl_mesh_create(BOX_V, 8, oob, 1, &m2), "out-of-range face index rejected");
        }
        {
            const double inf_dir[3] = {1e308*10, 0, 1};
            expect_err(fl_project(box, inf_dir, NULL, 0, &o, &r), "non-finite view direction rejected");
        }

        bo = o; bo.resolution = 0;
        expect_err(fl_project(box, view, NULL, 0, &bo, &r), "zero resolution rejected");
        bo = o; bo.resolution = -1;
        expect_err(fl_project(box, view, NULL, 0, &bo, &r), "negative resolution rejected");
        bo = o; bo.resolution = 1e-9;
        expect_err(fl_project(box, view, NULL, 0, &bo, &r), "an impossibly fine resolution is refused");
        bo = o; bo.precision = (fl_precision)42;
        expect_err(fl_project(box, view, NULL, 0, &bo, &r), "invalid precision rejected");

        {   /* A successful call must clear the previous error. */
            fl_status good = fl_project(box, view, NULL, 0, &o, &r);
            check(good == FL_OK && fl_last_error()[0] == '\0',
                  "a successful call clears the error message");
        }
        {
            fl_batch_desc d;
            fl_result out1;
            int32_t badcol[1] = {9};
            double v1[3] = {1,0,0};
            double mat[8];
            int i;
            for (i = 0; i < 8; ++i) mat[i] = 1.0;
            fl_batch_desc_init(&d);
            d.views = v1; d.view_count = 1;
            d.field_matrix = mat; d.field_rows = 8; d.field_cols = 1;
            d.field_columns = badcol;
            expect_err(fl_project_batch(box, &d, &o, &out1), "out-of-range field column rejected");
        }
    }

    printf("\n== file-backed mesh and field ==\n");
    if (argc > 1) {
        fl_mesh* fm = NULL;
        s = fl_mesh_load(argv[1], &fm);
        check(s == FL_OK && fm != NULL, "fl_mesh_load reads an OBJ");
        if (fm) {
            const double view[3] = {1,0,0};
            fl_options fo = o;
            fo.resolution = 0.005;
            check(fl_mesh_vertex_count(fm) > 0 && fl_mesh_face_count(fm) > 0,
                  "the loaded mesh has geometry");
            s = fl_project(fm, view, NULL, 0, &fo, &r);
            check(s == FL_OK, "the loaded mesh projects");
            near_(r.area, 1.0, 0.02, "the cube example still projects to area 1");

            if (argc > 2) {
                fl_field_matrix* fmx = NULL;
                s = fl_field_matrix_load(argv[2], fm, FL_FIELD_AUTO, &fmx);
                check(s == FL_OK && fmx != NULL, "fl_field_matrix_load reads a field file");
                if (fmx) {
                    const double* data = fl_field_matrix_data(fmx);
                    size_t rows = fl_field_matrix_rows(fmx);
                    size_t cols = fl_field_matrix_cols(fmx);
                    check(rows == fl_mesh_vertex_count(fm), "the matrix has one row per vertex");
                    check(cols >= 1 && data != NULL, "the matrix exposes its data");
                    {
                        double* col = (double*)malloc(rows * sizeof(double));
                        s = fl_field_matrix_copy_column(fmx, 0, col, rows);
                        check(s == FL_OK, "fl_field_matrix_copy_column succeeds");
                        expect_err(fl_field_matrix_copy_column(fmx, cols + 5, col, rows),
                                   "copying an out-of-range column is refused");
                        free(col);
                    }
                    fl_field_matrix_destroy(fmx);
                }
            }
            fl_mesh_destroy(fm);
        }
    } else {
        printf("  (skipped: no mesh path given)\n");
    }

    fl_mesh_destroy(box);
    fl_mesh_destroy(NULL);   /* destroying NULL must be safe */
    ok_("destroying a NULL handle is safe");

    printf("\nC API: %d passed, %d failed\n", passed, failed);
    return failed == 0 ? 0 : 1;
}
