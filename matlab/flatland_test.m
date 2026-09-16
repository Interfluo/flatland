function failures = flatland_test()
%FLATLAND_TEST  Self-check for the FlatLand MATLAB binding.
%
%   flatland_test            run every check, print a summary, error if any fail
%   n = flatland_test()      as above but return the failure count instead of
%                            erroring, for scripted use
%
%   Run it from the matlab/ folder, or with that folder on the MATLAB path:
%       >> cd /path/to/flatland/matlab
%       >> flatland_test
%
%   WHERE THE EXPECTED NUMBERS COME FROM
%   Every numeric expectation below was produced by compiling and running the
%   equivalent C against the same libflatland.so, so these are cross-checks
%   against the C ABI rather than against this binding's own behaviour. The
%   headline values:
%
%     unit box, view [1 0 0], resolution 0.002, double precision
%         area = 1            covered = 250000     raster 502 x 502
%     node field = vertex x, same view
%         average = min = max = 0      (the culled view resolves to x = 0)
%     3-view batch, 8x3 field matrix holding 10 / 20 / 30 per column
%         averages = 10 / 20 / 30      areas = 1 / 1 / 1
%     the same matrix sent WITHOUT the row-major transpose
%         averages = 20 / 18.33336 / 15    <- the bug this suite guards against
%     sliver thinner than one pixel
%         covered = 0, has_field = 1, has_stats = 0
%
%   See also FLATLAND, FLATLAND_LOAD.

    fprintf('\nFlatLand MATLAB binding self-check\n');
    fprintf('==================================\n');

    state = struct('passed', 0, 'failed', 0, 'messages', {{}});

    % A closed unit box spanning [0,1]^3 with CCW outward normals. 8 vertices
    % and 12 faces, so node and face fields are never ambiguous.
    % FACES ARE 1-BASED here, as this binding requires.
    BOX_V = [0 0 0; 1 0 0; 1 1 0; 0 1 0; ...
             0 0 1; 1 0 1; 1 1 1; 0 1 1];
    BOX_F = [1 4 3; 1 3 2;  5 6 7; 5 7 8; ...
             1 2 6; 1 6 5;  2 3 7; 2 7 6; ...
             3 4 8; 3 8 7;  4 1 5; 4 5 8];

    % Handle objects created below are released when this function's workspace
    % is destroyed, including on an early error, so there is no cleanup object
    % here: an anonymous function would capture the variable's value at
    % construction time (still empty) rather than the object assigned later.
    box = [];

    try
        % =================================================================
        section('library and version');

        info = flatland_load();
        state = check(state, ~isempty(info.alias), 'the library loads');
        state = check(state, exist(info.libfile, 'file') == 2 || strcmp(info.libfile, 'libflatland'), ...
                      sprintf('loaded %s', info.libfile));
        state = check(state, info.version(1) == 3, ...
                      sprintf('library reports version %d.%d.%d', info.version));
        state = check(state, ~isempty(Flatland.versionString()), ...
                      'fl_version_string is non-empty');
        state = check(state, strcmp(Flatland.statusString(0), 'ok'), ...
                      'fl_status_string(FL_OK) is "ok"');

        if ~isempty(info.warnings)
            fprintf('  note: loadlibrary emitted %d warning(s):\n', numel(info.warnings));
            for k = 1:min(numel(info.warnings), 5)
                fprintf('        %s\n', info.warnings{k});
            end
        end

        % The three structs that cross the boundary must all be registered,
        % or every marshalled call below is guesswork.
        for want = {'fl_options', 'fl_result', 'fl_batch_desc'}
            state = check(state, any(strcmp(info.structs, want{1})), ...
                          sprintf('loadlibrary registered struct %s', want{1}));
        end

        % =================================================================
        section('struct layout');

        % Field names and classes must match the C structs exactly, or MATLAB
        % marshals into the wrong offsets. Measured on the build host:
        %   fl_options 56 bytes, fl_result 96 bytes, fl_batch_desc 88 bytes.
        r0 = libstruct('fl_result');
        resultFields = {'area', 'average', 'integral', 'min', 'max', ...
                        'covered_pixels', 'width', 'height', 'has_field', 'has_stats'};
        missing = {};
        for k = 1:numel(resultFields)
            try
                v = r0.(resultFields{k});           %#ok<NASGU>
            catch
                missing{end+1} = resultFields{k};   %#ok<AGROW>
            end
        end
        state = check(state, isempty(missing), ...
                      sprintf('fl_result exposes all 10 fields%s', fmtMissing(missing)));

        o0 = libstruct('fl_options');
        calllib(info.alias, 'fl_options_init', o0);
        state = check(state, abs(double(o0.resolution) - 0.001) < 1e-15, ...
                      sprintf('fl_options_init default resolution = %g', double(o0.resolution)));
        state = check(state, double(o0.cull) == 1, 'fl_options_init default cull = 1');
        state = check(state, double(o0.threads) == 0, 'fl_options_init default threads = 0');
        state = check(state, all(double(o0.reserved) == 0), ...
                      'fl_options_init zeroes the reserved words');

        % =================================================================
        section('mesh from MATLAB arrays');

        box = Flatland(BOX_V, BOX_F);
        state = check(state, box.vertexCount == 8, 'vertex count round-trips');
        state = check(state, box.faceCount == 12, 'face count round-trips');

        % This is the real test of the column-major -> interleaved transpose:
        % V(:) instead of V.'(:) would come back permuted.
        V = box.vertices();
        state = check(state, isequal(size(V), [8 3]) && max(abs(V(:) - BOX_V(:))) < 1e-12, ...
                      'vertices() round-trips in the original frame and orientation');

        % And this is the real test of the 1-based <-> 0-based conversion.
        F = box.faces();
        state = check(state, isequal(size(F), [12 3]) && isequal(F, BOX_F), ...
                      'faces() round-trips, still 1-based');
        state = check(state, min(F(:)) == 1 && max(F(:)) == 8, ...
                      sprintf('face indices span 1..%d, not 0..%d', max(F(:)), max(F(:)) - 1));

        % =================================================================
        section('projection');

        r = box.project([1 0 0], 'Resolution', 0.002, 'Precision', 'double');
        state = near(state, r.area, 1.0, 0.01, 'unit box projects to area 1');
        state = check(state, r.coveredPixels == 250000, ...
                      sprintf('covered pixels = %d (C measured 250000)', r.coveredPixels));
        state = check(state, r.width == 502 && r.height == 502, ...
                      sprintf('raster is %dx%d (C measured 502x502)', r.width, r.height));
        state = check(state, ~r.hasField && ~r.hasStats, 'no field means no statistics');
        state = check(state, isnan(r.average) && isnan(r.integral) && ...
                             isnan(r.min) && isnan(r.max), ...
                      'a fieldless view reports NaN statistics, not zeros');

        % Field = vertex x. Looking along +X the visible surface is x = 0, so
        % the mean must be 0. This is what catches an inverted cull, and it
        % needs a field that VARIES across the mesh.
        fx = BOX_V(:, 1);
        r = box.project([1 0 0], 'Field', fx, 'Resolution', 0.002, 'Precision', 'double');
        state = check(state, r.hasField && r.hasStats, 'a covered field view has statistics');
        state = near(state, r.average, 0.0, 1e-6, 'culled view resolves to the NEAR surface');
        state = near(state, r.min, 0.0, 1e-6, 'field min over the near surface');
        state = near(state, r.max, 0.0, 1e-6, 'field max over the near surface');

        r2 = box.project([1 0 0], 'Field', fx, 'Resolution', 0.002, ...
                         'Precision', 'double', 'Cull', false);
        state = near(state, r2.average, 0.0, 1e-6, 'culling agrees with no culling');

        % A constant face field: the area integral is the constant times the area.
        ff = 4.0 * ones(12, 1);
        r = box.project([0 0 1], 'Field', ff, 'Resolution', 0.002, 'Precision', 'double');
        state = near(state, r.average, 4.0, 1e-9, 'constant face field average');
        state = near(state, r.integral, 4.0 * r.area, 1e-9, 'integral == constant * area');

        % Node vs face is chosen by length here, since 8 ~= 12.
        state = check(state, r.hasField, 'a 12-row field was taken as a face field');

        rf = box.project([1 0 0], 'Resolution', 0.002, 'Precision', 'single');
        rd = box.project([1 0 0], 'Resolution', 0.002, 'Precision', 'double');
        state = near(state, rf.area, rd.area, 1e-3, 'single and double agree on area');

        % =================================================================
        section('NaN, not zero, when has_stats is false');

        % A sliver thinner than one pixel covers nothing, so its field
        % statistics are not measurements. The C struct leaves them at zero;
        % this binding must report NaN. THIS IS THE MOST IMPORTANT CHECK HERE.
        sliver = Flatland([0 0 0; 1 0 0; 1 0.000001 0], [1 2 3]);
        cleanupSliver = onCleanup(@() deleteIfLive(sliver)); %#ok<NASGU>
        rs = sliver.project([0 0 -1], 'Field', [-5; -3; -1], ...
                            'Resolution', 0.1, 'Cull', false, 'Precision', 'double');

        state = check(state, rs.coveredPixels == 0, 'the sliver covers no pixels');
        state = check(state, rs.hasField, 'hasField is true (a field was supplied)');
        state = check(state, ~rs.hasStats, 'hasStats is false (nothing was measured)');
        state = check(state, isnan(rs.average),  'average  is NaN, not 0');
        state = check(state, isnan(rs.integral), 'integral is NaN, not 0');
        state = check(state, isnan(rs.min),      'min      is NaN, not 0');
        state = check(state, isnan(rs.max),      'max      is NaN, not 0');
        % The two values that ARE measurements must survive as real zeros.
        state = check(state, rs.area == 0 && ~isnan(rs.area), ...
                      'area stays a real measurement (0), not NaN');
        state = check(state, rs.coveredPixels == 0 && ~isnan(rs.coveredPixels), ...
                      'coveredPixels stays a real measurement (0), not NaN');

        % =================================================================
        section('batch');

        views = [1 0 0; 0 1 0; 0 0 1];
        % 8 rows (one per vertex) x 3 columns (one per timestep).
        M = repmat([10 20 30], 8, 1);

        R = box.projectBatch(views, 'FieldMatrix', M, ...
                             'FieldColumns', [1 2 3], ...
                             'Resolutions', [0.004 0.004 0.004], ...
                             'Precision', 'double');

        state = check(state, numel(R) == 3, 'the batch returns one result per view');
        state = near(state, R(1).average, 10.0, 1e-9, 'view 1 reads column 1');
        state = near(state, R(2).average, 20.0, 1e-9, 'view 2 reads column 2');
        state = near(state, R(3).average, 30.0, 1e-9, 'view 3 reads column 3');
        state = near(state, R(1).area, 1.0, 0.02, 'batch view 1 area');

        % Regression guard on the row-major/column-major transpose. Sending M
        % unchanged instead of M.' does not error, it just returns these
        % numbers (measured in C): 20 / 18.33336 / 15. If the transpose is
        % ever dropped, this is what the averages become.
        state = check(state, abs(R(1).average - 20.0) > 1e-6, ...
                      'field matrix was transposed to row-major (view 1 is not 20)');
        state = check(state, abs(R(2).average - 18.33336) > 1e-3, ...
                      'field matrix was transposed to row-major (view 2 is not 18.3)');
        state = check(state, abs(R(3).average - 15.0) > 1e-6, ...
                      'field matrix was transposed to row-major (view 3 is not 15)');

        % Omitting FieldColumns means every view reads column 1.
        R0 = box.projectBatch(views, 'FieldMatrix', M, ...
                              'Resolutions', [0.004 0.004 0.004], 'Precision', 'double');
        state = near(state, R0(3).average, 10.0, 1e-9, ...
                     'no FieldColumns means every view reads column 1');

        % Results must not depend on how many workers ran them.
        R1 = box.projectBatch(views, 'FieldMatrix', M, 'FieldColumns', [1 2 3], ...
                              'Resolutions', [0.004 0.004 0.004], ...
                              'Precision', 'double', 'Threads', 1);
        state = check(state, isequal([R.average], [R1.average]), ...
                      'batch results are independent of thread count');

        % Geometry-only batch.
        Rg = box.projectBatch(views, 'Resolution', 0.004, 'Precision', 'double');
        state = near(state, Rg(1).area, 1.0, 0.02, 'geometry-only batch still measures area');
        state = check(state, all(isnan([Rg.average])), ...
                      'a geometry-only batch reports NaN averages, not zeros');

        % A single view is accepted as a 1-by-3 row.
        Rs = box.projectBatch([1 0 0], 'Resolution', 0.004, 'Precision', 'double');
        state = check(state, numel(Rs) == 1, 'a single view row is a valid batch');

        % =================================================================
        section('render');

        fy  = BOX_V(:, 2);
        img = box.render([1 0 0], 'Field', fy, 'Resolution', 0.002, 'Precision', 'double');
        state = check(state, img.width == img.result.width && ...
                             img.height == img.result.height, ...
                      'image dimensions match the result');
        state = check(state, isequal(size(img.mask), [img.height img.width]), ...
                      sprintf('mask is height-by-width (%dx%d)', img.height, img.width));
        state = check(state, nnz(img.mask) == img.result.coveredPixels, ...
                      'mask agrees with coveredPixels');
        state = check(state, isequal(size(img.values), [img.height img.width]), ...
                      'values is height-by-width');
        state = check(state, all(isnan(img.values(~img.mask))), ...
                      'uncovered pixels in values are NaN');
        state = check(state, ~any(isnan(img.values(img.mask))), ...
                      'covered pixels in values are real numbers');

        imgNoField = box.render([1 0 0], 'Resolution', 0.004, 'Precision', 'double');
        state = check(state, isempty(imgNoField.values), ...
                      'a fieldless render has empty values');
        state = check(state, nnz(imgNoField.mask) == imgNoField.result.coveredPixels, ...
                      'a fieldless render still has a mask');

        % A view that covers nothing must yield a valid, empty image.
        tri = Flatland([0 0 0; 1 0 0; 1 1 0], [1 2 3]);
        cleanupTri = onCleanup(@() deleteIfLive(tri)); %#ok<NASGU>
        empty = tri.render([1 0 0], 'Resolution', 0.002, 'Precision', 'double');
        state = check(state, empty.width == 0 && empty.height == 0, ...
                      'an edge-on view yields a 0x0 image, not a stale one');
        state = check(state, isempty(empty.mask), 'the empty image has an empty mask');
        state = check(state, empty.result.coveredPixels == 0 && ~empty.result.hasStats, ...
                      'the empty view reports no coverage and no statistics');

        % =================================================================
        section('angle conversion');

        d = Flatland.angleToDir(0, 0);
        state = near(state, norm(d - [1 0 0]), 0, 1e-12, 'angleToDir(0,0) is +X');
        d = Flatland.angleToDir(90, 0);
        state = near(state, norm(d - [0 1 0]), 0, 1e-12, 'angleToDir(90,0) is +Y');
        d = Flatland.angleToDir(0, 90);
        state = near(state, norm(d - [0 0 1]), 0, 1e-12, 'angleToDir(0,90) is +Z');
        D = Flatland.angleToDir([0 90 0], [0 0 90]);
        state = check(state, isequal(size(D), [3 3]), 'angleToDir vectorises to N-by-3');

        % =================================================================
        section('input validation');

        state = throwsId(state, @() box.project([1 0 0], 'Field', zeros(7, 1)), ...
                         'Flatland:InvalidField', 'a wrong-length field is rejected');
        state = throwsId(state, @() box.project([0 0 0]), ...
                         'Flatland:InvalidView', 'a zero view direction is rejected');
        state = throwsId(state, @() box.project([1 0]), ...
                         'Flatland:InvalidView', 'a 2-element view is rejected');
        state = throwsId(state, @() box.project([1 0 0], 'Nonsense', 1), ...
                         'Flatland:InvalidOption', 'an unknown option is rejected');
        state = throwsId(state, @() box.project([1 0 0], 'Resolution'), ...
                         'Flatland:InvalidOption', 'a dangling option name is rejected');
        state = throwsId(state, @() box.project([1 0 0], 'Resolution', 0), ...
                         'Flatland:InvalidOption', 'a zero resolution is rejected');
        state = throwsId(state, @() box.project([1 0 0], 'Resolution', -1), ...
                         'Flatland:InvalidOption', 'a negative resolution is rejected');
        state = throwsId(state, @() box.project([1 0 0], 'Precision', 'quad'), ...
                         'Flatland:InvalidOption', 'an unknown precision is rejected');
        state = throwsId(state, @() box.project([1 0 0], 'Field', [NaN; zeros(7,1)]), ...
                         'Flatland:InvalidField', 'a non-finite field is rejected');
        state = throwsId(state, @() Flatland(BOX_V, BOX_F - 1), ...
                         'Flatland:InvalidFaces', '0-based faces are rejected with advice');
        state = throwsId(state, @() Flatland(BOX_V, [1 2 99]), ...
                         'Flatland:InvalidFaces', 'an out-of-range face index is rejected');
        state = throwsId(state, @() Flatland(BOX_V.', BOX_F), ...
                         'Flatland:InvalidVertices', '3-by-N vertices are rejected with advice');
        state = throwsId(state, @() Flatland.load('no_such_file_here.obj'), ...
                         'Flatland:FileNotFound', 'a missing mesh file is rejected');
        state = throwsId(state, @() box.projectBatch(views, 'FieldMatrix', M, ...
                                                     'FieldColumns', [1 2 9]), ...
                         'Flatland:InvalidField', 'an out-of-range field column is rejected');
        state = throwsId(state, @() box.projectBatch(views, 'FieldColumns', [1 2 3]), ...
                         'Flatland:InvalidField', 'FieldColumns without FieldMatrix is rejected');
        state = throwsId(state, @() box.projectBatch(views, 'FieldMatrix', M, ...
                                                     'FieldColumns', [1 2]), ...
                         'Flatland:InvalidField', 'a short FieldColumns is rejected');

        % =================================================================
        section('errors from the library itself');

        % An impossibly fine resolution is refused by the engine, not by the
        % MATLAB-side validation, so this exercises the fl_status path and
        % checks that fl_last_error text reaches the MATLAB message.
        msg = '';
        try
            box.project([1 0 0], 'Resolution', 1e-9);
            state = fail(state, 'an impossibly fine resolution is refused');
        catch err
            msg = err.message;
            state = check(state, ~isempty(strfind(err.identifier, 'Flatland:')), ...
                          sprintf('library error uses a Flatland identifier (%s)', err.identifier));
            state = check(state, ~isempty(strfind(msg, 'fl_status')), ...
                          'the message names the fl_status code');
            state = check(state, length(msg) > 40, ...
                          'the message carries fl_last_error detail');
        end
        if ~isempty(msg)
            fprintf('       -> %s\n', truncate(msg, 130));
        end

        % =================================================================
        section('lifetime');

        tmp = Flatland(BOX_V, BOX_F);
        state = check(state, tmp.hasMesh(), 'a fresh object owns a mesh');
        delete(tmp);
        state = check(state, ~tmp.hasMesh(), 'delete() releases the handle');
        state = throwsId(state, @() tmp.project([1 0 0]), ...
                         'Flatland:NoMesh', 'using a deleted object is rejected');
        delete(tmp);   % must be safe twice
        state = check(state, true, 'delete() is safe to call twice');

        % =================================================================
        section('mesh from file');

        objPath = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                           'examples', 'cube_area', 'cube.obj');
        if exist(objPath, 'file') == 2
            cube = Flatland.load(objPath);
            cleanupCube = onCleanup(@() deleteIfLive(cube)); %#ok<NASGU>
            state = check(state, cube.vertexCount == 26 && cube.faceCount == 48, ...
                          sprintf('cube.obj loads (%d vertices, %d faces; C measured 26/48)', ...
                                  cube.vertexCount, cube.faceCount));
            rc = cube.project([1 0 0], 'Resolution', 0.005, 'Precision', 'double');
            state = near(state, rc.area, 1.0, 0.02, 'the cube example projects to area 1');
            Fc = cube.faces();
            state = check(state, min(Fc(:)) >= 1, ...
                          'faces from an OBJ come back 1-based');
        else
            fprintf('  (skipped: %s not found)\n', objPath);
        end

    catch err
        state.failed = state.failed + 1;
        fprintf(2, '\n  ABORTED: %s\n       %s\n', err.identifier, err.message);
        for k = 1:min(numel(err.stack), 6)
            fprintf(2, '       at %s line %d\n', err.stack(k).name, err.stack(k).line);
        end
    end

    deleteIfLive(box);

    fprintf('\n----------------------------------\n');
    fprintf('MATLAB binding: %d passed, %d failed\n\n', state.passed, state.failed);

    if nargout > 0
        failures = state.failed;
    elseif state.failed > 0
        error('Flatland:TestFailed', '%d check(s) failed.', state.failed);
    end
end

% =========================================================================
% helpers

function section(name)
    fprintf('\n== %s ==\n', name);
end

function s = check(s, cond, what)
    if cond
        fprintf('  PASS %s\n', what);
        s.passed = s.passed + 1;
    else
        fprintf(2, '  FAIL %s\n', what);
        s.failed = s.failed + 1;
    end
end

function s = fail(s, what)
    s = check(s, false, what);
end

function s = near(s, got, want, tol, what)
    if ~isscalar(got) || ~isnumeric(got)
        s = check(s, false, sprintf('%s (not a numeric scalar)', what));
        return;
    end
    if isnan(got) && ~isnan(want)
        s = check(s, false, sprintf('%s (got NaN, want %.9g)', what, want));
        return;
    end
    ok = abs(got - want) <= tol;
    if ok
        s = check(s, true, sprintf('%s (%.9g ~= %.9g)', what, got, want));
    else
        s = check(s, false, sprintf('%s (got %.9g, want %.9g)', what, got, want));
    end
end

function s = throwsId(s, fn, wantId, what)
    try
        fn();
        s = check(s, false, sprintf('%s (no error was raised)', what));
    catch err
        if strcmp(err.identifier, wantId)
            s = check(s, true, sprintf('%s [%s]', what, wantId));
        else
            s = check(s, false, sprintf('%s (got %s, want %s)', ...
                                        what, err.identifier, wantId));
        end
    end
end

function deleteIfLive(obj)
    try
        if ~isempty(obj) && isa(obj, 'Flatland') && isvalid(obj)
            delete(obj);
        end
    catch
    end
end

function t = fmtMissing(missing)
    if isempty(missing)
        t = '';
    else
        t = sprintf(' (missing: %s)', strjoin(missing, ', '));
    end
end

function t = truncate(s, n)
    s = s(:).';
    if length(s) > n
        t = [s(1:n) '...'];
    else
        t = s;
    end
end
