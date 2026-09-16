classdef Flatland < handle
%FLATLAND  Projected area and field statistics for a triangle mesh.
%
%   FlatLand computes the visible (projected) surface area of a triangle mesh
%   from an arbitrary view direction, and, given a scalar field, the mean,
%   extremes and area integral of that field over the visible projection.
%
%   This class is a thin MATLAB wrapper over the FlatLand C ABI, reached with
%   loadlibrary/calllib. It owns a mesh handle and frees it in delete().
%
%   CONSTRUCTION
%     fl = Flatland(vertices, faces)   vertices N-by-3 double, faces M-by-3
%     fl = Flatland.load('part.obj')   OBJ or STL, ASCII or binary
%
%   USE
%     r   = fl.project([1 0 0])
%     r   = fl.project([1 0 0], 'Field', values, 'Resolution', 1e-3)
%     R   = fl.projectBatch(views, 'FieldMatrix', M, 'FieldColumns', cols)
%     img = fl.render([1 0 0], 'Field', values)
%     delete(fl)
%
%   ============================ 1-BASED INDEXING ==========================
%   Everything this class takes and returns is 1-BASED, the MATLAB way:
%
%     * FACES are 1-based. faces(k,:) = [1 2 3] refers to vertices(1:3,:).
%       The underlying C API is 0-based throughout; this class subtracts one
%       on the way in and adds one back in faces(). Do NOT pre-subtract.
%     * FIELDCOLUMNS are 1-based. 'FieldColumns', [1 2 3] reads the first,
%       second and third columns of the field matrix.
%     * IMAGE rows are 1-based, and row 1 is the BOTTOM row in mesh space
%       (the C API's row 0). Use flipud() for the usual top-down orientation.
%
%   This is the single easiest thing to get wrong when moving between the C
%   API, the CLI and this binding. If your projected areas look plausible but
%   your field values look scrambled, check this first.
%   ========================================================================
%
%   ============================ NaN, NOT ZERO =============================
%   A view can legitimately cover no pixels: geometry thinner than one pixel,
%   or a mesh seen exactly edge-on. The C struct then reports has_stats = 0
%   and leaves average/integral/min/max at ZERO. Zero is a perfectly plausible
%   field value, so propagating it would be indistinguishable from a real
%   measurement. This class reports those four as NaN whenever hasStats is
%   false. A view with no coverage has no statistics; it does not have
%   statistics equal to zero. area and coveredPixels stay real measurements
%   (both 0) and are never NaN'd.
%   ========================================================================
%
%   MATRIX LAYOUT
%   MATLAB is column-major and the C API reads its arrays row-major, so this
%   class transposes on the way through. Pass the MATLAB-natural shapes --
%   vertices N-by-3, faces M-by-3, views K-by-3, field matrix
%   nEntities-by-nTimesteps -- and it handles the rest.
%
%   See also FLATLAND_LOAD, FLATLAND_TEST.

    properties (SetAccess = private)
        % Plain defaults rather than typed property declarations, so the class
        % loads on releases predating property validation syntax.
        vertexCount = 0           % number of vertices in the mesh
        faceCount   = 0           % number of triangles in the mesh
    end

    properties (Access = private)
        handle = []               % lib.pointer to the C fl_mesh, or [] if none
        info   = []               % cached struct from flatland_load()
    end

    properties (Constant, Access = private)
        % Measured with sizeof() on the build host (x86-64 Linux): fl_result is
        % 96 bytes. Used only by the byte-buffer fallback in allocResults.
        RESULT_BYTES = 96
    end

    % =====================================================================
    methods
        % -----------------------------------------------------------------
        function obj = Flatland(vertices, faces)
        %FLATLAND  Build a mesh from vertex and face arrays.
        %
        %   FL = FLATLAND(VERTICES, FACES) where VERTICES is N-by-3 double
        %   (x y z per row) and FACES is M-by-3 of 1-BASED vertex indices.
        %
        %   Both arrays are copied by the library, so you may modify or clear
        %   your own immediately afterwards.
        %
        %   FL = FLATLAND() returns an object with no mesh attached. Only
        %   Flatland.load uses this; every method rejects it until a mesh
        %   exists.

            obj.info = flatland_load();

            if nargin == 0
                return;
            end
            if nargin ~= 2
                error('Flatland:NotEnoughInputs', ...
                      ['Flatland requires both vertices and faces. ' ...
                       'Use Flatland.load(path) to read a mesh file.']);
            end

            vertices = Flatland.validateVertices(vertices);
            faces    = Flatland.validateFaces(faces, size(vertices, 1));

            nv = size(vertices, 1);
            nf = size(faces, 1);

            % MATLAB stores column-major but the C API wants xyz interleaved,
            % so transpose first: V.'(:) gives x1 y1 z1 x2 y2 z2 ... whereas
            % V(:) would give every x, then every y, then every z.
            vflat = vertices.';
            fflat = int32(faces.' - 1);        % ... and 0-based for C

            pv = libpointer('doublePtr', vflat(:));
            pf = libpointer('int32Ptr',  fflat(:));

            pp  = Flatland.newOutHandle(obj.info.meshPtrType);
            out = obj.raw('fl_mesh_create', pv, uint64(nv), pf, uint64(nf), pp);
            [status, h] = Flatland.takeHandle(out, pp, 'fl_mesh_create');
            obj.check(status, 'fl_mesh_create');

            obj.adoptHandle(h);
        end

        % -----------------------------------------------------------------
        function delete(obj)
        %DELETE  Release the underlying C mesh handle.
        %
        %   Called automatically when the object goes out of scope, or
        %   explicitly with delete(fl). Safe to call more than once.

            if ~isempty(obj.handle)
                h = obj.handle;
                obj.handle = [];
                try
                    obj.raw('fl_mesh_destroy', h);
                catch
                    % A destructor must never throw. The usual cause is the
                    % library having been unloaded already, in which case there
                    % is nothing left to free.
                end
            end
        end

        % -----------------------------------------------------------------
        function tf = hasMesh(obj)
        %HASMESH  True while this object still owns a live mesh handle.
            tf = ~isempty(obj.handle);
        end

        % -----------------------------------------------------------------
        function r = project(obj, view, varargin)
        %PROJECT  Project one view and return area and field statistics.
        %
        %   R = FL.PROJECT(VIEW) projects along the 3-vector VIEW, which is the
        %   direction the camera LOOKS ALONG; the visible surface is the one
        %   facing back toward -VIEW. Magnitude is irrelevant, the vector is
        %   normalised internally.
        %
        %   R = FL.PROJECT(VIEW, 'Name', value, ...) accepts:
        %     'Field'       vector of field values, one per vertex (node field)
        %                   or one per face (face field). Default: none.
        %     'Resolution'  pixel edge length in mesh units. Default 1e-3.
        %     'Precision'   'single' (or 'float') | 'double'. Default 'single'.
        %                   Selects the engine's internal working type only;
        %                   every value crossing the boundary is double either
        %                   way.
        %     'Cull'        true (default) to backface cull, false to render
        %                   every face.
        %     'FieldMode'   'auto' (default) | 'node' | 'face'. Needed only when
        %                   the mesh has equally many vertices and faces, which
        %                   makes the field length ambiguous.
        %
        %   R has fields:
        %     area           covered area, in mesh units squared
        %     average        mean field value over covered pixels, or NaN
        %     integral       integral of f dA over the projection, or NaN
        %     min, max       field extremes over covered pixels, or NaN
        %     coveredPixels  number of covered pixels
        %     width, height  raster dimensions actually used
        %     hasField       true if a field was supplied for this view
        %     hasStats       true if the four statistics are real measurements
        %
        %   When hasStats is false, average/integral/min/max are NaN rather
        %   than 0. See the class help.

            obj.assertLive();
            view = Flatland.validateView(view);
            o    = obj.parseOptions(varargin, {'Field'});

            [pField, fieldLen] = obj.prepareField(o.Field, o.FieldMode);
            popts = obj.buildOptions(o);
            pres  = libstruct('fl_result');

            out = obj.raw('fl_project', obj.handle, view, pField, ...
                          uint64(fieldLen), popts, pres);
            obj.check(out{1}, 'fl_project');

            r = Flatland.decodeResult(pres);
        end

        % -----------------------------------------------------------------
        function R = projectBatch(obj, views, varargin)
        %PROJECTBATCH  Project many views in one parallel call.
        %
        %   R = FL.PROJECTBATCH(VIEWS) where VIEWS is K-by-3, one view
        %   direction per row. Returns a 1-by-K struct array with the same
        %   fields PROJECT returns. Gather a column with [R.area].
        %
        %   R = FL.PROJECTBATCH(VIEWS, 'Name', value, ...) accepts everything
        %   PROJECT does except 'Field', plus:
        %     'FieldMatrix'   nEntities-by-nTimesteps double. One ROW per mesh
        %                     entity (vertex or face), one COLUMN per timestep
        %                     -- the MATLAB-natural orientation.
        %     'FieldColumns'  K 1-BASED column indices, one per view, saying
        %                     which column of FieldMatrix each view reads.
        %                     Default: every view reads column 1.
        %     'Resolutions'   K per-view pixel sizes, overriding 'Resolution'.
        %     'Threads'       worker threads; 0 (default) means one per core.
        %
        %   The whole view array and the whole field matrix cross the boundary
        %   in a SINGLE calllib. Looping per view from MATLAB would re-validate
        %   and re-centre the mesh every time and throw away the parallelism
        %   that makes this path worth having.
        %
        %   ROW-MAJOR vs COLUMN-MAJOR: the C API reads FieldMatrix row-major,
        %   element (r,c) at data[r*nCols + c], and MATLAB stores column-major,
        %   so this method sends M.' and lets MATLAB's column-major flattening
        %   produce exactly C's row-major order. Sending M unchanged does not
        %   error -- it silently scrambles which entity holds which value.

            obj.assertLive();
            views = Flatland.validateViews(views);
            K     = size(views, 1);

            o = obj.parseOptions(varargin, ...
                    {'FieldMatrix', 'FieldColumns', 'Resolutions', 'Threads'});

            popts = obj.buildOptions(o);
            pdesc = libstruct('fl_batch_desc');
            obj.raw('fl_batch_desc_init', pdesc);   % zeroes the reserved words

            % Every libpointer below must stay in scope until after the call:
            % fl_batch_desc stores raw addresses into their buffers, and MATLAB
            % may release a libpointer that nothing references any more.
            vflat  = views.';
            pViews = libpointer('doublePtr', vflat(:));
            pdesc.views      = pViews;
            pdesc.view_count = uint64(K);

            pMat = []; pCols = []; pRes = [];

            if ~isempty(o.FieldMatrix)
                M = o.FieldMatrix;
                if ~isnumeric(M) || ~ismatrix(M) || isempty(M)
                    error('Flatland:InvalidField', ...
                          '''FieldMatrix'' must be a non-empty numeric matrix.');
                end
                if ~all(isfinite(M(:)))
                    error('Flatland:InvalidField', ...
                          '''FieldMatrix'' contains %d non-finite value(s).', ...
                          sum(~isfinite(M(:))));
                end
                rows = size(M, 1);
                cols = size(M, 2);
                obj.assertFieldLength(rows, o.FieldMode, 'FieldMatrix');

                Mt   = double(M).';          % the column-major -> row-major step
                pMat = libpointer('doublePtr', Mt(:));
                pdesc.field_matrix = pMat;
                pdesc.field_rows   = uint64(rows);
                pdesc.field_cols   = uint64(cols);

                if ~isempty(o.FieldColumns)
                    c = o.FieldColumns(:);
                    if ~isnumeric(c) || numel(c) ~= K
                        error('Flatland:InvalidField', ...
                              '''FieldColumns'' needs one entry per view (%d), got %d.', ...
                              K, numel(c));
                    end
                    if any(c ~= fix(c)) || any(c < 1) || any(c > cols)
                        error('Flatland:InvalidField', ...
                              ['''FieldColumns'' entries must be integers in 1..%d ' ...
                               '(1-based); got values spanning %g..%g.'], ...
                              cols, min(c), max(c));
                    end
                    pCols = libpointer('int32Ptr', int32(c - 1));   % C is 0-based
                    pdesc.field_columns = pCols;
                end
            elseif ~isempty(o.FieldColumns)
                error('Flatland:InvalidField', ...
                      '''FieldColumns'' was given without ''FieldMatrix''.');
            end

            if ~isempty(o.Resolutions)
                res = o.Resolutions(:);
                if ~isnumeric(res) || numel(res) ~= K
                    error('Flatland:InvalidOption', ...
                          '''Resolutions'' needs one entry per view (%d), got %d.', ...
                          K, numel(res));
                end
                if any(~isfinite(res)) || any(res <= 0)
                    error('Flatland:InvalidOption', ...
                          '''Resolutions'' entries must be positive and finite.');
                end
                pRes = libpointer('doublePtr', double(res));
                pdesc.resolutions = pRes;
            end

            pOut = Flatland.allocResults(K);

            out = obj.raw('fl_project_batch', obj.handle, pdesc, popts, pOut);
            obj.check(out{1}, 'fl_project_batch');

            R = Flatland.decodeResultArray(pOut.Value, K);

            % pViews/pMat/pCols/pRes are deliberately still in scope here.
            % fl_batch_desc holds their raw addresses, and MATLAB does not see
            % that as a reference, so the named locals are what keep the
            % buffers alive across the call.
        end

        % -----------------------------------------------------------------
        function img = render(obj, view, varargin)
        %RENDER  Project one view and keep the raster.
        %
        %   IMG = FL.RENDER(VIEW, ...) takes the same options as PROJECT and
        %   returns a struct with:
        %     width, height  raster size (both 0 if the view covers nothing)
        %     mask           height-by-width logical, true where covered
        %     values         height-by-width double of interpolated field
        %                    values, NaN at uncovered pixels; [] when no field
        %                    was supplied
        %     result         the same struct PROJECT returns
        %
        %   ROW 1 IS THE BOTTOM ROW in mesh space, matching the C API. For the
        %   conventional top-down image orientation:
        %       imagesc(flipud(img.values)); axis image off
        %
        %   The library exposes its raster as borrowed memory valid only until
        %   the image handle is destroyed, so this method copies out of it and
        %   destroys the handle before returning. IMG is plain MATLAB data.

            obj.assertLive();
            view = Flatland.validateView(view);
            o    = obj.parseOptions(varargin, {'Field'});

            [pField, fieldLen] = obj.prepareField(o.Field, o.FieldMode);
            popts = obj.buildOptions(o);
            pres  = libstruct('fl_result');

            % NOTE the argument order: out_image comes BEFORE out_result in
            % fl_project_image, so the out-handle is second to last.
            pp  = Flatland.newOutHandle(obj.info.imagePtrType);
            out = obj.raw('fl_project_image', obj.handle, view, pField, ...
                          uint64(fieldLen), popts, pp, pres);
            [status, hImg] = Flatland.takeHandle(out, pp, 'fl_project_image');
            obj.check(status, 'fl_project_image');

            % Frees the image however this function exits, including on error.
            guard = onCleanup(@() obj.destroyImage(hImg)); %#ok<NASGU>

            wc = obj.raw('fl_image_width',  hImg);
            hc = obj.raw('fl_image_height', hImg);
            w  = double(wc{1});
            h  = double(hc{1});

            img = struct('width', w, 'height', h, ...
                         'mask', false(0, 0), 'values', [], ...
                         'result', Flatland.decodeResult(pres));

            if w <= 0 || h <= 0
                % An empty view yields a valid 0-by-0 image, not a stale raster.
                return;
            end

            n = w * h;

            mc = obj.raw('fl_image_mask', hImg);
            pm = mc{1};
            if ~Flatland.isNullPtr(pm)
                % The C function hands back a bare pointer with no length, so
                % tell MATLAB how much to read before touching .Value.
                setdatatype(pm, 'uint8Ptr', 1, n);
                % C order is row-major with row 0 at the bottom, indexed
                % [y*width + x]. reshape(...,[w h]).' recovers a
                % height-by-width MATLAB matrix whose row 1 is that bottom row.
                img.mask = reshape(double(pm.Value), w, h).' ~= 0;
            else
                img.mask = false(h, w);
            end

            vc = obj.raw('fl_image_values', hImg);
            pv = vc{1};
            if ~Flatland.isNullPtr(pv)
                setdatatype(pv, 'doublePtr', 1, n);
                vals = reshape(double(pv.Value), w, h).';
                % Uncovered pixels hold no measured value; say so with NaN
                % rather than passing on whatever is in the buffer.
                vals(~img.mask) = NaN;
                img.values = vals;
            end
        end

        % -----------------------------------------------------------------
        function V = vertices(obj)
        %VERTICES  Copy the mesh vertices back out as an N-by-3 matrix.
        %
        %   Returned in the original input frame. (FlatLand recentres geometry
        %   internally for numerical conditioning; that is not visible here.)

            obj.assertLive();
            n  = obj.vertexCount;
            pv = libpointer('doublePtr', zeros(n * 3, 1));
            out = obj.raw('fl_mesh_copy_vertices', obj.handle, pv, uint64(n * 3));
            obj.check(out{1}, 'fl_mesh_copy_vertices');
            V = reshape(double(pv.Value), 3, n).';
        end

        % -----------------------------------------------------------------
        function F = faces(obj)
        %FACES  Copy the mesh faces back out as an M-by-3 matrix.
        %
        %   Indices are 1-BASED, converted up from the C API's 0-based storage.

            obj.assertLive();
            m  = obj.faceCount;
            pf = libpointer('int32Ptr', zeros(m * 3, 1, 'int32'));
            out = obj.raw('fl_mesh_copy_faces', obj.handle, pf, uint64(m * 3));
            obj.check(out{1}, 'fl_mesh_copy_faces');
            F = reshape(double(pf.Value), 3, m).' + 1;
        end
    end

    % =====================================================================
    methods (Static)
        % -----------------------------------------------------------------
        function fl = load(path)
        %FLATLAND.LOAD  Build a mesh from an OBJ or STL file.
        %
        %   FL = FLATLAND.LOAD('part.obj') reads OBJ, or STL in either ASCII or
        %   binary form (auto-detected). Polygonal OBJ faces are fan-triangulated.
        %
        %   OBJ files on disk are 1-based and the loader converts; faces() then
        %   gives them back 1-based too, so the round trip is consistent.
        %
        %   STL carries no shared vertices, so a mesh loaded from STL has three
        %   vertices per triangle and its natural field mode is 'face'.

            if nargin < 1 || ~(ischar(path) || isstring(path))
                error('Flatland:InvalidPath', ...
                      'Flatland.load requires a file path.');
            end
            path = char(path);
            if isempty(path)
                error('Flatland:InvalidPath', ...
                      'Flatland.load requires a non-empty file path.');
            end
            % Checked here so the message names the path the user typed,
            % rather than surfacing a generic IO error from the library.
            if exist(path, 'file') ~= 2
                error('Flatland:FileNotFound', 'No such mesh file: %s', path);
            end

            fl = Flatland();                       % shell with no mesh yet
            info = flatland_load();
            pp   = Flatland.newOutHandle(info.meshPtrType);
            out  = Flatland.rawStatic(info, 'fl_mesh_load', path, pp);
            [status, h] = Flatland.takeHandle(out, pp, 'fl_mesh_load');
            fl.check(status, 'fl_mesh_load');
            fl.adoptHandle(h);
        end

        % -----------------------------------------------------------------
        function v = version()
        %FLATLAND.VERSION  Version of the loaded library, [major minor patch].
            info = flatland_load();
            v = info.version;
        end

        % -----------------------------------------------------------------
        function s = versionString()
        %FLATLAND.VERSIONSTRING  Version of the loaded library, as text.
            info = flatland_load();
            out  = Flatland.rawStatic(info, 'fl_version_string');
            s    = Flatland.toChar(out{1});
        end

        % -----------------------------------------------------------------
        function d = angleToDir(azimuthDeg, elevationDeg)
        %FLATLAND.ANGLETODIR  Azimuth/elevation in degrees to a unit direction.
        %
        %   D = FLATLAND.ANGLETODIR(AZ, EL) returns a 1-by-3 unit vector.
        %   Azimuth sweeps around +Z measured from +X; elevation rises from the
        %   XY plane toward +Z. AZ and EL may be equal-length vectors (or one
        %   of them scalar), in which case D is N-by-3 and feeds straight into
        %   projectBatch.

            if nargin < 2
                error('Flatland:NotEnoughInputs', ...
                      'Flatland.angleToDir requires azimuth and elevation, in degrees.');
            end
            if ~isnumeric(azimuthDeg) || ~isnumeric(elevationDeg)
                error('Flatland:InvalidView', 'azimuth and elevation must be numeric.');
            end
            az = double(azimuthDeg(:));
            el = double(elevationDeg(:));
            if isscalar(az) && ~isscalar(el), az = repmat(az, size(el)); end
            if isscalar(el) && ~isscalar(az), el = repmat(el, size(az)); end
            if numel(az) ~= numel(el)
                error('Flatland:InvalidView', ...
                      'azimuth and elevation must have the same number of elements.');
            end
            if ~all(isfinite(az)) || ~all(isfinite(el))
                error('Flatland:InvalidView', 'azimuth and elevation must be finite.');
            end

            info = flatland_load();
            d = zeros(numel(az), 3);
            for k = 1:numel(az)
                pd = libpointer('doublePtr', zeros(3, 1));
                Flatland.rawStatic(info, 'fl_angle_to_dir', az(k), el(k), pd);
                d(k, :) = double(pd.Value).';
            end
        end

        % -----------------------------------------------------------------
        function s = statusString(code)
        %FLATLAND.STATUSSTRING  Human-readable name of an fl_status code.
            info = flatland_load();
            out  = Flatland.rawStatic(info, 'fl_status_string', int32(code));
            s    = Flatland.toChar(out{1});
        end
    end

    % =====================================================================
    methods (Access = private)

        % -----------------------------------------------------------------
        function adoptHandle(obj, h)
        %ADOPTHANDLE  Take ownership of a fresh mesh handle and read its counts.
            obj.handle = h;
            vc = obj.raw('fl_mesh_vertex_count', obj.handle);
            fc = obj.raw('fl_mesh_face_count',   obj.handle);
            obj.vertexCount = double(vc{1});
            obj.faceCount   = double(fc{1});
        end

        % -----------------------------------------------------------------
        function destroyImage(obj, hImg)
            try
                obj.raw('fl_image_destroy', hImg);
            catch
            end
        end

        % -----------------------------------------------------------------
        function out = raw(obj, fname, varargin)
        %RAW  calllib with the right number of outputs, returned as a cell.
            out = Flatland.rawStatic(obj.info, fname, varargin{:});
        end

        % -----------------------------------------------------------------
        function assertLive(obj)
            if isempty(obj.handle)
                error('Flatland:NoMesh', ...
                      ['This Flatland object has no mesh. It was either deleted ' ...
                       'with delete(), or constructed with no arguments.']);
            end
        end

        % -----------------------------------------------------------------
        function check(obj, status, what)
        %CHECK  Turn a non-zero fl_status into a MATLAB error.
        %
        %   Every failing call leaves detail behind in fl_last_error(), which is
        %   thread-local and stays valid only until the next FlatLand call, so
        %   it is read immediately and folded into the message.

            code = Flatland.statusCode(status);
            if code == 0
                return;
            end

            name = '';
            detail = '';
            try
                sn   = obj.raw('fl_status_string', int32(code));
                name = Flatland.toChar(sn{1});
            catch
            end
            try
                le     = obj.raw('fl_last_error');
                detail = Flatland.toChar(le{1});
            catch
            end
            if isempty(name)
                name = sprintf('status %d', code);
            end

            id = Flatland.errorId(code);
            if isempty(detail)
                error(id, '%s failed: %s (fl_status %d).', what, name, code);
            else
                error(id, '%s failed: %s (fl_status %d): %s', what, name, code, detail);
            end
        end

        % -----------------------------------------------------------------
        function o = parseOptions(~, args, extraNames)
        %PARSEOPTIONS  Name/value pairs shared by project, projectBatch, render.

            o = struct('Field', [], 'FieldMatrix', [], 'FieldColumns', [], ...
                       'Resolution', 1e-3, 'Resolutions', [], ...
                       'Precision', 'single', 'Cull', true, ...
                       'FieldMode', 'auto', 'Threads', 0);

            allowed = [{'Resolution', 'Precision', 'Cull', 'FieldMode'}, extraNames];

            if mod(numel(args), 2) ~= 0
                error('Flatland:InvalidOption', ...
                      ['Options must be given as name/value pairs; got %d ' ...
                       'argument(s) after the view.'], numel(args));
            end

            for k = 1:2:numel(args)
                name = args{k};
                if ~(ischar(name) || isstring(name))
                    error('Flatland:InvalidOption', ...
                          'Expected an option name at argument %d, got a %s.', ...
                          k, class(name));
                end
                name = char(name);
                idx = find(strcmpi(allowed, name), 1);
                if isempty(idx)
                    error('Flatland:InvalidOption', ...
                          'Unknown option "%s". Accepted here: %s.', ...
                          name, strjoin(allowed, ', '));
                end
                o.(allowed{idx}) = args{k+1};
            end

            if ~isnumeric(o.Resolution) || ~isscalar(o.Resolution) || ...
               ~isfinite(o.Resolution) || o.Resolution <= 0
                error('Flatland:InvalidOption', ...
                      '''Resolution'' must be a positive finite scalar; got %s.', ...
                      Flatland.describe(o.Resolution));
            end

            o.Precision = Flatland.mapEnum(o.Precision, ...
                {'single', 'float', 'double'}, [0 0 1], 'Precision');

            o.FieldMode = Flatland.mapEnum(o.FieldMode, ...
                {'auto', 'node', 'vertex', 'face'}, [0 1 1 2], 'FieldMode');

            if ~(islogical(o.Cull) || isnumeric(o.Cull)) || ~isscalar(o.Cull)
                error('Flatland:InvalidOption', ...
                      '''Cull'' must be a logical scalar (true or false).');
            end
            o.Cull = logical(o.Cull);

            if ~isnumeric(o.Threads) || ~isscalar(o.Threads) || ...
               ~isfinite(o.Threads) || o.Threads < 0 || o.Threads ~= fix(o.Threads)
                error('Flatland:InvalidOption', ...
                      '''Threads'' must be a non-negative integer (0 = one per core).');
            end
        end

        % -----------------------------------------------------------------
        function popts = buildOptions(obj, o)
        %BUILDOPTIONS  An fl_options the library initialised, then overridden.
        %
        %   The header is explicit that fl_options_init must be called before
        %   any field is set, because it zeroes the reserved words that future
        %   versions depend on for ABI compatibility. So the struct is always
        %   initialised by the library rather than assembled here from scratch.

            popts = libstruct('fl_options');
            obj.raw('fl_options_init', popts);

            popts.resolution = double(o.Resolution);
            popts.cull       = int32(o.Cull);
            popts.threads    = int32(o.Threads);

            Flatland.setEnum(popts, 'precision',  o.Precision, ...
                             {'FL_PRECISION_FLOAT', 'FL_PRECISION_DOUBLE'});
            Flatland.setEnum(popts, 'field_mode', o.FieldMode, ...
                             {'FL_FIELD_AUTO', 'FL_FIELD_NODE', 'FL_FIELD_FACE'});
        end

        % -----------------------------------------------------------------
        function [pField, fieldLen] = prepareField(obj, field, fieldMode)
        %PREPAREFIELD  Validate a per-view field and wrap it for the call.

            if isempty(field)
                pField   = libpointer('doublePtr');   % NULL: geometry-only run
                fieldLen = 0;
                return;
            end

            if ~isnumeric(field) || ~isvector(field)
                error('Flatland:InvalidField', ...
                      '''Field'' must be a numeric vector; got %s.', ...
                      Flatland.describe(field));
            end
            f = double(field(:));
            if ~all(isfinite(f))
                error('Flatland:InvalidField', ...
                      '''Field'' contains %d non-finite value(s).', sum(~isfinite(f)));
            end

            obj.assertFieldLength(numel(f), fieldMode, 'Field');

            pField   = libpointer('doublePtr', f);
            fieldLen = numel(f);
        end

        % -----------------------------------------------------------------
        function assertFieldLength(obj, n, fieldMode, what)
        %ASSERTFIELDLENGTH  Catch a mis-sized field before it reaches C.
        %
        %   The library would reject it anyway, but the message here can name
        %   both counts and the option that resolves an ambiguity.

            nv = obj.vertexCount;
            nf = obj.faceCount;

            switch fieldMode
                case 1      % node
                    if n ~= nv
                        error('Flatland:InvalidField', ...
                              ['''%s'' has %d row(s) but ''FieldMode'' is ''node'' ' ...
                               'and the mesh has %d vertices.'], what, n, nv);
                    end
                case 2      % face
                    if n ~= nf
                        error('Flatland:InvalidField', ...
                              ['''%s'' has %d row(s) but ''FieldMode'' is ''face'' ' ...
                               'and the mesh has %d faces.'], what, n, nf);
                    end
                otherwise   % auto
                    if nv == nf && n == nv
                        error('Flatland:AmbiguousFieldMode', ...
                              ['This mesh has %d vertices and %d faces, so a field ' ...
                               'of %d row(s) could be either. Set ''FieldMode'' to ' ...
                               '''node'' or ''face''.'], nv, nf, n);
                    end
                    if n ~= nv && n ~= nf
                        error('Flatland:InvalidField', ...
                              ['''%s'' has %d row(s), which matches neither the ' ...
                               'mesh''s %d vertices nor its %d faces.'], what, n, nv, nf);
                    end
            end
        end
    end

    % =====================================================================
    methods (Static, Access = private)

        % -----------------------------------------------------------------
        function out = rawStatic(info, fname, varargin)
        %RAWSTATIC  calllib asking for exactly the outputs MATLAB declares.
        %
        %   out{1} is the C return value for a non-void function, followed by
        %   one entry per pointer argument calllib chose to hand back. The
        %   count comes from flatland_load's introspection of libfunctions,
        %   not from an assumption about MATLAB's const-pointer rules.

            if ~isfield(info.nout, fname)
                error('Flatland:UnknownFunction', ...
                      'The loaded library does not expose %s.', fname);
            end
            n = info.nout.(fname);

            if n <= 0
                calllib(info.alias, fname, varargin{:});
                out = {[]};
                return;
            end

            cells = cell(1, n);
            try
                [cells{:}] = calllib(info.alias, fname, varargin{:});
            catch err
                % Safety net: introspection said n outputs, this release
                % disagrees. Retry asking for just the return value.
                if n > 1
                    one = cell(1, 1);
                    [one{:}] = calllib(info.alias, fname, varargin{:});
                    out = one;
                    return;
                end
                rethrow(err);
            end
            out = cells;
        end

        % -----------------------------------------------------------------
        function pp = newOutHandle(ptrType)
        %NEWOUTHANDLE  An allocated slot for a T** output argument.
        %
        %   The library writes a fresh handle through this pointer, so the slot
        %   must already exist: a bare libpointer('...PtrPtr') is NULL, which
        %   the library would dereference. Giving it an initial value is what
        %   makes MATLAB allocate the storage.

            forms = { @() libpointer([ptrType 'Ptr'], libpointer(ptrType)), ...
                      @() libpointer('voidPtrPtr',    libpointer('voidPtr')), ...
                      @() libpointer([ptrType 'Ptr']), ...
                      @() libpointer('voidPtrPtr') };
            lastErr = [];
            for k = 1:numel(forms)
                try
                    pp = forms{k}();
                    return;
                catch err
                    lastErr = err;
                end
            end
            msg = 'no diagnostic available';
            if ~isempty(lastErr), msg = lastErr.message; end
            error('Flatland:OutHandleAlloc', ...
                  'Could not allocate an output handle of type %sPtr: %s', ptrType, msg);
        end

        % -----------------------------------------------------------------
        function [status, h] = takeHandle(out, pp, fname)
        %TAKEHANDLE  Recover the handle a T** output argument was filled with.
        %
        %   Preferred source is the slot we allocated, which the library wrote
        %   through directly. Some MATLAB releases instead hand a modified
        %   pointer argument back as an extra calllib output, so that is
        %   checked as a fallback.

            status = out{1};
            h = [];

            try
                v = pp.Value;
                if ~Flatland.isNullPtr(v)
                    h = v;
                end
            catch
            end

            if isempty(h)
                for k = numel(out):-1:2
                    if isa(out{k}, 'lib.pointer') && ~Flatland.isNullPtr(out{k})
                        h = out{k};
                        break;
                    end
                end
            end

            if Flatland.statusCode(status) == 0 && isempty(h)
                error('Flatland:HandleNotReturned', ...
                      ['%s reported success but no handle came back. This is the ' ...
                       'one place the binding depends on how your MATLAB release ' ...
                       'marshals a T** output argument. Please report the output ' ...
                       'of  libfunctions(''flatland'', ''-full'')'], fname);
            end
        end

        % -----------------------------------------------------------------
        function pOut = allocResults(K)
        %ALLOCRESULTS  Room for K fl_result records, filled by one C call.

            proto = repmat(Flatland.blankResult(), 1, K);
            try
                pOut = libpointer('fl_resultPtr', proto);
            catch
                % Fallback for releases that will not build a libpointer from a
                % struct array: allocate the raw bytes and reinterpret them.
                % fl_result is 96 bytes, measured with sizeof on the build host.
                pOut = libpointer('uint8Ptr', ...
                                  zeros(1, K * Flatland.RESULT_BYTES, 'uint8'));
                setdatatype(pOut, 'fl_resultPtr', 1, K);
            end
        end

        % -----------------------------------------------------------------
        function r = blankResult()
        %BLANKRESULT  A zeroed MATLAB mirror of the C fl_result struct.
        %
        %   Field names, order and classes match fl_result exactly. Measured
        %   layout on x86-64: 96 bytes total; doubles at offsets 0/8/16/24/32,
        %   int64 at 40, int32s at 48/52/56/60, uint32[8] at 64.
        %
        %   The reserved array is wrapped in a cell so struct() reads it as one
        %   field value rather than as instructions to build an 8-element
        %   struct array.

            r = struct('area', 0, 'average', 0, 'integral', 0, ...
                       'min', 0, 'max', 0, ...
                       'covered_pixels', int64(0), ...
                       'width', int32(0), 'height', int32(0), ...
                       'has_field', int32(0), 'has_stats', int32(0), ...
                       'reserved', {zeros(1, 8, 'uint32')});
        end

        % -----------------------------------------------------------------
        function R = decodeResultArray(s, K)
            if numel(s) < K
                error('Flatland:BatchTruncated', ...
                      'The batch returned %d result(s) for %d view(s).', numel(s), K);
            end
            R = repmat(Flatland.decodeResult(s(1)), 1, K);
            for k = 2:K
                R(k) = Flatland.decodeResult(s(k));
            end
        end

        % -----------------------------------------------------------------
        function r = decodeResult(s)
        %DECODERESULT  C fl_result -> MATLAB result struct, with NaN semantics.
        %
        %   Accepts either a libstruct handle or a plain MATLAB struct.
        %
        %   has_stats is 1 only when the view covered at least one pixel AND a
        %   field was supplied. When it is 0 the C struct still contains
        %   average/integral/min/max, but they are not measurements: the
        %   library leaves them at zero. Zero is a perfectly ordinary field
        %   value, so passing it on would be indistinguishable from a real
        %   result. NaN says "not measured", and propagates through arithmetic
        %   instead of quietly dragging a mean toward zero.

            hasStats = logical(double(s.has_stats));

            r.area = double(s.area);
            if hasStats
                r.average  = double(s.average);
                r.integral = double(s.integral);
                r.min      = double(s.min);
                r.max      = double(s.max);
            else
                r.average  = NaN;
                r.integral = NaN;
                r.min      = NaN;
                r.max      = NaN;
            end
            r.coveredPixels = double(s.covered_pixels);
            r.width         = double(s.width);
            r.height        = double(s.height);
            r.hasField      = logical(double(s.has_field));
            r.hasStats      = hasStats;
        end

        % -----------------------------------------------------------------
        function V = validateVertices(V)
            if ~isnumeric(V) || isempty(V)
                error('Flatland:InvalidVertices', ...
                      'vertices must be a non-empty numeric N-by-3 matrix; got %s.', ...
                      Flatland.describe(V));
            end
            if ~ismatrix(V) || size(V, 2) ~= 3
                error('Flatland:InvalidVertices', ...
                      ['vertices must be N-by-3 (one xyz triple per ROW); got %s. ' ...
                       'If yours is 3-by-N, transpose it.'], Flatland.describe(V));
            end
            V = double(V);
            if ~all(isfinite(V(:)))
                error('Flatland:InvalidVertices', ...
                      'vertices contain %d non-finite value(s).', sum(~isfinite(V(:))));
            end
        end

        % -----------------------------------------------------------------
        function F = validateFaces(F, nv)
            if ~isnumeric(F) || isempty(F)
                error('Flatland:InvalidFaces', ...
                      'faces must be a non-empty numeric M-by-3 matrix; got %s.', ...
                      Flatland.describe(F));
            end
            if ~ismatrix(F) || size(F, 2) ~= 3
                error('Flatland:InvalidFaces', ...
                      ['faces must be M-by-3 (one triangle per ROW); got %s. ' ...
                       'If yours is 3-by-M, transpose it.'], Flatland.describe(F));
            end
            F = double(F);
            if ~all(isfinite(F(:))) || any(F(:) ~= fix(F(:)))
                error('Flatland:InvalidFaces', 'faces must contain integer indices.');
            end
            if any(F(:) < 1)
                error('Flatland:InvalidFaces', ...
                      ['faces indices are 1-BASED in this binding, but %d index(es) ' ...
                       'are below 1. Do not subtract one yourself; Flatland converts ' ...
                       'to the C API''s 0-based indices internally.'], sum(F(:) < 1));
            end
            if any(F(:) > nv)
                error('Flatland:InvalidFaces', ...
                      'faces reference vertex %d but there are only %d vertices.', ...
                      max(F(:)), nv);
            end
        end

        % -----------------------------------------------------------------
        function v = validateView(v)
            if ~isnumeric(v) || numel(v) ~= 3
                error('Flatland:InvalidView', ...
                      'The view direction must be a 3-element numeric vector; got %s.', ...
                      Flatland.describe(v));
            end
            v = double(v(:).');
            if ~all(isfinite(v))
                error('Flatland:InvalidView', ...
                      'The view direction must be finite; got [%g %g %g].', v(1), v(2), v(3));
            end
            if all(v == 0)
                error('Flatland:InvalidView', ...
                      'The view direction must be non-zero.');
            end
        end

        % -----------------------------------------------------------------
        function V = validateViews(V)
            if ~isnumeric(V) || isempty(V)
                error('Flatland:InvalidView', ...
                      'views must be a non-empty numeric K-by-3 matrix; got %s.', ...
                      Flatland.describe(V));
            end
            if isvector(V) && numel(V) == 3
                V = V(:).';                 % a single view is fine
            end
            if ~ismatrix(V) || size(V, 2) ~= 3
                error('Flatland:InvalidView', ...
                      ['views must be K-by-3 (one direction per ROW); got %s. ' ...
                       'If yours is 3-by-K, transpose it.'], Flatland.describe(V));
            end
            V = double(V);
            if ~all(isfinite(V(:)))
                error('Flatland:InvalidView', ...
                      'views contain %d non-finite value(s).', sum(~isfinite(V(:))));
            end
            zero = all(V == 0, 2);
            if any(zero)
                error('Flatland:InvalidView', ...
                      'view %d is the zero vector; every direction must be non-zero.', ...
                      find(zero, 1));
            end
        end

        % -----------------------------------------------------------------
        function code = statusCode(status)
        %STATUSCODE  Numeric fl_status, whether MATLAB gave a number or a name.

            if isnumeric(status) && ~isempty(status)
                code = double(status(1));
            elseif ischar(status) || isstring(status)
                names = {'FL_OK', 'FL_ERR_INVALID_ARGUMENT', 'FL_ERR_OUT_OF_MEMORY', ...
                         'FL_ERR_IO', 'FL_ERR_PARSE', 'FL_ERR_DIMENSION', ...
                         'FL_ERR_NUMERIC', 'FL_ERR_UNSUPPORTED', 'FL_ERR_INTERNAL'};
                idx = find(strcmpi(names, char(status)), 1);
                if isempty(idx)
                    code = -1;
                else
                    code = idx - 1;
                end
            elseif islogical(status)
                code = double(status);
            else
                code = -1;
            end
        end

        % -----------------------------------------------------------------
        function id = errorId(code)
        %ERRORID  A stable MATLAB error identifier per fl_status code.
            switch code
                case 1,    id = 'Flatland:InvalidArgument';
                case 2,    id = 'Flatland:OutOfMemory';
                case 3,    id = 'Flatland:IOError';
                case 4,    id = 'Flatland:ParseError';
                case 5,    id = 'Flatland:DimensionMismatch';
                case 6,    id = 'Flatland:NumericError';
                case 7,    id = 'Flatland:Unsupported';
                case 8,    id = 'Flatland:InternalError';
                otherwise, id = 'Flatland:Error';
            end
        end

        % -----------------------------------------------------------------
        function v = mapEnum(value, names, codes, what)
        %MAPENUM  Accept a name or an already-numeric enum value.

            if isnumeric(value) && isscalar(value)
                if ~ismember(double(value), unique(codes))
                    error('Flatland:InvalidOption', ...
                          '''%s'' value %g is not valid. Accepted: %s.', ...
                          what, double(value), strjoin(names, ', '));
                end
                v = double(value);
                return;
            end
            if ~(ischar(value) || isstring(value))
                error('Flatland:InvalidOption', ...
                      '''%s'' must be one of: %s.', what, strjoin(names, ', '));
            end
            value = char(value);
            idx = find(strcmpi(names, value), 1);
            if isempty(idx)
                error('Flatland:InvalidOption', ...
                      'Unknown ''%s'' value "%s". Accepted: %s.', ...
                      what, value, strjoin(names, ', '));
            end
            v = codes(idx);
        end

        % -----------------------------------------------------------------
        function setEnum(s, field, code, names)
        %SETENUM  Write an enum field, by numeric value or by enumerator name.
            try
                s.(field) = int32(code);
            catch
                s.(field) = names{code + 1};
            end
        end

        % -----------------------------------------------------------------
        function tf = isNullPtr(p)
        %ISNULLPTR  True for [], a NULL lib.pointer, or a zero address.
        %
        %   The numeric case matters: reading .Value off a T** can yield a raw
        %   address rather than a lib.pointer, and address zero is NULL. Without
        %   that arm a failed allocation would be mistaken for a live handle.

            tf = true;
            if isempty(p)
                return;
            end
            if isa(p, 'lib.pointer')
                try
                    tf = p.isNull;
                catch
                    tf = false;
                end
                return;
            end
            if isnumeric(p)
                tf = all(p(:) == 0);
                return;
            end
            tf = false;
        end

        % -----------------------------------------------------------------
        function s = toChar(v)
        %TOCHAR  calllib returns a C string as char; normalise anything else.
            if ischar(v)
                s = v;
            elseif isstring(v)
                s = char(v);
            elseif isa(v, 'lib.pointer')
                try
                    s = char(v.Value);
                catch
                    s = '';
                end
            elseif isnumeric(v)
                s = char(v(:).');
            else
                s = '';
            end
            s = s(:).';
            z = find(s == char(0), 1);
            if ~isempty(z)
                s = s(1:z-1);
            end
        end

        % -----------------------------------------------------------------
        function d = describe(v)
        %DESCRIBE  Short human description of a rejected value.
            try
                d = sprintf('%s %s', mat2str(size(v)), class(v));
            catch
                d = class(v);
            end
        end
    end
end
