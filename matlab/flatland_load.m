function varargout = flatland_load(cmd)
%FLATLAND_LOAD  Load, unload and introspect the FlatLand shared library.
%
%   INFO = FLATLAND_LOAD() loads libflatland (if it is not already loaded) and
%   returns a struct describing it. Calling it repeatedly is cheap: the
%   introspection is cached and the library is only loaded once per session.
%
%   FLATLAND_LOAD('unload')   unloads the library and clears the cache.
%   FLATLAND_LOAD('reload')   unloads and loads again (useful after rebuilding).
%   TF = FLATLAND_LOAD('isloaded')
%   INFO = FLATLAND_LOAD('info')   loads if needed, same as FLATLAND_LOAD().
%
%   INFO fields:
%     alias       the loadlibrary alias, 'flatland'
%     libfile     absolute path of the shared library that was loaded
%     header      absolute path of the header loadlibrary parsed
%     version     [major minor patch] reported by the library itself
%     nout        struct mapping each C function name to the number of outputs
%                 calllib produces for it (see "Output counts" below)
%     structs     cell array of the struct type names loadlibrary registered
%     meshPtrType name of the type loadlibrary chose for 'fl_mesh *'
%     warnings    cell array of warnings loadlibrary emitted while parsing
%
%   WHICH LIBRARY IS LOADED
%   Searched in order:
%     1. the path in the FLATLAND_LIB environment variable, if set
%     2. the repository root next to this file (../libflatland.so|.dylib|.dll)
%     3. libflatland on the system loader path
%
%   WHICH HEADER IS PARSED
%   flatland_matlab.h in this folder, never include/flatland.h directly. That
%   shim exists because every declaration in the real header is prefixed with
%   FL_API, which expands to __attribute__((visibility("default"))) on
%   Linux/macOS and __declspec(dllimport) on Windows -- neither of which is in
%   the C subset loadlibrary parses. See the comments in flatland_matlab.h.
%
%   OUTPUT COUNTS
%   calllib returns the C return value first, then an extra output for each
%   pointer argument it decides the callee may have written through. The exact
%   set depends on the MATLAB release, so rather than assume it, this function
%   asks MATLAB via libfunctions(alias,'-full') and records the real count.
%   Flatland.m uses those counts so it never mis-indexes a call's outputs.
%
%   See also FLATLAND, FLATLAND_TEST.

    if nargin < 1
        cmd = 'info';
    end
    cmd = lower(strtrim(char(cmd)));

    persistent INFO
    alias = 'flatland';

    switch cmd
        case 'isloaded'
            varargout{1} = libisloaded(alias);
            return;

        case 'unload'
            if libisloaded(alias)
                unloadlibrary(alias);
            end
            INFO = [];
            if nargout > 0, varargout{1} = false; end
            return;

        case 'reload'
            if libisloaded(alias)
                unloadlibrary(alias);
            end
            INFO = [];

        case {'info', 'load'}
            % fall through

        otherwise
            error('Flatland:BadCommand', ...
                  'Unknown command "%s". Use load, unload, reload, isloaded or info.', cmd);
    end

    % If another part of the session unloaded the library behind our back, the
    % cached introspection is stale.
    if ~isempty(INFO) && ~libisloaded(alias)
        INFO = [];
    end

    if isempty(INFO)
        INFO = doLoad(alias);
    end

    varargout{1} = INFO;
end

% -------------------------------------------------------------------------

function info = doLoad(alias)

    here   = fileparts(mfilename('fullpath'));
    header = fullfile(here, 'flatland_matlab.h');
    stubs  = fullfile(here, 'stubinc');

    if exist(header, 'file') ~= 2
        error('Flatland:HeaderMissing', ...
              'Cannot find the loadlibrary shim header at:\n  %s', header);
    end

    libfile = findLibrary(here);

    if ~libisloaded(alias)
        % loadlibrary needs a configured C compiler to preprocess the header.
        % Check first so the failure is a clear message instead of a parser
        % error forty lines deep.
        assertCompiler();

        % stubinc shadows <stddef.h> and <stdint.h> for the parse only. GCC's
        % real <stddef.h> defines max_align_t with long double, __attribute__
        % and an unnamed struct, none of which loadlibrary's parser accepts.
        [notfound, warns] = loadlibrary(libfile, header, ...
                                        'alias', alias, ...
                                        'includepath', stubs);

        notfound = toCellOfChar(notfound);
        if ~isempty(notfound)
            % Declared in the header but missing from the binary: almost always
            % a stale library next to a newer header.
            unloadlibrary(alias);
            error('Flatland:SymbolsMissing', ...
                  ['The library at\n  %s\ndoes not export %d function(s) the ' ...
                   'header declares:\n  %s\nRebuild it with "make lib" from the ' ...
                   'repository root.'], ...
                  libfile, numel(notfound), strjoin(notfound, ', '));
        end
    else
        warns = {};
    end

    info = struct();
    info.alias    = alias;
    info.libfile  = libfile;
    info.header   = header;
    info.warnings = toCellOfChar(warns);
    info.structs  = sort(toCellOfChar(libstructs(alias)));
    info.nout     = signatureOutputCounts(alias);

    % loadlibrary is free to name the type for an incomplete struct such as
    % 'fl_mesh' either after the tag (fl_meshPtr) or generically (voidPtr).
    % Both work; we just have to construct out-parameters with the name MATLAB
    % actually chose, so read it back off the real signature.
    info.meshPtrType  = argTypeName(alias, 'fl_mesh_destroy',  1, 'voidPtr');
    info.imagePtrType = argTypeName(alias, 'fl_image_destroy', 1, 'voidPtr');
    info.fmatPtrType  = argTypeName(alias, 'fl_field_matrix_destroy', 1, 'voidPtr');

    % Ask the library for its own version. This is the first real call across
    % the boundary, so it doubles as a smoke test that marshalling works.
    %
    % fl_version returns void and writes through three int* arguments, so the
    % outputs calllib produces are exactly those three -- but ask MATLAB rather
    % than assume, and cope with it disagreeing.
    n = 0;
    if isfield(info.nout, 'fl_version')
        n = info.nout.fl_version;
    end

    info.version = [NaN NaN NaN];
    if n >= 3
        out = cell(1, n);
        try
            [out{:}] = calllib(alias, 'fl_version', int32(0), int32(0), int32(0));
            % The three int* values are the last three outputs.
            info.version = double([out{n-2}, out{n-1}, out{n}]);
        catch
        end
    end

    if any(isnan(info.version))
        warning('Flatland:VersionUnavailable', ...
                ['calllib reported %d output(s) for fl_version, expected 3, so the ' ...
                 'library version could not be confirmed. The binding will still ' ...
                 'work.'], n);
    elseif info.version(1) ~= 3
        warning('Flatland:VersionMismatch', ...
                ['This binding was written against FlatLand 3.x; the loaded ' ...
                 'library reports %d.%d.%d.'], info.version(1), info.version(2), info.version(3));
    end
end

% -------------------------------------------------------------------------

function libfile = findLibrary(here)
%FINDLIBRARY  Locate the FlatLand shared library.

    if ispc
        names = {'flatland.dll', 'libflatland.dll'};
    elseif ismac
        names = {'libflatland.dylib'};
    else
        names = {'libflatland.so'};
    end

    candidates = {};

    env = getenv('FLATLAND_LIB');
    if ~isempty(env)
        candidates{end+1} = env;                              %#ok<AGROW>
    end

    root = fileparts(here);          % matlab/.. == repository root
    for k = 1:numel(names)
        candidates{end+1} = fullfile(root, names{k});         %#ok<AGROW>
        candidates{end+1} = fullfile(here, names{k});         %#ok<AGROW>
        candidates{end+1} = fullfile(root, 'build', names{k});%#ok<AGROW>
    end

    for k = 1:numel(candidates)
        if exist(candidates{k}, 'file') == 2
            libfile = candidates{k};
            return;
        end
    end

    % Last resort: let the system loader find it by name.
    libfile = 'libflatland';
    if ~isempty(env)
        error('Flatland:LibraryMissing', ...
              'FLATLAND_LIB is set to "%s" but no file is there.', env);
    end
    warning('Flatland:LibraryNotFound', ...
            ['Could not find %s next to the repository root (%s). Falling back ' ...
             'to the system loader path. Build it with "make lib", or set the ' ...
             'FLATLAND_LIB environment variable to its full path.'], ...
            names{1}, root);
end

% -------------------------------------------------------------------------

function assertCompiler()
%ASSERTCOMPILER  loadlibrary preprocesses the header with a C compiler.

    try
        cc = mex.getCompilerConfigurations('C', 'Selected');
    catch
        return;   % older MATLAB without this API; let loadlibrary speak for itself
    end
    if isempty(cc)
        error('Flatland:NoCompiler', ...
              ['loadlibrary must preprocess flatland_matlab.h with a C compiler, ' ...
               'and MATLAB has none configured.\n' ...
               'Run "mex -setup C" and pick a supported compiler:\n' ...
               '  Linux   gcc\n' ...
               '  macOS   Xcode command line tools\n' ...
               '  Windows MinGW-w64 (free, install from the Add-On Explorer) or MSVC']);
    end
end

% -------------------------------------------------------------------------

function nout = signatureOutputCounts(alias)
%SIGNATUREOUTPUTCOUNTS  Number of outputs calllib produces for each function.
%
%   libfunctions(alias,'-full') returns strings shaped like
%       '[fl_status, doublePtr, voidPtrPtr] fl_mesh_create(doublePtr, ...)'
%   when there is more than one output, or
%       'int32 fl_image_width(voidPtr)'
%   when the return value is the only one. Parsing this is how the binding
%   learns the real output arity instead of guessing at MATLAB's const rules.

    sigs = libfunctions(alias, '-full');
    nout = struct();

    for k = 1:numel(sigs)
        s = strtrim(sigs{k});
        n = 0;

        if ~isempty(s) && s(1) == '['
            close = find(s == ']', 1, 'first');
            if isempty(close), continue; end
            inner = strtrim(s(2:close-1));
            if ~isempty(inner)
                n = numel(strsplit(inner, ','));
            end
            rest = strtrim(s(close+1:end));
        else
            rest = s;
            % A leading identifier followed by whitespace is the return type.
            tok = regexp(rest, '^([A-Za-z_]\w*)\s+(?=[A-Za-z_])', 'tokens', 'once');
            if ~isempty(tok) && ~strcmp(tok{1}, 'void')
                n = 1;
            end
        end

        name = regexp(rest, '([A-Za-z_]\w*)\s*\(', 'tokens', 'once');
        if isempty(name), continue; end
        nout.(name{1}) = n;
    end

    if ~isfield(nout, 'fl_project') || ~isfield(nout, 'fl_mesh_create')
        error('Flatland:SignatureParse', ...
              ['Could not read the function signatures back from MATLAB. ' ...
               'Run libfunctions(''%s'',''-full'') and report what it prints.'], alias);
    end
end

% -------------------------------------------------------------------------

function c = toCellOfChar(v)
%TOCELLOFCHAR  Normalise loadlibrary's assorted return shapes to a cell array.
%
%   loadlibrary's notfound/warnings and libstructs return a cell array of char
%   in the common case, but an empty char or a string array in others.

    if isempty(v)
        c = {};
        return;
    end
    if ischar(v)
        c = {v};
    elseif iscell(v)
        c = v(:).';
        keep = true(1, numel(c));
        for k = 1:numel(c)
            if isstring(c{k})
                c{k} = char(c{k});
            elseif ~ischar(c{k})
                keep(k) = false;
            end
            if keep(k) && isempty(c{k})
                keep(k) = false;
            end
        end
        c = c(keep);
    elseif isstring(v)
        c = cellstr(v(:).');
    else
        c = {};
    end
end

% -------------------------------------------------------------------------

function t = argTypeName(alias, fname, argIndex, fallback)
%ARGTYPENAME  Type loadlibrary assigned to one argument of one function.

    t = fallback;
    try
        sigs = libfunctions(alias, '-full');
        for k = 1:numel(sigs)
            s = strtrim(sigs{k});
            m = regexp(s, ['(?<![A-Za-z_])' fname '\s*\(([^)]*)\)'], 'tokens', 'once');
            if isempty(m), continue; end
            args = strtrim(strsplit(m{1}, ','));
            if numel(args) >= argIndex && ~isempty(args{argIndex})
                t = args{argIndex};
            end
            return;
        end
    catch
        % Keep the fallback; voidPtr is accepted wherever an opaque handle goes.
    end
end
