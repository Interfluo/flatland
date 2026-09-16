% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% FlatLand from MATLAB - end-to-end example
%
% Run this from the matlab/ folder (or with that folder on the MATLAB path):
%
%     >> cd /path/to/flatland/matlab
%     >> example
%
% It needs libflatland built first:
%
%     $ cd /path/to/flatland && make lib
%
% Everything here is 1-BASED: faces index vertices from 1, and FieldColumns
% index the field matrix from 1. The binding converts to the C API's 0-based
% indices internally. See "Indexing" in README.md.

clear; close all; %#ok<CLALL>

%% 0. Check the library is there
info = flatland_load();
fprintf('FlatLand %d.%d.%d\n', info.version);
fprintf('library: %s\n', info.libfile);
fprintf('header : %s\n\n', info.header);

%% 1. A mesh from MATLAB arrays
% A closed unit box spanning [0,1]^3, with counter-clockwise outward normals.
% 8 vertices, 12 triangles.
V = [0 0 0; 1 0 0; 1 1 0; 0 1 0; ...
     0 0 1; 1 0 1; 1 1 1; 0 1 1];

% 1-BASED face indices. Do not subtract one; Flatland does that for you.
F = [1 4 3; 1 3 2;   5 6 7; 5 7 8; ...
     1 2 6; 1 6 5;   2 3 7; 2 7 6; ...
     3 4 8; 3 8 7;   4 1 5; 4 5 8];

box = Flatland(V, F);
fprintf('mesh: %d vertices, %d faces\n', box.vertexCount, box.faceCount);

%% 2. Geometry only: how much area does the box present along +X?
r = box.project([1 0 0], 'Resolution', 0.002, 'Precision', 'double');

fprintf('\n-- projected area along +X --\n');
fprintf('  area          %.6f   (a unit box should give 1)\n', r.area);
fprintf('  raster        %d x %d\n', r.width, r.height);
fprintf('  coveredPixels %d\n', r.coveredPixels);
fprintf('  hasStats      %d   <- no field was supplied\n', r.hasStats);
fprintf('  average       %g   <- NaN, because nothing was measured\n', r.average);

%% 3. With a per-vertex field
% The view direction is the direction the camera LOOKS ALONG, so looking along
% +X the visible surface is the one at x = 0. A field equal to the vertex x
% coordinate therefore averages to 0 over that surface.
fieldX = V(:, 1);                       % one value per vertex (node field)

r = box.project([1 0 0], 'Field', fieldX, ...
                'Resolution', 0.002, 'Precision', 'double');

fprintf('\n-- field = vertex x, viewed along +X --\n');
fprintf('  hasStats      %d\n', r.hasStats);
fprintf('  average       %.6f   (the near surface is x = 0)\n', r.average);
fprintf('  min / max     %.6f / %.6f\n', r.min, r.max);
fprintf('  integral      %.6f   (integral of f dA)\n', r.integral);

%% 4. A per-face field
% 12 values, one per triangle. The mesh has 8 vertices and 12 faces, so the
% length alone picks the mode; 'FieldMode' is only needed when they tie.
pressure = 4.0 * ones(box.faceCount, 1);
r = box.project([0 0 1], 'Field', pressure, ...
                'Resolution', 0.002, 'Precision', 'double');

fprintf('\n-- constant face field of 4.0, viewed along +Z --\n');
fprintf('  average       %.6f\n', r.average);
fprintf('  integral      %.6f   (== 4.0 * area = %.6f)\n', r.integral, 4.0 * r.area);

%% 5. A view that covers nothing: NaN, not zero
% A triangle thinner than one pixel covers no pixels at all, so its field
% statistics were never measured. They come back NaN rather than 0, because a
% zero here would be indistinguishable from a real measurement of zero.
sliver = Flatland([0 0 0; 1 0 0; 1 1e-6 0], [1 2 3]);
rs = sliver.project([0 0 -1], 'Field', [-5; -3; -1], ...
                    'Resolution', 0.1, 'Cull', false);

fprintf('\n-- a sliver thinner than one pixel --\n');
fprintf('  coveredPixels %d\n', rs.coveredPixels);
fprintf('  hasField      %d   <- a field WAS supplied\n', rs.hasField);
fprintf('  hasStats      %d   <- but nothing was measured\n', rs.hasStats);
fprintf('  average       %g\n', rs.average);
fprintf('  area          %g   <- still a real measurement, so not NaN\n', rs.area);
delete(sliver);

%% 6. A batch: a time series of views and fields in ONE call
% 36 views sweeping the azimuth, and a field matrix with one ROW per vertex
% and one COLUMN per timestep. That is the MATLAB-natural orientation; the
% binding transposes it into the row-major layout the C API reads.
nViews = 36;
az = linspace(0, 350, nViews).';
views = Flatland.angleToDir(az, zeros(nViews, 1));      % nViews-by-3

% A field that grows linearly with the timestep, so the per-view averages are
% easy to eyeball: column t is the vertex x coordinate scaled by t.
nSteps = nViews;
M = V(:, 1) * (1:nSteps);                % 8-by-36, entity-by-timestep

tic;
R = box.projectBatch(views, ...
                     'FieldMatrix',  M, ...
                     'FieldColumns', 1:nViews, ...   % 1-BASED, one per view
                     'Resolution',   0.004, ...
                     'Precision',    'double', ...
                     'Threads',      0);             % 0 = one worker per core
elapsed = toc;

areas    = [R.area];
averages = [R.average];

fprintf('\n-- batch of %d views in one calllib --\n', nViews);
fprintf('  elapsed       %.3f s\n', elapsed);
fprintf('  area   min/max  %.4f / %.4f   (a box is 1 face-on, sqrt(2) at 45 deg)\n', ...
        min(areas), max(areas));
fprintf('  every view has statistics: %d\n', all([R.hasStats]));

% Anything with hasStats false would carry NaN, so guard aggregates.
good = [R.hasStats];
if any(good)
    fprintf('  mean of the measured averages: %.4f\n', mean(averages(good)));
end

%% 7. Render a raster
img = box.render([1 1 0.4], 'Field', fieldX, ...
                 'Resolution', 0.004, 'Precision', 'double');

fprintf('\n-- rendered raster --\n');
fprintf('  size          %d x %d\n', img.width, img.height);
fprintf('  covered       %d of %d pixels\n', nnz(img.mask), numel(img.mask));
fprintf('  values range  %.4f .. %.4f (NaN off the silhouette)\n', ...
        min(img.values(:)), max(img.values(:)));

if usejava('jvm') && ~isempty(img.values)
    figure('Name', 'FlatLand render', 'Color', 'w');

    subplot(1, 2, 1);
    % Row 1 of mask/values is the BOTTOM row in mesh space, matching the C
    % API, so flip for the usual top-down image orientation.
    imagesc(flipud(double(img.mask)));
    axis image off; colormap(gca, gray); title('coverage mask');

    subplot(1, 2, 2);
    h = imagesc(flipud(img.values));
    set(h, 'AlphaData', flipud(img.mask));   % leave uncovered pixels blank
    axis image off; colorbar; title('field = vertex x');
end

%% 8. A mesh from a file
repoRoot = fileparts(fileparts(mfilename('fullpath')));   % matlab/.. == repo root
objPath  = fullfile(repoRoot, 'examples', 'cube_area', 'cube.obj');
if exist(objPath, 'file') == 2
    cube = Flatland.load(objPath);
    rc = cube.project([1 0 0], 'Resolution', 0.005, 'Precision', 'double');
    fprintf('\n-- %s --\n', objPath);
    fprintf('  %d vertices, %d faces\n', cube.vertexCount, cube.faceCount);
    fprintf('  area along +X: %.6f\n', rc.area);

    % OBJ files on disk are 1-based, the C API is 0-based, and this binding
    % gives them back 1-based. So the round trip is consistent.
    Fc = cube.faces();
    fprintf('  face indices span %d..%d (1-based)\n', min(Fc(:)), max(Fc(:)));
    delete(cube);
else
    fprintf('\n(skipped the file example: %s not found)\n', objPath);
end

%% 9. Clean up
% delete() frees the C mesh handle. Flatland is a handle class, so this is
% deterministic rather than waiting on garbage collection. Letting the
% variable go out of scope would free it too.
delete(box);

fprintf('\nDone.\n');
