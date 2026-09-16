% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% Example 1 - Projected area, and finding the orientation that minimises it
%
% The question this answers: *how much of my part does a given direction see,
% and which orientation exposes the least (or most) of it?* That shows up as
% frontal area for drag, as an aperture for radiative load, as a footprint for
% packing.
%
%     >> cd /path/to/flatland/examples/matlab
%     >> ex01_projected_area
%
% Needs the library built first:  cd /path/to/flatland && make lib
%
% Covered:
%   - loading an OBJ and projecting along a few directions
%   - an azimuth/elevation sweep in a single parallel batch call
%   - checking the answer against a closed form, because for a cube there is one
%   - choosing a resolution by watching the answer converge

clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matlab'));
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

info = flatland_load();
fprintf('FlatLand %d.%d.%d\n', info.version);

fl = Flatland.load(fullfile(root, 'examples', 'cube_area', 'cube.obj'));
fprintf('loaded a mesh with %d vertices and %d triangles\n\n', ...
        fl.vertexCount, fl.faceCount);

%% A few directions
% A view direction is the direction the camera LOOKS ALONG. Magnitude is
% irrelevant, so [2 0 0] and [1 0 0] are the same view.
fprintf('A few directions\n----------------\n');
views = [1 0 0; 0 1 0; 1 1 0; 1 1 1];
for k = 1:size(views, 1)
    r = fl.project(views(k,:), 'Resolution', 1e-3);
    fprintf('  view [%4.1f %4.1f %4.1f]   area %.6f   (%d pixels, %dx%d raster)\n', ...
            views(k,:), r.area, r.coveredPixels, r.width, r.height);
end

%% Is that right?
% This mesh is a unit cube, and a convex polyhedron's projected area has a
% closed form:  A(n) = 1/2 sum_i A_i |n . n_i|.  For a unit cube that reduces
% to |nx| + |ny| + |nz|. Worth checking against, because a silently wrong area
% looks exactly like a right one.
fprintf('\nIs that right?\n--------------\n');
cubeExact = @(v) sum(abs(v)) / norm(v);

checks = [1 0 0; 1 1 0; 1 1 1; 2 -1 3; 0.3 0.9 -0.31];
worst = 0;
for k = 1:size(checks, 1)
    rk   = fl.project(checks(k,:), 'Resolution', 5e-4, 'Precision', 'double');
    got  = rk.area;
    want = cubeExact(checks(k,:));
    rel  = abs(got - want) / want;
    worst = max(worst, rel);
    fprintf('  [%5.2f %5.2f %5.2f]  measured %.7f   exact %.7f   rel err %.1e\n', ...
            checks(k,:), got, want, rel);
end
fprintf('  worst relative error: %.1e\n', worst);

%% Sweeping orientation to find the extremes
% 60 directions around the object, evaluated in ONE call. projectBatch parses
% and centres the mesh once and runs the views across all cores; looping over
% project() would redo that work every time.
fprintf('\nSweeping orientation to find the extremes\n');
fprintf('-----------------------------------------\n');

az = 0:30:330;
el = [-60 -30 0 30 60];
[AZ, EL] = meshgrid(az, el);
AZ = AZ(:); EL = EL(:);

sweepViews = Flatland.angleToDir(AZ, EL);     % K-by-3 unit directions
R = fl.projectBatch(sweepViews, 'Resolution', 2e-3, 'Threads', 0);

areas = [R.area];
[minArea, iMin] = min(areas);
[maxArea, iMax] = max(areas);

fprintf('  %d orientations evaluated in one batch\n', numel(R));
fprintf('  minimum  %.5f at azimuth %3d, elevation %3d\n', minArea, AZ(iMin), EL(iMin));
fprintf('  maximum  %.5f at azimuth %3d, elevation %3d\n', maxArea, AZ(iMax), EL(iMax));
fprintf('  mean     %.5f over all orientations\n', mean(areas));

% For a convex body the mean projected area over ALL directions is exactly a
% quarter of the surface area (Cauchy). A unit cube has S = 6, so the true mean
% is 1.5 - this grid is coarse and biased toward the equator, so it only lands
% near it, but a full sweep converges.
fprintf('  (Cauchy''s formula gives S/4 = 1.5 for the true all-direction mean)\n');

%% Choosing a resolution
% Coverage is decided by testing whether a pixel CENTRE falls inside a triangle,
% so the answer converges as the pixel shrinks. Sweep it and stop when the digits
% you care about stop moving. Convergence is first order and not monotonic, so
% look at the trend rather than one step.
fprintf('\nChoosing a resolution\n---------------------\n');
probe = [0.3 0.9 -0.31];
exact = cubeExact(probe);
fprintf('  %-10s %-12s %-10s\n', 'pixel', 'area', 'rel err');
for res = [2e-2 1e-2 5e-3 2e-3 1e-3 5e-4]
    rr = fl.project(probe, 'Resolution', res, 'Precision', 'double');
    a  = rr.area;
    fprintf('  %-10.5f %-12.7f %.2e\n', res, a, abs(a - exact) / exact);
end

%% Plot the sweep, if this session can show a figure
try
    figure('Name', 'Projected area vs orientation');
    A = reshape(areas, numel(el), numel(az));
    imagesc(az, el, A); axis xy; colorbar;
    xlabel('azimuth (deg)'); ylabel('elevation (deg)');
    title('Projected area of the unit cube');
catch
    fprintf('\n  (no display available, skipping the figure)\n');
end

delete(fl);
