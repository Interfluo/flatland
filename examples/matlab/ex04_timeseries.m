% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% Example 4 - A whole time series in one call
%
% The question this answers: *I have hundreds of timesteps, each with its own
% view direction and its own field. How do I not wait all afternoon?*
%
% projectBatch takes every view and the entire field matrix across the boundary
% in ONE calllib, then runs the views across all cores. The mesh is validated and
% recentred a single time no matter how many timesteps there are. Looping over
% project() from MATLAB instead pays the marshalling cost per step and throws
% away the parallelism.
%
%     >> cd /path/to/flatland/examples/matlab
%     >> ex04_timeseries
%
% Covered:
%   - one field matrix, one column per timestep
%   - 'FieldColumns' is 1-BASED
%   - what the batch path actually buys, measured rather than asserted
%   - per-view resolutions
%   - guarding on hasStats before aggregating

clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matlab'));
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

flatland_load();
fprintf('FlatLand %s\n', Flatland.versionString());

fl = Flatland.load(fullfile(root, 'examples', 'with_fields', 'bunny.obj'));
nv = fl.vertexCount;
fprintf('loaded the bunny: %d vertices, %d triangles\n\n', nv, fl.faceCount);

STEPS = 240;

%% Building the time series
fprintf('Building the time series\n------------------------\n');

% A camera orbiting in azimuth at a fixed elevation.
az    = 360 * (0:STEPS-1)' / STEPS;
views = Flatland.angleToDir(az, 20 * ones(STEPS, 1));    % STEPS-by-3

% The field matrix is ONE array: one ROW per mesh entity, one COLUMN per
% timestep. A travelling wave along x here - substitute your own solver output.
% Shape it exactly like this and the binding handles the row-major transpose.
V    = fl.vertices();
x    = V(:,1);
nx   = (x - min(x)) / max(max(x) - min(x), eps);
tIdx = 0:STEPS-1;
fieldMatrix = sin(4*pi*nx - 2*pi*tIdx/STEPS);            % nv-by-STEPS

fprintf('  %d views, field matrix %d x %d (%d values)\n', ...
        STEPS, size(fieldMatrix,1), size(fieldMatrix,2), numel(fieldMatrix));

%% Running it
fprintf('\nRunning it\n----------\n');

% 'FieldColumns' is 1-BASED: timestep t reads column t.
tic;
R = fl.projectBatch(views, ...
                    'FieldMatrix',  fieldMatrix, ...
                    'FieldColumns', 1:STEPS, ...
                    'Resolution',   1e-3, ...
                    'Threads',      0);          % 0 = one worker per core
batchTime = toc;

fprintf('  %d timesteps in %.3f s  (%.0f steps/s, %.2f ms/step)\n', ...
        numel(R), batchTime, numel(R)/batchTime, 1e3*batchTime/numel(R));

%% What the batch call bought
% The same views one at a time, with the field columns pulled out BEFORE the
% timed region so this measures FlatLand rather than array slicing.
fprintf('\nWhat the batch call bought\n--------------------------\n');

nCmp = 24;
cols = cell(nCmp, 1);
for t = 1:nCmp
    cols{t} = fieldMatrix(:, t);
end

tic;
for t = 1:nCmp
    fl.project(views(t,:), 'Field', cols{t}, 'Resolution', 1e-3);
end
loopTime = toc;

% Single-threaded batch, to separate "crossing the boundary once" from
% "running on every core".
tic;
fl.projectBatch(views(1:nCmp,:), 'FieldMatrix', fieldMatrix, ...
                'FieldColumns', 1:nCmp, 'Resolution', 1e-3, 'Threads', 1);
serialBatch = toc;

fprintf('  one at a time        %6.2f ms/step\n', 1e3*loopTime/nCmp);
fprintf('  batched, 1 thread    %6.2f ms/step\n', 1e3*serialBatch/nCmp);
fprintf('  batched, all cores   %6.2f ms/step\n', 1e3*batchTime/STEPS);
fprintf('  speedup              %6.2fx\n', (loopTime/nCmp) / (batchTime/STEPS));

%% The results
% Always check hasStats before aggregating. A view that covered no pixels has no
% field statistics, and this binding reports them as NaN rather than 0 precisely
% so they cannot be averaged in by accident.
fprintf('\nThe results\n-----------\n');

good = logical([R.hasStats]);
if ~all(good)
    fprintf('  %d of %d views covered no pixels and were excluded\n', ...
            sum(~good), numel(R));
end

areas     = [R.area];
integrals = [R(good).integral];

fprintf('  area      min %.6f  max %.6f  mean %.6f\n', ...
        min(areas), max(areas), mean(areas));
fprintf('  integral  min %+.6f  max %+.6f  mean %+.6f\n', ...
        min(integrals), max(integrals), mean(integrals));

[peakArea, iPeak] = max(areas);
fprintf('  largest silhouette at timestep %d (azimuth %.0f deg): %.6f\n', ...
        iPeak, az(iPeak), peakArea);

fprintf('\n  first six timesteps:\n');
fprintf('  %-5s %-11s %-13s %-13s %-8s\n', 'step', 'area', 'mean f', 'integral', 'pixels');
for t = 1:6
    fprintf('  %-5d %-11.6f %-+13.6f %-+13.6f %-8d\n', ...
            t, R(t).area, R(t).average, R(t).integral, R(t).coveredPixels);
end

%% Per-view resolution
% Each view can carry its own pixel size - useful for a cheap coarse pass before
% a fine one, or when some orientations need more detail than others.
fprintf('\nPer-view resolution\n-------------------\n');
resPerView = [repmat(4e-3, 1, 4), repmat(1e-3, 1, 4)];
S = fl.projectBatch(views(1:8,:), 'Resolutions', resPerView);
for t = 1:8
    fprintf('  step %d  res %.4f  area %.6f  raster %dx%d\n', ...
            t, resPerView(t), S(t).area, S(t).width, S(t).height);
end

%% Plot the series, if this session can show a figure
try
    figure('Name', 'FlatLand time series');
    subplot(2,1,1);
    plot(az, areas, 'LineWidth', 1.2); grid on;
    xlabel('azimuth (deg)'); ylabel('projected area'); xlim([0 360]);
    title('Silhouette area as the camera orbits');

    subplot(2,1,2);
    allIntegrals = nan(1, numel(R));
    allIntegrals(good) = integrals;
    plot(az, allIntegrals, 'LineWidth', 1.2); grid on;
    xlabel('azimuth (deg)'); ylabel('\int f dA');
    title('Area integral of the evolving field');
catch
    fprintf('\n  (no display available, skipping the figure)\n');
end

delete(fl);
