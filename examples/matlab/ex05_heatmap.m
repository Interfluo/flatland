% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% Example 5 - Getting at the raster itself
%
% The question this answers: *I do not just want the number, I want to see where
% it came from.* render() hands back the coverage mask and the interpolated
% per-pixel field values as ordinary MATLAB matrices, ready for imagesc,
% imwrite, or whatever comes next in your pipeline.
%
%     >> cd /path/to/flatland/examples/matlab
%     >> ex05_heatmap
%
% Covered:
%   - the coverage mask, checked against coveredPixels
%   - per-pixel field values, and recomputing the integral from them by hand
%   - displaying with imagesc, and the ROW ORDER that catches people out
%   - writing a PNG with imwrite
%   - what an empty view gives you

clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matlab'));
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

flatland_load();
fprintf('FlatLand %s\n\n', Flatland.versionString());

fl   = Flatland.load(fullfile(root, 'examples', 'sphere_areas', 'sphere_fine.obj'));
viewDir = [0.4 0.3 1.0];
res  = 4e-3;

% f = cos(incidence): brightest facing the camera, falling to zero at the limb.
V     = fl.vertices();
u     = viewDir / norm(viewDir);
field = -(V * u(:)) ./ sqrt(sum(V.^2, 2));

%% Rendering
fprintf('Rendering\n---------\n');
img = fl.render(viewDir, 'Field', field, 'Resolution', res, 'Precision', 'double');
r   = img.result;

fprintf('  raster      %d x %d\n', img.width, img.height);
fprintf('  covered     %d pixels of %d\n', r.coveredPixels, img.width * img.height);
fprintf('  area        %.6f\n', r.area);
fprintf('  field range %.4f .. %.4f, mean %.4f\n', r.min, r.max, r.average);

%% The mask
% mask is height-by-width logical, true where a pixel is covered.
fprintf('\nThe mask\n--------\n');
fprintf('  mask size   %d x %d  (matches the raster: %d)\n', ...
        size(img.mask,1), size(img.mask,2), ...
        isequal(size(img.mask), [img.height img.width]));
fprintf('  covered in the mask: %d  (matches coveredPixels: %d)\n', ...
        nnz(img.mask), nnz(img.mask) == r.coveredPixels);

%% The values, and recomputing the integral by hand
% values holds the interpolated field at each covered pixel, with NaN off the
% silhouette. Summing it and multiplying by the pixel area reproduces the
% integral FlatLand reported - which is exactly how it is defined.
fprintf('\nThe values, and recomputing the integral by hand\n');
fprintf('------------------------------------------------\n');

pixelArea = res^2;
manual    = sum(img.values(img.mask)) * pixelArea;
fprintf('  FlatLand''s integral : %.8f\n', r.integral);
fprintf('  recomputed by hand  : %.8f\n', manual);
fprintf('  agree to 1e-9       : %d\n', abs(manual - r.integral) < 1e-9);

%% Row order - the one that catches people out
% ROW 1 IS THE BOTTOM ROW in mesh space, matching the C API. Most MATLAB image
% functions put row 1 at the top, so either use axis xy with imagesc, or flipud
% before imwrite. Get it wrong and your heatmap is mirrored vertically - which
% looks entirely plausible on a symmetric object.
fprintf('\nRow order\n---------\n');
fprintf('  img.values(1,:)   is the BOTTOM row in mesh space\n');
fprintf('  img.values(end,:) is the TOP row\n');
fprintf('  imagesc(img.values); axis xy      %% correct orientation\n');
fprintf('  imwrite(flipud(...), ''out.png'')   %% flip for top-down formats\n');

%% Display
try
    figure('Name', 'FlatLand raster');

    subplot(1,2,1);
    imagesc(img.values); axis image xy off;
    colormap(jet); colorbar;
    title('cos(incidence)');

    subplot(1,2,2);
    imagesc(double(img.mask)); axis image xy off;
    colormap(gray);
    title('coverage mask');
catch
    fprintf('\n  (no display available, skipping the figure)\n');
end

%% Writing a PNG
% Map the field through a colormap, paint the background, flip to top-down.
fprintf('\nWriting a PNG\n-------------\n');

lo   = r.min;
span = max(r.max - r.min, eps);
norm01 = (img.values - lo) / span;
norm01(~img.mask) = 0;

cmap = jet(256);
% Keep the indices double: uint8 saturates at 255, so uint8(255)+1 would fold
% the top of the colour ramp back onto itself instead of reaching index 256.
idx  = round(min(max(norm01, 0), 1) * 255) + 1;      % 1 .. 256
rgb  = ind2rgb(idx, cmap);                            % height-by-width-by-3

% Background for the uncovered pixels, matching FlatLand's own PNG output.
bg   = reshape([30 30 35] / 255, 1, 1, 3);
keep = double(img.mask);
rgb  = rgb .* keep + bg .* (1 - keep);                % implicit expansion, R2016b+

outPng = fullfile(tempdir, 'sphere_heatmap.png');
imwrite(flipud(rgb), outPng);                % flipud: PNG is top-down
fprintf('  wrote %s\n', outPng);

%% An empty view
% A view that covers nothing yields a valid 0-by-0 image rather than a stale
% raster from a previous call.
fprintf('\nAn empty viewDir\n-------------\n');
flatV = [0 0 0; 1 0 0; 1 1 0];
flatF = [1 2 3];
flatMesh = Flatland(flatV, flatF);
edgeOn = flatMesh.render([1 0 0], 'Resolution', 1e-2);   % exactly edge-on
fprintf('  edge-on raster: %dx%d, mask numel %d\n', ...
        edgeOn.width, edgeOn.height, numel(edgeOn.mask));
fprintf('  coveredPixels : %d, hasStats: %d\n', ...
        edgeOn.result.coveredPixels, edgeOn.result.hasStats);
delete(flatMesh);

delete(fl);
