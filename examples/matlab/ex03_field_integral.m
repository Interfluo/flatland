% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% Example 3 - The area integral of a scalar field over the visible projection
%
% The question this answers: *I have a quantity defined on my mesh. What is its
% integral over the part a given direction can actually see?*
%
%     I = integral f dA = sum over covered pixels of (value x pixel_area)
%
% FlatLand attaches no physical meaning to f. Supply a radiance and I is a
% radiant intensity; supply a pressure and it is a force; supply an
% emissivity-weighted temperature and it is a thermal signature. The arithmetic
% is the same and the interpretation is yours.
%
%     >> cd /path/to/flatland/examples/matlab
%     >> ex03_field_integral
%
% Covered:
%   - node fields (interpolated) and face fields (piecewise constant)
%   - the Lambertian sphere, checked against its closed form
%   - which surface the statistics come from, and why that matters
%   - hasStats: statistics that do not exist, reported as NaN rather than 0

clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matlab'));

flatland_load();
fprintf('FlatLand %s\n\n', Flatland.versionString());

%% A node field: the Lambertian sphere
% f = cos(angle between the surface normal and the view direction). For a sphere
% centred at the origin the outward normal at a vertex is just the vertex
% direction, so f = -vhat . n on the visible side.
%
% This has a closed form. At projected radius rho the cosine is
% sqrt(1 - (rho/r)^2), so
%
%   integral f dA = int_0^r sqrt(1-(rho/r)^2) 2 pi rho drho = (2/3) pi r^2
%   mean          = 2/3
%
% Read f as a radiance and that integral is the radiant intensity of a
% Lambertian sphere. FlatLand does not know or care.
fprintf('A node field: the Lambertian sphere\n');
fprintf('-----------------------------------\n');

viewDir       = [0 0 1];
exactIntegral = (2/3) * pi;
exactMean     = 2/3;

fprintf('  %-8s %-8s %-13s %-13s %-10s\n', 'subdiv', 'tris', 'integral', 'exact', 'rel err');
for k = 1:4
    [V, F] = icosphere(k, 1.0);
    % One field value per VERTEX, in mesh vertex order.
    field = -(V * viewDir(:)) ./ sqrt(sum(V.^2, 2));

    m = Flatland(V, F);
    r = m.project(viewDir, 'Field', field, 'Resolution', 2e-3, 'Precision', 'double');
    delete(m);

    rel = abs(r.integral - exactIntegral) / exactIntegral;
    fprintf('  %-8d %-8d %-13.7f %-13.7f %.2e\n', k, size(F,1), r.integral, exactIntegral, rel);
end
fprintf('  mean over the disc: %.6f   (exact 2/3 = %.6f)\n', r.average, exactMean);
fprintf('  min %.4f at the limb, max %.4f at the sub-observer point\n', r.min, r.max);

%% A face field: constant per triangle
% A face field has one value per TRIANGLE and is constant across it - right for
% something measured per facet, like a per-panel emissivity or a per-element CFD
% cell value. The length decides which kind it is; when a mesh has as many faces
% as vertices, pass 'FieldMode' explicitly.
fprintf('\nA face field: constant per triangle\n');
fprintf('-----------------------------------\n');

[V, F] = icosphere(2, 1.0);
faceField = 2.5 * ones(size(F, 1), 1);
m = Flatland(V, F);
r = m.project(viewDir, 'Field', faceField, 'Resolution', 2e-3, 'Precision', 'double');
fprintf('  %d faces, all carrying 2.5\n', size(F,1));
fprintf('  mean     %.6f  (a constant field must come back exactly)\n', r.average);
fprintf('  integral %.6f  = 2.5 x area %.6f = %.6f\n', r.integral, r.area, 2.5 * r.area);
delete(m);

%% Which surface do the statistics describe?
% The visible one - the surface facing the camera. A field that varies through
% the object makes that concrete. On a unit box with f = x, looking along +x
% means the near face is x = 0, so the mean must be 0.
fprintf('\nWhich surface do the statistics describe?\n');
fprintf('-----------------------------------------\n');

[BV, BF] = unitBox();
fx  = BV(:,1);
box = Flatland(BV, BF);
near = box.project([ 1 0 0], 'Field', fx, 'Resolution', 2e-3, 'Precision', 'double');
far  = box.project([-1 0 0], 'Field', fx, 'Resolution', 2e-3, 'Precision', 'double');
fprintf('  looking along +x: mean f = %.4f   (the x = 0 face)\n', near.average);
fprintf('  looking along -x: mean f = %.4f   (the x = 1 face)\n', far.average);
fprintf('  the areas are identical (%.4f), only the field differs\n', near.area);
delete(box);

%% When there are no statistics to report
% A view can legitimately cover nothing: geometry thinner than a pixel, or a
% mesh seen exactly edge-on. The four field statistics are then ABSENT, and this
% binding reports them as NaN. Zero is an ordinary field value, so passing it
% through would be indistinguishable from a measurement and would quietly drag
% any mean or minimum you compute toward it.
fprintf('\nWhen there are no statistics to report\n');
fprintf('--------------------------------------\n');

sliverV = [0 0 0; 1 0 0; 1 1e-6 0];
sliverF = [1 2 3];
sliver  = Flatland(sliverV, sliverF);
r = sliver.project([0 0 -1], 'Field', [-5; -3; -1], 'Resolution', 0.1, 'Cull', false);
fprintf('  coveredPixels : %d\n',    r.coveredPixels);
fprintf('  hasField      : %d   (a field WAS supplied)\n', r.hasField);
fprintf('  hasStats      : %d   (but nothing was measured)\n', r.hasStats);
fprintf('  average       : %g\n',    r.average);
fprintf('  area          : %g   (still a real measurement, so not NaN)\n', r.area);
delete(sliver);

fprintf('\n  When aggregating a batch, guard on hasStats:\n');
fprintf('    good = [R.hasStats];\n');
fprintf('    m    = mean([R(good).average]);   %% or mean(..., ''omitnan'')\n');

%% ------------------------------------------------------------------------
function [V, F] = icosphere(subdivisions, radius)
%ICOSPHERE Icosahedron subdivided and projected onto a sphere. 1-based faces.
    t = (1 + sqrt(5)) / 2;
    V = [-1  t  0;  1  t  0; -1 -t  0;  1 -t  0;
          0 -1  t;  0  1  t;  0 -1 -t;  0  1 -t;
          t  0 -1;  t  0  1; -t  0 -1; -t  0  1];
    V = V ./ sqrt(sum(V.^2, 2));
    F = [ 1 12  6;  1  6  2;  1  2  8;  1  8 11;  1 11 12;
          2  6 10;  6 12  5; 12 11  3; 11  8  7;  8  2  9;
          4 10  5;  4  5  3;  4  3  7;  4  7  9;  4  9 10;
          5 10  6;  3  5 12;  7  3 11;  9  7  8; 10  9  2];

    for s = 1:subdivisions
        cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
        newF  = zeros(size(F,1) * 4, 3);
        n     = 0;
        for f = 1:size(F,1)
            a = F(f,1); b = F(f,2); c = F(f,3);
            [ab, V, cache] = midpoint(a, b, V, cache);
            [bc, V, cache] = midpoint(b, c, V, cache);
            [ca, V, cache] = midpoint(c, a, V, cache);
            newF(n+1:n+4, :) = [a ab ca; b bc ab; c ca bc; ab bc ca];
            n = n + 4;
        end
        F = newF;
    end
    V = V * radius;
end

function [idx, V, cache] = midpoint(a, b, V, cache)
%MIDPOINT Index of the vertex halfway between a and b, pushed onto the sphere.
    key = sprintf('%d_%d', min(a,b), max(a,b));
    if isKey(cache, key)
        idx = cache(key);
        return
    end
    m = (V(a,:) + V(b,:)) / 2;
    V(end+1, :) = m / norm(m);
    idx = size(V, 1);
    cache(key) = idx;
end

function [V, F] = unitBox()
%UNITBOX Closed unit box spanning [0,1]^3, 1-based faces, outward normals.
    V = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    F = [1 4 3; 1 3 2;   5 6 7; 5 7 8; ...
         1 2 6; 1 6 5;   2 3 7; 2 7 6; ...
         3 4 8; 3 8 7;   4 1 5; 4 5 8];
end
