% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% Example 2 - Geometry straight from MATLAB arrays, no files involved
%
% The question this answers: *my mesh is already in MATLAB - do I have to write
% it to disk?* No. Flatland() takes vertex and face matrices directly, which is
% the whole reason the C ABI exists. Nothing here touches the filesystem.
%
%     >> cd /path/to/flatland/examples/matlab
%     >> ex02_mesh_from_arrays
%
% Covered:
%   - building a mesh from MATLAB matrices
%   - 1-BASED faces, and the error you get for passing 0-based ones
%   - face winding, and why it decides which surface you measure
%   - sweeping a parametric shape without ever serialising it
%   - reading the geometry back out

clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matlab'));

flatland_load();
fprintf('FlatLand %s\n\n', Flatland.versionString());

%% A mesh from MATLAB matrices
fprintf('A mesh from MATLAB matrices\n---------------------------\n');

[V, F] = tessellatedCylinder(128, 1.0, 2.0);
fl = Flatland(V, F);
fprintf('  built a 128-sided prism in memory: %d vertices, %d triangles\n', ...
        fl.vertexCount, fl.faceCount);

% Viewed down the axis this is a disc of area pi r^2; side-on it is a 2r x h
% rectangle. Both have closed forms, so the numbers are checkable.
endOn  = fl.project([0 0 1], 'Resolution', 1e-3, 'Precision', 'double');
sideOn = fl.project([1 0 0], 'Resolution', 1e-3, 'Precision', 'double');
fprintf('  end-on  %.6f   (pi r^2 = %.6f)\n', endOn.area, pi);
fprintf('  side-on %.6f   (2 r h  = %.6f)\n', sideOn.area, 4.0);

%% Indexing: 1-based here, 0-based in C
% Faces index vertices from 1, as everything in MATLAB does. The binding
% converts to the C API's 0-based indices for you. Do NOT subtract one
% yourself - the binding checks for it and says so.
fprintf('\nIndexing\n--------\n');
try
    Flatland(V, F - 1);                          % 0-based: wrong here
    fprintf('  ERROR: 0-based faces should have been rejected\n');
catch err
    fprintf('  0-based faces rejected, as they should be:\n');
    fprintf('    %s: %s\n', err.identifier, strtrim(err.message));
end

%% Winding decides which surface you measure
% Each triangle's vertices must run counter-clockwise seen from OUTSIDE, so its
% normal points outward. FlatLand culls faces pointing away from the camera, so
% a mesh wound inside-out measures the FAR surface. The areas would look fine;
% any field statistic would be wrong. Here is that failure, made visible.
fprintf('\nWinding\n-------\n');

[BV, BF] = unitBox();
fx = BV(:,1);                                    % field = the x coordinate

flGood = Flatland(BV, BF);
flFlip = Flatland(BV, BF(:, [1 3 2]));           % every triangle reversed

rGood = flGood.project([1 0 0], 'Field', fx, 'Resolution', 2e-3, 'Precision', 'double');
rFlip = flFlip.project([1 0 0], 'Field', fx, 'Resolution', 2e-3, 'Precision', 'double');

fprintf('  correct winding: area %.4f, mean f %.4f  (the near face, x = 0)\n', ...
        rGood.area, rGood.average);
fprintf('  reversed       : area %.4f, mean f %.4f  (the FAR face, x = 1)\n', ...
        rFlip.area, rFlip.average);
fprintf('  the area is identical either way - only the field gives it away\n');
delete(flGood); delete(flFlip);

%% Refining a parametric shape without writing a single file
% The prism converges on a true cylinder as the facet count grows. Watching that
% convergence is painful through a CLI and trivial in process: the mesh never
% leaves memory.
fprintf('\nRefining a parametric shape\n---------------------------\n');
fprintf('  %-8s %-14s %-12s\n', 'sides', 'end-on area', 'rel err vs pi');
for sides = [8 16 32 64 128 256]
    [Vs, Fs] = tessellatedCylinder(sides, 1.0, 2.0);
    m = Flatland(Vs, Fs);
    rs = m.project([0 0 1], 'Resolution', 1e-3, 'Precision', 'double');
    fprintf('  %-8d %-14.7f %.2e\n', sides, rs.area, abs(rs.area - pi) / pi);
    delete(m);
end

%% Reading the geometry back
% vertices() and faces() return what you put in. FlatLand recentres meshes
% internally for numerical conditioning, but that is invisible here, and faces()
% comes back 1-based so the round trip is consistent.
fprintf('\nReading the geometry back\n-------------------------\n');
Vout = fl.vertices();
Fout = fl.faces();
fprintf('  vertices out: %d x %d, faces out: %d x %d\n', ...
        size(Vout, 1), size(Vout, 2), size(Fout, 1), size(Fout, 2));
fprintf('  vertices round trip exactly: %d\n', isequal(size(Vout), size(V)) && max(abs(Vout(:) - V(:))) < 1e-12);
fprintf('  faces round trip exactly   : %d\n', isequal(Fout, F));

delete(fl);

%% ------------------------------------------------------------------------
function [V, F] = tessellatedCylinder(sides, radius, height)
%TESSELLATEDCYLINDER Closed prism approximating a cylinder, axis along +z.
%   Vertices are 1-based in F, wound counter-clockwise seen from outside.
    hz = height / 2;
    a  = 2 * pi * (0:sides-1)' / sides;
    ring = [radius * cos(a), radius * sin(a)];

    V = [ring, -hz * ones(sides, 1);             % bottom ring: rows 1..sides
         ring,  hz * ones(sides, 1)];            % top ring:    rows sides+1..2*sides
    bottomCentre = 2 * sides + 1;
    topCentre    = 2 * sides + 2;
    V = [V; 0 0 -hz; 0 0 hz];

    i = (1:sides)';
    j = mod(i, sides) + 1;                       % next vertex around the ring
    lo_i = i;          lo_j = j;
    hi_i = i + sides;  hi_j = j + sides;

    F = [lo_i, lo_j, hi_j;                       % wall, two triangles per facet
         lo_i, hi_j, hi_i;
         repmat(bottomCentre, sides, 1), lo_j, lo_i;    % bottom cap, normal -z
         repmat(topCentre,    sides, 1), hi_i, hi_j];   % top cap,    normal +z
end

function [V, F] = unitBox()
%UNITBOX Closed unit box spanning [0,1]^3, 1-based faces, outward normals.
    V = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    F = [1 4 3; 1 3 2;   5 6 7; 5 7 8; ...
         1 2 6; 1 6 5;   2 3 7; 2 7 6; ...
         3 4 8; 3 8 7;   4 1 5; 4 5 8];
end
