% SPDX-License-Identifier: AGPL-3.0-or-later
% Copyright (C) 2026 Interfluo
%
% FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
% for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

%% Example 6 - Checking FlatLand against closed-form answers
%
% The question this answers: *how do I know the number is right?*
%
% The recipe is always the same. Find a geometry whose answer you can write
% down, run FlatLand on it at the resolution you intend to use for real, and
% compare. Then refine the mesh and the pixel size separately, because they are
% different error sources and only one of them is usually worth spending effort
% on.
%
%     >> cd /path/to/flatland/examples/matlab
%     >> ex06_verification
%
% Covered:
%   - projected area of a cube: A(n) = |nx| + |ny| + |nz|
%   - Cauchy's identity <A> = S/4, over many directions at once
%   - blackbody radiant intensity, and Stefan-Boltzmann recovered from it
%   - the phase curve of a sphere in radiative equilibrium
%   - separating rasterization error from mesh error, which is the point
%
% The same closed forms drive the project's own validation suite. See
% docs/VALIDATION.md, and validation/cases.py for the derivations.

clear; clc;
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'matlab'));
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

flatland_load();
fprintf('FlatLand %s\n', Flatland.versionString());

% sigma is exact in the 2019 SI: h, k and c are defined constants, so
% sigma = 2 pi^5 k^4 / (15 h^3 c^2) carries no experimental uncertainty.
SIGMA = 5.670374419184431e-8;                      % W m^-2 K^-4
allOk = true;

%% 1. Projected area of a cube
% For a unit cube the general convex formula collapses to something you can
% check in your head: A(n) = |nx| + |ny| + |nz| for a UNIT direction. Down a
% face normal that is 1; down a face diagonal sqrt(2); down the body diagonal
% sqrt(3), the regular hexagon.
fprintf('\n1. Projected area of a cube\n');
fprintf('---------------------------\n');
printHeader('measured', 'closed form');

[cv, cf] = unitCube();
cube = Flatland(cv, cf);
dirs  = [1 0 0; 1 1 0; 1 1 1; 0.3 0.9 -0.31];
names = {'face normal', 'face diagonal', 'body diagonal (hexagon)', 'generic direction'};
for i = 1:size(dirs, 1)
    u = dirs(i,:) / norm(dirs(i,:));
    r = cube.project(dirs(i,:), 'Resolution', 5e-4, 'Precision', 'double');
    allOk = report(names{i}, r.area, sum(abs(u)), 2e-3) && allOk;
end

%% 2. Cauchy's identity: <A> = S/4
% An integral identity, so it probes many directions at once rather than a few
% hand-picked ones. A cube has S = 6, so the mean projected area is 1.5 - and
% note that no single direction gives 1.5.
fprintf('\n2. Cauchy''s identity: <A> = S/4\n');
fprintf('-------------------------------\n');
printHeader('measured', 'closed form');

sweep = fibonacciDirections(160);
rs    = cube.projectBatch(sweep, 'Resolution', 1e-3, 'Precision', 'double');
areas = [rs.area];
allOk = report('cube, mean over 160 directions', mean(areas), 1.5, 3e-3) && allOk;
fprintf('    (the spread is %.4f .. %.4f - the identity is about the MEAN)\n', ...
        min(areas), max(areas));
delete(cube);

%% 3. Blackbody radiant intensity
% FlatLand integrates over the PROJECTED area, and dA_proj = cos(theta) dA. So
% if the field is a radiance L, the integral IS the radiant intensity:
%
%     I(n) = integral L cos(theta) dA   [W/sr]
%
% A blackbody is a Lambertian emitter with L = sigma T^4 / pi, hence for an
% isothermal body I = L x A_proj, and for a sphere the pi cancels:
%
%     I = (sigma T^4 / pi)(pi r^2) = sigma T^4 r^2,  the same from every side.
fprintf('\n3. Blackbody radiant intensity\n');
fprintf('------------------------------\n');

T        = 1200.0;
radiance = SIGMA * T^4 / pi;
sphere   = Flatland.load(fullfile(root, 'examples', 'sphere_areas', 'sphere_fine.obj'));
V        = sphere.vertices();
F        = sphere.faces();
field    = repmat(radiance, size(V,1), 1);
exactI   = SIGMA * T^4;                            % r = 1

fprintf('  isothermal sphere at %.0f K: L = %.2f W/m2/sr\n\n', T, radiance);
printHeader('FlatLand', 'sigma T^4 r^2');
for d = [0 0 1; 1 1 1; -0.4 0.9 0.2]'
    r = sphere.project(d', 'Field', field, 'Resolution', 1e-3, 'Precision', 'double');
    allOk = report(sprintf('I along (%.1f, %.1f, %.1f)', d), r.integral, exactI, 2e-3) && allOk;
end

% The MEAN is a different kind of check. A constant field must interpolate to
% exactly that constant at every covered pixel, whatever shape those pixels
% cover, so the mean carries no mesh error at all - it comes back to ten digits,
% while the integral above is held to a few parts in 1e4 by the inscribed mesh's
% projected-area deficit.
r = sphere.project([0 0 1], 'Field', field, 'Resolution', 1e-3, 'Precision', 'double');
allOk = report('mean radiance (interpolation)', r.average, radiance, 1e-9) && allOk;

%% 4. Stefan-Boltzmann, recovered from projected areas
% Integrate the intensity over all directions and apply Cauchy:
%
%     integral I dOmega = (sigma T^4/pi)(4 pi)(S/4) = sigma T^4 S
%
% The left side never uses the surface area. Getting sigma T^4 S back out of a
% pile of projected areas is the check - it is the Stefan-Boltzmann law,
% reassembled from geometry.
fprintf('\n4. Stefan-Boltzmann, recovered from projected areas\n');
fprintf('---------------------------------------------------\n');
printHeader('4 pi <I> [W]', 'sigma T^4 S');

sweep = fibonacciDirections(300);
rs    = sphere.projectBatch(sweep, 'FieldMatrix', field, ...
                            'Resolution', 2e-3, 'Precision', 'double');
power   = 4*pi*mean([rs.integral]);
C       = faceCross(V, F);
surface = 0.5 * sum(sqrt(sum(C.^2, 2)));
allOk = report('sphere_fine, 300 directions', power, SIGMA*T^4*surface, 5e-3) && allOk;
fprintf('    (mesh surface area %.6f, vs 4 pi = %.6f for a true sphere)\n', ...
        surface, 4*pi);

%% 5. A sphere in radiative equilibrium - the phase curve
% Absorbed flux goes as the incidence cosine, so balancing it against sigma T^4
% gives the subsolar law T(psi) = T_sub cos^{1/4}(psi) on the lit side and
% nothing beyond the terminator. The disc-integrated result is the Lambert-sphere
% phase function (Russell 1916):
%
%     I(alpha)   = (2/3) sigma T_sub^4 r^2 Phi(alpha)
%     Phi(alpha) = [sin(alpha) + (pi - alpha) cos(alpha)] / pi
%
% Two references are shown, and the gap between them is the whole lesson: the
% middle column is exact FOR THIS MESH, so the "raster" error is the tool's own;
% "vs sphere" adds the mesh's fidelity to a true sphere on top.
fprintf('\n5. A sphere in radiative equilibrium - the phase curve\n');
fprintf('------------------------------------------------------\n');

src  = [0 0 1];
vals = radiance * max(0, (V ./ sqrt(sum(V.^2, 2))) * src(:));
base = (2/3) * SIGMA * T^4;

fprintf('  sphere_fine, %d triangles, %.0f K subsolar\n\n', size(F,1), T);
fprintf('  %6s %14s %14s %14s %10s %10s\n', ...
        'phase', 'I FlatLand', 'exact (mesh)', 'exact (sphere)', 'raster', 'vs sphere');
fprintf('  %6s %14s %14s %14s %10s %10s\n', ...
        rep('-',6), rep('-',14), rep('-',14), rep('-',14), rep('-',10), rep('-',10));
for deg = [0 30 60 90 120]
    a = deg * pi/180;
    % The observer sits at angle alpha from the source; the camera looks the
    % other way, back at the body.
    viewDir   = -[sin(a) 0 cos(a)];
    r         = sphere.project(viewDir, 'Field', vals, 'Resolution', 1e-3, 'Precision', 'double');
    meshExact = fieldIntegralConvex(V, F, vals, viewDir);
    smooth    = base * (sin(a) + (pi - a)*cos(a)) / pi;
    raster    = abs(r.integral - meshExact) / meshExact;
    fprintf('  %5d%s %14.4f %14.4f %14.4f %10.2e %10.2e\n', deg, char(176), ...
            r.integral, meshExact, smooth, raster, abs(r.integral - smooth)/smooth);
    allOk = allOk && raster <= 2e-3;
end
delete(sphere);

%% 6. Which error should you spend effort on?
% Refine the mesh and the pixel size independently. Only one of them is usually
% worth the money, and this table tells you which.
fprintf('\n6. Which error should you spend effort on?\n');
fprintf('------------------------------------------\n');
fprintf('  %-18s %8s %8s %14s %12s %12s\n', ...
        'mesh', 'tris', 'res', 'area', 'vs mesh', 'vs pi');
fprintf('  %-18s %8s %8s %14s %12s %12s\n', ...
        rep('-',18), rep('-',8), rep('-',8), rep('-',14), rep('-',12), rep('-',12));
probe = [1 0.3 0.2];
for name = {'sphere_coarse', 'sphere_fine'}
    m  = Flatland.load(fullfile(root, 'examples', 'sphere_areas', [name{1} '.obj']));
    mv = m.vertices(); mf = m.faces();
    exact = projectedAreaConvex(mv, mf, probe);
    for res = [4e-3 1e-3]
        r = m.project(probe, 'Resolution', res, 'Precision', 'double');
        fprintf('  %-18s %8d %8.0e %14.7f %12.2e %12.2e\n', name{1}, size(mf,1), ...
                res, r.area, abs(r.area - exact)/exact, abs(r.area - pi)/pi);
    end
    delete(m);
end

fprintf(['\n  Read the last two columns against each other.\n\n' ...
         '  "vs mesh" is FlatLand against the exact answer for the polyhedron it was\n' ...
         '  actually handed, so it is the tool''s own error. It stays small and bounces\n' ...
         '  around rather than falling cleanly: coverage is decided by point-sampling\n' ...
         '  pixel centres, and the boundary band that error lives in is not a smooth\n' ...
         '  function of the pixel size. Expect a trend, not a rate.\n\n' ...
         '  "vs pi" is the same runs against a true sphere, so it adds the mesh''s own\n' ...
         '  fidelity, and it falls sharply with subdivision until it reaches the floor\n' ...
         '  the first column sets.\n\n' ...
         '  Run this table on your own geometry before deciding whether to add\n' ...
         '  triangles or shrink pixels - guessing usually picks wrong.\n\n']);

if allOk
    fprintf('every closed-form check agreed within tolerance.\n');
else
    fprintf('SOME CHECKS WERE OUTSIDE TOLERANCE - see the OFF rows above.\n');
end

%% ------------------------------------------------------------------------
% Helpers. Local functions in a script need R2016b or later.
%% ------------------------------------------------------------------------

function s = rep(c, n)
s = repmat(c, 1, n);
end

function printHeader(a, b)
fprintf('  %-34s %14s %14s %10s\n', 'case', a, b, 'rel err');
fprintf('  %-34s %14s %14s %10s\n', rep('-',34), rep('-',14), rep('-',14), rep('-',10));
end

function ok = report(label, measured, exact, tol)
rel = abs(measured - exact) / abs(exact);
ok  = rel <= tol;
verdict = 'OFF';
if ok, verdict = 'ok'; end
fprintf('  %-34s %14.6f %14.6f %10.2e  %s\n', label, measured, exact, rel, verdict);
end

function [V, F] = unitCube()
[x, y, z] = ndgrid([-0.5 0.5], [-0.5 0.5], [-0.5 0.5]);
V = [x(:) y(:) z(:)];
% 1-BASED faces, wound counter-clockwise seen from outside. Indices follow
% ndgrid order: the fastest-varying subscript is x.
F = [1 3 4; 1 4 2; 5 6 8; 5 8 7; 1 2 6; 1 6 5; ...
     3 7 8; 3 8 4; 1 5 7; 1 7 3; 2 4 8; 2 8 6];
end

function C = faceCross(V, F)
%FACECROSS  (v1-v0) x (v2-v0) per face. Twice the area, along the outward normal.
C = cross(V(F(:,2),:) - V(F(:,1),:), V(F(:,3),:) - V(F(:,1),:), 2);
end

function A = projectedAreaConvex(V, F, direction)
%PROJECTEDAREACONVEX  A(n) = 1/4 sum |n . c_i| - exact for a closed convex mesh.
%
%   Every line through a convex body crosses the surface twice, so each
%   direction sees half the total projected face area. This is THE reference to
%   compare against: the exact answer for the mesh you actually handed over,
%   with no reference to whatever smooth shape it approximates.
u = direction(:) / norm(direction);
A = 0.25 * sum(abs(faceCross(V, F) * u));
end

function I = fieldIntegralConvex(V, F, values, direction)
%FIELDINTEGRALCONVEX  Exact area integral of a node field over the visible side.
%
%   Barycentric interpolation reproduces a linear function exactly and the mean
%   of a linear function over a triangle is the mean of its vertex values, so a
%   front-facing triangle contributes (its projected area) x (that mean).
u    = direction(:) / norm(direction);
C    = faceCross(V, F);
proj = C * u;
vis  = proj < 0;                         % back-facing: the camera looks ALONG u
% The reshape is not decoration. Indexing a vector with a MATRIX gives the
% shape of the index, but indexing a vector with a VECTOR gives the
% orientation of the vector being indexed - so when exactly one face is
% visible, values(F(vis,:)) comes back 3-by-1 instead of 1-by-3 and the mean
% would be taken along the wrong dimension.
means = mean(reshape(values(F(vis,:)), [], 3), 2);
I = sum(0.5 * abs(proj(vis)) .* means);
end

function D = fibonacciDirections(count)
%FIBONACCIDIRECTIONS  Quasi-uniform directions on the sphere.
%
%   Deterministic, and converges far faster than random sampling for the
%   direction averages above.
i  = (0:count-1)';
ga = pi * (3 - sqrt(5));
z  = 1 - (2*i + 1)/count;
r  = sqrt(max(0, 1 - z.^2));
D  = [r.*cos(i*ga), r.*sin(i*ga), z];
end
