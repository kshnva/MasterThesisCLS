% smoke_test.m
%
% Quick end-to-end check that the whole pipeline runs: geometry build,
% penetrator sampling, electrical network, hemodynamic + electrical solves,
% and the adaptation ODE. Uses the smallest network and one coupling model,
% so it finishes in well under a minute. Not a scientific run.

PROJ = fileparts(mfilename('fullpath'));
addpath(PROJ); addpath(fullfile(PROJ,'ElectricModel'));

fprintf('Building a small network (ncombsx=1)...\n');
[S_ref, pial_idx, ~, Qreg, nloops_total] = sens_build(1, 2e-3, 3.35, 1234);
fprintf('  built: %d pial segments, %d possible loops, Qreg=%.3g\n', ...
    numel(pial_idx), nloops_total, Qreg);

fprintf('Running ArtCoupling to steady state (cap 50 chunks)...\n');
[nl, nsurv, conv] = sens_ss(S_ref, 'ArtCoupling', Qreg, 0.1, 86*133, ...
    pial_idx, 5e-6, 50, []);
fprintf('  result: %d/%d loops survived, %d segments alive, converged=%d\n', ...
    nl, nloops_total, nsurv, conv);

if conv
    fprintf('\nSMOKE TEST PASSED: build + solve + adapt all ran and converged.\n');
else
    fprintf('\nRan end-to-end without error (did not converge within 50 chunks, still fine for a smoke test).\n');
end
