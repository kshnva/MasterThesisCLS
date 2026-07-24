function [S_ref, pial_idx, sink_idx, Qreg, nloops_total, A_mm2, npen] = sens_build(ncombsx, ElementLength, density, seed)
% Build a Phase-2 leptomeningeal network and calibrate Qreg (NoCoupling SS).
% Shared builder for the sensitivity-analysis scripts.
%
% In:  ncombsx (rows of hexagons), ElementLength [m], density [pen/mm^2], seed
% Out: S_ref, pial_idx, sink_idx, Qreg, nloops_total, area [mm^2], n penetrators
%
% Notes: branchlevel=4 -> ncombsy=4; penetrator count = round(density*area).
% Penetrators are frozen at 15 um (kReg=0), same convention as all Phase-2 runs.

    if nargin < 4 || isempty(seed), seed = 1234; end
    if nargin < 3 || isempty(density), density = 3.35; end
    if nargin < 2 || isempty(ElementLength), ElementLength = 2e-3; end
    if nargin < 1 || isempty(ncombsx), ncombsx = 2; end
    rng(seed);

    P_op   = 86;
    r_dead = 5e-6;

    Sin = struct;
    Sin.npenetrator          = 1;
    Sin.ncombsx              = ncombsx;
    Sin.branchlevel          = 4;
    Sin.ElementLength        = ElementLength;
    Sin.balance_Gsin_Gsout   = 10;
    Sin.Vref                 = -0.05;
    Sin.gainRadius2Vref      = 10;
    Sin.ksens                = 0.4;
    Sin.kReg                 = 0.01;
    Sin.CapillaryPotential   = -0.06;
    Sin.CapillaryConductance = 2e-7;
    Sin.kSensQ               = 0.4;
    Sin.gCapBase             = 0.2;
    Sin.ECapBase             = -0.04;
    Sin.DoCapillaryCurrent   = 0;
    Sin.ECresistivity        = 0.1;
    Sin.KeepDynamics         = 0;
    Sin.GsSink               = 1.6e-16;

    S_ref = defaultPars('leptomeningeal_geometry', Sin);

    % penetrator count from THIS geometry's area
    A_hex = (3*sqrt(3)/2) * ElementLength^2;
    n_hex = ncombsx * 2^(4-2);                 % ncombsy = 2^(branchlevel-2) = 4
    A_mm2 = n_hex * A_hex * 1e6;
    npen  = round(density * A_mm2);

    pen = load_penetrator_data('penetrators.csv');
    idx = randi(length(pen.r), npen, 1);
    sp.r = pen.r(idx);  sp.l = pen.l(idx);  sp.R = pen.R(idx);  sp.G = 1./pen.R(idx);
    S_ref = integrate_sampled_penetrators(S_ref, sp);
    n_new = S_ref.nSE - length(S_ref.sources);
    S_ref.sources = [S_ref.sources, zeros(1, n_new)];

    pial_idx = find([S_ref.IE.pial] == 1);
    sink_idx = find(~S_ref.sources);
    nloops_total = sens_cycle_rank(S_ref, pial_idx);

    % calibrate Qreg from a NoCoupling steady state
    [~,~,~,r_cal] = sens_ss(S_ref, 'NoCoupling', 1, 0.1, P_op*133, pial_idx, r_dead, 50, []);
    St = S_ref; [St.IE.r] = vout(r_cal);
    St = calcConductance2025(St); St = solvehemodyn2025(St);
    Qreg = mean(abs([St.SE(sink_idx).Qs]));
end
