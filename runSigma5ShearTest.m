function runSigma5ShearTest(direction)
% direction = +1 (default), -1 (reverse shear), or 0 (unstrained control).
if nargin < 1, direction = 1; end
validateattributes(direction, {'numeric'}, {'scalar','real','finite'});
assert(ismember(direction,[-1 0 1]), 'Use direction -1, 0, or +1.');
params = makeSigma5ShearParams(direction);
[~, tag] = fileparts(tempname);
params.output_dir = fullfile(pwd, sprintf('sigma5_shear_%+d_%s',direction,tag));
fprintf('Sigma5: +/- %.6f deg, L=%.6g, dx=%.6g\n', ...
    params.theta_1,L,params.dx);
fprintf('Output directory: %s\n',params.output_dir);
simulateXPFC(params);
end
