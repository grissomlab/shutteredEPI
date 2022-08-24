function out = ssq(in,dim)

% sum-of-squares across last non-singleton dim
out = sqrt(sum(abs(in).^2,dim));
