function mlderiv = MLmax4Dline2(av, datameas, covmats205)

% fit 204/205, 206/205, 207/205, 208/205 data to line
% no overdispersion calc

a = [0 av(1) av(2) av(3)]';  %intercept with 204/205 = 0 plane
v = [1 av(4) av(5) av(6)]';  %direction vector with 206/204 207/204 208/204
n = size(covmats205,3);

mlderiv = 0;
for i = 1:n
    ri = datameas(:,i);  %measured ratios i  (a blank IC run)
    xi = ri - a;         %used often
    si = covmats205(:,:,i);
    dLdai = -2*(ri-a)'*(inv(si) - (si\v*v'/si)/(v'/si*v));
    dLdvi = -2*(v'/si*v * v'/si*xi * (si\xi) - (v'/si*xi)^2*(si\v))/((v'/si*v)^2);
%    sinv = inv(si);
%    dLdvi2 = -2*(v'*sinv*v * v'*sinv*xi *sinv*xi - (v'*sinv*xi)^2*sinv*v)/((v'*sinv*v)^2);
%    dLdai_short = dLdai(1:4)';    %since we know a(1)
%    dLdvi_short = dLdvi(1:4);    %likewise with v
    mlderiv = mlderiv + [dLdai'; dLdvi];
end
