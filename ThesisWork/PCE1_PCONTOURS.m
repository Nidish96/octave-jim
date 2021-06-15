function [z1s, z2s, Zsamps] = PCE1_PCONTOURS(zxps, perc, pdist, lims, Nsamp)
%PCE1_PCONTOURS returns the Percentile contours of the set of PCEs
%
%   USAGE:
%       [z1s, z2s] = PCE1_CONTOURS(zxps, alpha, pdist, lims, Nsamp)
%   INPUTS:
%       zxps
%       alpha
%       pdist
%       lims
%       Nsamp
%   Outputs
%       z1s
%       z2s
    z1s = zeros(size(zxps(:,1)));
    z2s = zeros(size(zxps(:,1)));
    
    if Nsamp~=0
        Zsamps = zeros(length(z1s), Nsamp);
    end
    for iq=1:length(z1s)
%         p1 = (1-perc)/2;
%         p2 = (1+perc)/2;
        p1 = perc;
        p2 = 1-perc;
        switch Nsamp
            case 0  % Analytical
                z10 = polyval(zxps(iq,:), icdf(pdist, p1));
                dzdx10 = polyval(polyder(zxps(iq,:)), icdf(pdist, p1));
                if sign(dzdx10)==-1
                    z20 = z10;
                    z10 = polyval(zxps(iq,:), icdf(pdist, p2));
                else
                    z20 = polyval(zxps(iq,:), icdf(pdist, p1));
                end    
                
                z1s(iq) = fzero(@(z) POLYCDF(zxps(iq,:), z, @(z) cdf(pdist, z), ...
                    @(z) pdf(pdist, z), lims, p1), z10);
                z2s(iq) = fzero(@(z) POLYCDF(zxps(iq,:), z, @(z) cdf(pdist, z), ...
                    @(z) pdf(pdist, z), lims, p2), z20);
                Zsamps = [];
            otherwise  % Sampling
               rx = random(pdist, Nsamp, 1);
               rz = polyval(zxps(iq,:), rx);

               [f, z] = ecdf(rz);
               [~, ui] = unique(z);

               z1s(iq) = interp1(f(ui), z(ui), p1);
               z2s(iq) = interp1(f(ui), z(ui), p2);
               
               Zsamps(iq, :) = rz;
        end
        
        fprintf('%d\n', iq);
    end
end