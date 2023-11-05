function [eta_electrolysis, eta_methanation, etaGFU, GCVng, CDFe, CDFg] = GEopfParas()

eta_electrolysis = 0.7;    % from the energy perspective, the effciency is about 80%
eta_methanation = 0.8;
etaGFU = 0.4211;           % from the energy perspective, 从1/200换算而来
GCVng = 41.04 * 1e6;     % J/m3
CDFe = 1e4; % MW/hour, 大概数值，从jia文章中拿的
CDFg = CDFe * GCVng / 3600 / 24;

end