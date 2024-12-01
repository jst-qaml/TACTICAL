function sps_score = calSpsScore(ap, np, af, nf, sps_metric)
% spsCalculator function returns the suspiciousness score based on ap, np, af, nf and sps_metric.
%
% Insuts:
%   ap: successed and activated
%   np: successed and unactivated
%   af: failed and activated
%   nf: failed and unactivated
%   sps_metric: suspiciousness metrics, including alpha, beta and minus (alpha-beta) version
% Outputs:
%   sps_score: suspiciousness score

if strcmp(sps_metric, 'Tarantula')
    sps_score = (af./(af+nf)) ./ ((af./(af+nf)) + (ap./(ap+np)));
elseif strcmp(sps_metric, 'Ochiai')
    sps_score = af ./ sqrt((af+nf) .* (af+ap));
elseif strcmp(sps_metric, 'D2')
    sps_score = (af.^2) ./ (ap+nf);
elseif strcmp(sps_metric, 'D3')
    sps_score = (af.^3) ./ (ap+nf);
elseif strcmp(sps_metric, 'Jaccard')
    sps_score = af ./ (af+nf+ap);
elseif strcmp(sps_metric, 'Kulczynski1')
    sps_score = af ./ (nf+ap);
elseif strcmp(sps_metric, 'Kulczynski2')
    sps_score = 0.5 * (af./(af+nf) + af./(af+ap));
elseif strcmp(sps_metric, 'Op2')
    sps_score = af - (ap ./ (ap+np+1));
% elseif strcmp(sps_metric, 'betaTarantula')
%     sps_score = (nf/(nf+af))/((nf/(nf+af))+(np/(np+ap)));
% elseif strcmp(sps_metric, 'betaOchiai')
%     sps_score = nf/sqrt((nf+af)*(nf+np));
% elseif strcmp(sps_metric, 'betaDstar2')
%     sps_score = (nf^2)/(np+af);
% elseif strcmp(sps_metric, 'betaDstar3')
%     sps_score = (nf^3)/(np+af);
% elseif strcmp(sps_metric, 'betaJaccard')
%     sps_score = nf/(nf+af+np);
% elseif strcmp(sps_metric, 'betaKulczynski1')
%     sps_score = nf/(af+np);
% elseif strcmp(sps_metric, 'betaKulczynski2')
%     sps_score = 0.5 * (nf/(nf+af) + nf/(nf+np));
% elseif strcmp(sps_metric, 'betaOp2')
%     sps_score = nf-(np/(np+ap+1));
% elseif strcmp(sps_metric, 'minusTarantula')
%     sps_score = (af/(af+nf))/((af/(af+nf))+(ap/(ap+np))) - (nf/(nf+af))/((nf/(nf+af))+(np/(np+ap)));
% elseif strcmp(sps_metric, 'minusOchiai')
%     sps_score = af/sqrt((af+nf)*(af+ap)) - nf/sqrt((nf+af)*(nf+np));
% elseif strcmp(sps_metric, 'minusDstar2')
%     sps_score = (af^2)/(ap+nf) - (nf^2)/(np+af);
% elseif strcmp(sps_metric, 'minusDstar3')
%     sps_score = (af^3)/(ap+nf) - (nf^3)/(np+af);
% elseif strcmp(sps_metric, 'minusJaccard')
%     sps_score = af/(af+nf+ap) - nf/(nf+af+np);
% elseif strcmp(sps_metric, 'minusKulczynski1')
%     sps_score = af/(nf+ap) - nf/(af+np);
% elseif strcmp(sps_metric, 'minusKulczynski2')
%     sps_score = 0.5 * (af/(af+nf) + af/(af+ap)) - 0.5 * (nf/(nf+af) + nf/(nf+np));
% elseif strcmp(sps_metric, 'minusOp2')
%     sps_score = af-(ap/(ap+np+1)) - (nf-(np/(np+ap+1)));
else
    error('Check your suspiciousness metrics!');
end
end