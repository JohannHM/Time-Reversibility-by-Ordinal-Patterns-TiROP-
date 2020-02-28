function div = f_divergences_JS_KL(p, q, varargin)
% % % Kulback-Leibler and Jensen-Shannon Divergences.
%
% % Inputs: P and Q as vector of probabilities
% % Output: DIV as a divergence of probabilities
%           Flags: js or kl
%           'js': Jensen-Shannon Divergence between P and Q
%           'kl': Kullback-Leibler Divergence of P respect to Q
%           default (js)
%
% % Kulback-Leibler Divergence (KLD) of probability p respect to probability q.
% % KLD(p|q) = sum[pi*log(pi/qi)]
% % KLD(p|q) = sum[pi*[log(pi) - log(qi)]]
% % This is also known as relative entropy but itsn't a metric or distance since
% % itsn't symmetric: [(KLD(p|q) different from KLD(q|p))] nor the the triangle
% % unnequality is accomplished.
% % ""It is the amount of information lost when Q is used to approximate P.[4] In
% % applications, P typically represents the "true" distribution of data,
% % observations, or a precisely calculated theoretical distribution, while Q
% % typically represents a theory, model, description, or approximation of P.
% % Taken from Wiki""

% % Jensen-Shannon Divergence (JSD) between two probabilities p and q
% % JSD(p|q) = [KLD(p|M) + KLD(q|M)]/2  with    M = (p+q)/2
% % It's also considered as a symmetric divergence, then the root mean squared
% % of this divergence can be considered as a DISTANCE. sqrt(jsd(p,q)).

% JohannM.
% Paris (2017)
% ------------------------------------------------------------------------------

if ~isequal(size(p), size(q))
    error('Fool!... Both inputs must have the same dimension.')
end
if (abs(sum(p) - 1 > 0.00001) || abs(sum(q) - 1 > 0.00001))
    error('Fool!... Probabilities do not sum to 1.')
end

lnp = log(p);
lnp(isinf(lnp)) = eps;  %forcing log(0) elements which are Inf to be MATLAB eps
lnq = log(q);
lnq(isinf(lnq)) = eps;  %forcing log(0) elements which are Inf to be MATLAB eps
resta_PQ = lnp - lnq;   %for kulback-leibler
resta_PQ(resta_PQ <= 1e-10) = 0;

M = (p+q)/2;                %creating the common probability (for jensen-shannon)
lnM = log(M);
lnM(isinf(lnM)) = eps;      %forcing log(0) elements which are Inf to be MATLAB eps
resta_PM = (lnp - lnM);     %substracting log(P) - Log(M)
resta_QM = (lnq - lnM);     %substracting log(Q) - Log(M)
resta_PM(resta_PM <= 1e-10) = 0;  %2 closer probabilities could gives 0.000i (I avoid it)
resta_QM(resta_QM <= 1e-10) = 0;  %2 closer probabilities could gives 0.000i (I avoid it)

if ~isempty(varargin)   %si hay opciones extras haga cosas
    
    switch varargin{1}
        case 'js'
            % jensen-shannon divergence
            jsd = (sum(p.*resta_PM) + sum(q.*resta_QM))/2;
            div = jsd;
        case 'kl'
            % kullback-leibler divergence
            kld = sum(p.*resta_PQ);
            div = kld;
        otherwise
            error(['Last argument' ' " ' varargin{1} ' " ' 'not recognized.']);
    end
    
else                    %si NO hay opciones extra calcule por defecto esto
    % jensen-shannon divergence
    jsd = (sum(p.*resta_PM) + sum(q.*resta_QM))/2;
    div = jsd;
end
end