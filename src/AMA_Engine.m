classdef AMA_Engine < handle
%fMM X
%FLL X
%GSS X
%GRD
%SGD
% XXX bFourier in P
% XXX bAnalytic & complex
%fXXX bPhaseParam
properties
    objType
    errType
    bFourier

    S
    N
    P

    f
    w % neuron not stim

    fIndex % filterIndex

    plotLL

    nvar
    Lambda
    ns

% OUTPUT
    r
    R
    Rm

    LAll % [ nStim x nCtg ] @ correct category
    LCor % [ nStim x 1 ]

    PAll
    PCorr

    XCtg
    XHat
    XDev

    E
    Em


end
properties(Hidden)
    expn
    bAbs
    estMeth %  1 median,   2 mean, 3 circ
    devMeth %  1 subtract,         3 circ
    errMeth % -1 negative log, 1 log, 2

    NcmpMax
end
methods()
%- MAIN
    function AMA_Engine(S,N,P,f0,f,fIndex)
        obj.parse(S,N,P);
        if nargin >= 4
            obj.f0=f0;
        end
        if nargin >= 5
            obj.f=f;
        end
        if nargin >= 6
            obj.fIndex=fIndex;
        end
    end
    function Em = objective(obj,typ,f, fminOpts)
        if nargin < 2 || isempty(typ)
            typ=objType;
        else
            obj.objType=typ;
        end

        obj.unpack(f);

        switch typ
        case {'GSS','GMM','FLL','FLLGPU'}
            obj.respond();
            obj.posterior();
            Em=obj.error();
        case 'SGD'
            F=fminOpts;
            P=[obj.nF P obj.nFset,F.btchSz,F.MaxIter,F.stpSzMax,F.stpSzMin,F.stpSzEta];
            Em=amaR01objFuncSGD(obj.f,[],P{:});
        case 'GRD'
            P={obj.S.s, obj.S.ctgInd, obj.S.X, obj.N.rMax, obj.N.fano, obj.N.var0, obj.errorType};
            Em=amaR01objFuncFLLgrd(obj.f,[],P{:});
        case 'GRDGPU'
            % TODO
        end
        obj.Em=Em;
    end
    function [Em,E]=test(obj,S)
        if nargin > 2 && isempty(S)
            obj.parse_stim(S);
        end
        obj.respond();
        obj.posterior();
        [Em,E]=obj.error();
    end
    function XHat=estimate(obj,S);
        if nargin > 2 && isempty(S)
            obj.parse_stim(S);
        end
        obj.respond();
        obj.posterior();
        XHat=obj.getXHat();
    end
end
methods(Hidden)
%- PARSE

    function parse(obj,S,N,P)
        obj.parse_stim(S);
        obj.parse_N(N);
        obj.parwse_P(P);
    end
    function parse_stim(obj,S)
        if ~isefield(s,'sW') || isempty(S.sW) || (size(S.ctgInd,1) ~= size(obj.S.sW,1))
            S.sW=1;
        end
        if ~isfield(S,'ctgInd') || isempty(S.ctgInd);
            S.ctgInd=ones(size(obj.f,1),1);
        end
        if ~isfield(S,'X') || isempty(S.X)
            S.X=1;
        end
        obj.XCtg=obj.S.X(obj.ctgInd)';
        obj.S=S;
    end
    function parse_neuron(obj,N)
        if isfield(N,'w') && ~isempty(N.w)
            obj.w=N.w;
        end
        obj.N=N;
    end
    function parse_P(obj,P)
        errorType='';
        if isstruct(P)
            flds=fieldnames(P);
            fldsE={'estMeth','devMeth','errMeth','bAbs','expn'};
            for i = 1:length(fldsE)
                if ismember(flds{i},flds);
                    errType='custom';
                    obj.(flds{i})=P.(flds{i});
                end
            end
            if isfield(P,'type')
                obj.objType=P.type;
            end
            if isfield(P,'NcmpMax')
                obj.NcmpMax=P.NcmpMax;
            end
        elseif ischar(P)
            errType=P;
        end
        if isempty(errType)
            errorType='MSE2'
        end
        if ~strcmp(errorType,'custom')
            obj.parse_error_type(P.errorType);
        end
    end
    function parse_error_type(obj,errorType)

    % expn
    % bAbs
    % estMeth %0 mode  1 median,   2 mean, 3 circ
    % devMeth %  1 subtract,         3 circ
    % errMeth % -1 norm, 1 log, 2 -log

        % DEFAULTS
        obj.estMeth= 2;
        obj.devMeth= 1;
        obj.errMeth=-1;
        obj.bAbs=true;

        switch errorType
            case 'ME'
                obj.bAbs=false;
                obj.expn=1;
            % MSE BASED
            case {'MSE0','MEM'}
                obj.expn=0;
            case {'MSE1','MAE'}
                obj.expn=1;
            case {'MSE2','MES'}
                obj.expn=2;
            case {'MSEcirc'}
                obj.expn=2;
                obj.estMeth=3;
                obj.devMeth=3;
            % REG BASED
            case {'MLE'}
                % XXX obj.estMeth=0 w/ P
                obj.errMeth=1;
            case {'MAP','ARE'}
                % XXX obj.estMeth=0 w/ L
                obj.errMeth=2;
            case {'MED','LIN'}
                obj.estMeth=1;
                obj.expn=1;
                obj.bAbs=true;
            otherwise
                error(['amaError: WARNING! unhandled error type: ' P.errorType '!!!']);
        end
        if ismember(obj.errorType,{'MED','LIN','LLK'})
            disp(['amaError: WARNING! untested code for errorType: ' P.errorType]);
        end
    end
%- PACK
    function unpack(obj,fIn)
        [obj.f,obj.w]=obj.unpackLrn(fIn);
        obj.nF = size(obj.f,2);
    end
    function [f,w]=unpackLrn(obj,F,fIndex);
        % F(:,:,1) = filters
        % F(:,:,2) = misc params
        if nargin < 2 || isempty(fIndex)
            fIndex=obj.fIndex;
        end

        % filters
        f=F(:,:,1);

        % weights
        if size(F,3) > 1
            w=F(1:obj.nF,1:obj.nF,2);
        else
            w=[];
        end

        if isempty(fIndex)
            return
        end

        % SET FIXED f
        if ~fIndex.skipFFix && ~isempty(f);
            f(fIndex.fFix)=fIndex.f0(fIndex.fFix);
        end
        % SET FIXED w
        if ~fIndex.skipWFix && ~isempty(w);
            w(fIndex.wFix)=fIndex.w0(fIndex.wFix);
        end

        % APPLY W SYMMETRY CONSTRAINTS
        if fIndex.bCmpl
            w=Arr.triuOnTrilCmpl(w);
        elseif fIndex.bSym
            w=Arr.triuOnTril(w);
        end
        if fIndex.bEye
            w=Arr.onEye(w,0);
        end

        if nargout < 1
            obj.f=f;
            obj.w=w;
        end
    end
    function [f,w]=unpackLrn2(obj,F,fIndex)
        if nargin < 2 || isempty(fIndex)
            fIndex=obj.fIndex;
        end

        nFLrn=sum(fIndex.fLrn);
        nWLrn=sum(fIndex.wLrn);
        f=fIndex.f0;
        w=fIndex.w0;
        if nFLrn > 0
            f(fIndex.fLrn)=F(1:fLrn);
        end
        if nWLrn > 0
            w(fIndex.wLrn)=F(fLrn+1:end);
        end
        if fIndex.bCmpl
            w=Arr.triuOnTrilCmpl(w);
        elseif t0.bSym
            w=Arr.triuOnTril(w);
        end
        if nargout < 1
            obj.f=f;
            obj.w=w;
        end
    end
%- RESPOND
    function respond(obj)

        % RESPOND
        if obj.bFourier
            r =(obj.S.FFTs'*obj.f);
        else
            r =(obj.S.s'*obj.f);
        end

        % R NORM
        Rm = RNorm.get(r,obj.f,obj.N.normType,obj.N,obj.S);

        % NOISE
        obj.getNoise();
        if N.bNoise
            nu=obj.ns;
        else
            nu=0;
        end

        % MAXR
        if obj.N.MaximR==1
            Rm=Rm./max(Rm(:))*max(r(:));
        elseif obj.N.MaximR==2
            Rm=R./max(Rm(:));
        elseif obj.N.MaximR==3
            ind=abs(Rm)>1;
            Rm(ind)=1*sign(Rm(ind));
        end

        obj.r  = r;
        obj.R  = obj.N.rMax * (Rm+nu);
        obj.Rm = obj.N.rMax *  Rm;

    end
%- NOISE
    function getNoise(obj)
        obj.nvar=sqrt( obj.N.fano .* abs(obj.r) + obj.N.var0 );
        obj.ns = normrnd(0,obj.nVar,size(obj.r));
    end
    function getLambda(obj)
        % noise variance
        % MATCH CONSTANT ADDITIVE TO AVERAGE SCALED ADDITIVE NOISE VARIANCE
        varAvg     = mean(obj.nvar.^2);
        % INTERNAL FILTER RESPONSE COVARIANCE MATRIX (ASSUMING UNCORRELATED NOISE)
        obj.Lambda     = diag( varAvg.*ones(1,obj.nF) );

        if obj.N.rho==0
            return
        end

        % CORRELATED INTERNAL NEURAL NOISE %%
        C = corr2cov(varAvg,obj.N.rho);
        obj.Lambda(logical(triu(ones(nFobj.), 1))) = C;
        obj.Lambda(logical(tril(ones(obj.nF),-1))) = C;

    end
%- POSTERIOR
    function posterior(obj)
        switch obj.objType
        case 'GMM'
            [obj.PCorr,obj.PAll,obj.LCor,obj.LAll] = AMAengineGMM(obj.R,obj.Rm,obj.lambda,obj.S.ctgInd,obj.S.X,obj.NcmpMax);
        case 'GSS'
            [obj.PCorr,obj.PAll]=obj.posterior_GSS();
        case 'FLL'
                [PCorr,PAll]=        AMAengine(obj.R,obj.Rm,obj.nvar, obj.S.ctgInd);
        case 'FLLGPU'
            try
                [obj.PCorr,obj.PAll]=AMAengineGPU_comp(obj.R,obj.Rm,obj.nvar, obj.S.ctgInd,max(obj.S.ctgInd));
            catch
                [obj.PCorr,obj.PAll]=        AMAengine(obj.R,obj.Rm,obj.nvar,obj.S.ctgInd);
                disp(['amaR01objFuncFLL: WARNING! bGPU = 1, but GPU is not being utilized!!! Quit and type gpuDeviceReset() at command line!']);
            end
        end
    end
    function posterior_GSS(obj)
        obj.getLambda();
        for c = 1:length(obj.S.X)
            bInd       = obj.S.ctgInd == c;
            MU(c,:)    = mean(obj.Rm(bInd,:));
            COV(:,:,c) = cov( obj.Rm(bInd,:)) + obj.Lambda;
            try
                obj.LAll(:,c)  = mvnpdf(R,MU(c,:),COV(:,:,c));
            catch ME
                figure(909)
                colorbar;
                figure(910)
                rr=r(bInd,:);
                rr=r(:);
                hist(rr);
                rethrow(ME);
            end
            % @ correct category
            obj.LCor(bInd,1)  = obj.LAll(bInd,c);
        end
        %obj.LLm = mean(log(LCor));

        N=sum(obj.LAll,2);
        obj.PAll = bsxfun(@rdivide,obj.LAll,N);
        obj.PCor = bsxfun(@rdivide,obj.LCor   ,N); % @ correct category
        % PCorr    [ nStim x 1]
        % PAll [ nStim x nCtg]
        % LCor
        % LAll
        %
        % LLm
    end
%- ERROR
    function [Em,E]=error(obj)
        obj.getXHat();
        obj.getXDev();
        obj.getE();
        if nargout > 0
            Em=obj.Em;
            E=obj.E;
        end
    end
    function Xhat=getXHat(obj)
        switch obj.estMeth
        case 0
            % MODE
            % TODO
        case 1
            % MEDIAN
            obj.XHat  = zeros(size(obj.PAll,1),1);
            for i = 1:size(obj.PAll,1)
                obj.XHat(i,1) = interp1(cumsum(obj.PAll(i,:)),obj.X,0.5,'linear',min(X));
            end
        case 2
            % MEAN
            obj.XHat =  obj.PAll*obj.X';
        case 3
            % CIRC
            obj.XHat=circ_mean(obj.X',obj.PAll');
        end

        if nargout > 0
            XHat=obj.Xhat;
        end
    end
    function getXDev(obj)
        switch obj.devMeth
        case 1
            % SUBTRACTIVE
            obj.XDev=obj.XHat-obj.XCtg;
        case 2
            % DIVISIVE
            obj.XDev=obj.XHat/obj.XCtg;
        case 3
            % CIRC
            obj.XDev=circ_dist(obj.XCtg,obj.XHat);
        end
    end
    function getE(obj)
        switch obj.errMeth
        case -1
            if obj.expn==0
                % TODO
            else
                E=aexp( obj.XDev, obj.bAbs+1, obj.expn);
            end
        case 0
            % L0
            % TODO
        case 1
            % MLE
            E= log(obj.LCor);
        case 2
            % MAP
            E=-log(obj.PCorr);
        end
        obj.E=E .* obj.S.sW;
        obj.Em = mean(obj.E);
    end
end
end
