classdef AMAObjOpts < handle & AMAChld & AMAOpts
properties
    ppAlg
    errType

    expn
    bAbs
    estMeth
    devMeth
    errMeth

    regularization
    rlambda
end
properties(Hidden)
    bSetError=false
    bSettingErr=false
    AMADeps={''}
end
methods(Static)
    function P=getP()
        % XXX defaults
        P={...
            'ppAlg','GSS','';
            'errType','','';
            ...
            'expn',[],'';
            'bAbs',[],'';
            'estMeth',[],'';
            'devMeth',[],'';
            'errMeth',[],'';
            'regularization',[],'';
            'rlambda',0.1,'';
        };
    end
end
methods
    function obj=AMAObjOpts()
        if nargin < 1
            return
        end
    end
    function init(obj)
        fldsE={'estMeth','devMeth','errMeth','bAbs','expn'};
        if isempty(obj.errType)
            for i = 1:length(fldsE)
                fld=fldsE{i};
                if ~isempty(obj.(fld))
                    errType='custom';
                    break
                end
            end
        end
        if isempty(obj.errType)
            obj.setErr('KLP');
        end
        obj.setErr();
    end
    %function checkErrType(obj)
    %end
    function set.errType(obj,errType)
        if obj.bSettingErr
            obj.errType=errType;
        else
            obj.setErr(errType);
        end
    end
    function setErr(obj,errType)
    % expn
    % bAbs
    % estMeth %0 mode  1 median,   2 mean, 3 circ
    % devMeth %  1 subtract,         3 circ
    % errMeth % -1 norm, 1 -log LL, 2 -log PP 3 LL 4 PP


        obj.bSettingErr=true;
        if nargin < 2 || isempty(errType)
            if isempty(obj.errType)
                obj.errType='KLP';
                errType=obj.errType;
            else
                errType=obj.errType;
            end
        else
            obj.errType=errType;
        end
        obj.bSettingErr=false;


        % DEFAULTS
        obj.estMeth= 2;
        obj.devMeth= 1;
        obj.errMeth=-1;
        obj.bAbs=true;

        bSetError=true;
        switch errType
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
            case {'KLL'}
                % TODO obj.estMeth=0 w/ L
                obj.errMeth=1;
            case {'KLP'}
                % TODO obj.estMeth=0 w/ L
                obj.errMeth=2;
            case {'MLE'}
                % TODO obj.estMeth=0 w/ P
                obj.errMeth=3;
            case {'MAP','ARE'}
                % TODO obj.estMeth=0 w/ L
                obj.errMeth=4;
            case {'MED','LIN'}
                obj.estMeth=1;
                obj.expn=1;
                obj.bAbs=true;
            case 'custom'
                ;
            otherwise
                bSetError=false;
        end
        obj.bSetError=bSetError;
        obj.estMeth=2;
    end
    function checkErrType(obj)
        if isempty(obj.errType)
            obj.errType='KLP';
        end
    end
    function assert(obj)
        if strcmp(obj.ppAlg,'GMM')
            Error.warnSoft('amaR01: WARNING! amaType = GMM RUNS VERY SLOWLY. EXAMINE AMAengineGMM.m AND mixgaussMLEfit.m TO ATTEMPT SPEED IMPROVEMENTS!!!');
        end
        if ~obj.bSetError
            error(['amaError: WARNING! unhandled error type: ' obj.errType '!!!']);
        elseif ismember(obj.errType,{'MED','LIN','LLK'})
            disp(['amaError: WARNING! untested code for errType: ' obj.errType]);
        end
    end
end
end
