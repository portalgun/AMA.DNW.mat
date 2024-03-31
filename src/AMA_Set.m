classdef AMA_Set < handle
% XXX children have DEFAULTS as a parameter

%- OPTIONS FOR OTHER THINGSp
properties
    % NRN
    bFourierF
    rMax
    RMaxType
    fano
    normType
    bNoise
    var0
    rho

    %STIM
    sFName % XXX set method in stim
    sAlias % XXX set method in stim
    trainORtestS % XXX set method in stim & base, alias
    measure
    units
end
%- READ ONLY FOR OTHER THINGS
properties
    % NRN
    nF
    f
    w
    %r
    %Rm
    %Lambda
    %ns

    %STIM
    nStim
    nPix
    PszRC
    nDim
    nSplit
    bBino
    X
    %bFourierS
    %XCtg
    %ctgInd
    %S
    %w
    %As

    %TRN
    f0
    w0
end
methods(Static)
end
methods
%-MAIN
    function out=get_set_default(obj,cls,fld)
        out=obj.(cls).DEF.(fld);
    end
%- SPECIAL FIELDS
    function fld=getSFld(obj)
        if strcmp(obj.trainORtest,'train')
            fld='STRN';
        elseif strcmp(obj.trainORtst,'test')
            fld='STST';
        else
            error('invalid');
        end
    end
    function fld=getOFld(obj)
        if strcmp(obj.trainORtest,'train')
            fld='TRN';
        elseif strcmp(obj.trainORtst,'test')
            fld='TST';
        else
            error('invalid');
        end
    end
%- GET/SET
    function out=get_fld(obj,cls,fld)
        if strcmp(cls,'Stim')
            cls=obj.getSFld();
        elseif strcmp(cls,'OUT')
            cls=obj.getOFld();
        end
        if isempty(fld)
            out=obj.(cls);
        else
            out=obj.(cls).(fld);
        end
    end
    function set_fld(obj,cls,fld,val)
        if strcmp(cls,'Stim')
            cls=obj.getSFld();
        elseif strcmp(cls,'OUT')
            cls=obj.getOFld();
        end
        % default
        def=obj.get_set_default(cls,fld);
        if isempty(val)
            val=def;
        end

        % lval
        bLVal=false;
        if ~isfield(obj.last,cls)
            obj.last.(cls)=struct();
            lval=[];
        elseif isfield(obj.last.(cls),fld)
            lval=obj.last.(cls).(fld);
            bLVal=true;
        end

        % cval
        cval=obj.get_fld(cls,fld);

        % set
        obj.(cls).flds=val;

        % STORE
        % XXX OR COPY OBJECTS?
        if ~bLval && obj.bHasFit
            obj.last.(cls).(fld)=val;
        end
        obj.bSet=true;
    end
%- BASE CLASSES
    %- NRN
    function obj=get_nrn_fld(obj,varargin)
        obj.get_fld('Nrn',varargin{:});
    end
    function set_nrn_fld(obj,fld,val)
        obj.setFld('Nrn',fld,val);
    end
    %- STM
    function set_stimfld(obj,fld,val)
        obj.set_fld('Nrn',fld,val);
    end
    function obj=get_stim_fld(obj,varargin)
        obj.get_fld('Stim',varargin{:});
    end
%-Nrn
    %- SET
    function set.bFourierF(obj,val)
        obj.set_nrn_fld('bFourierF');
    end
    function set.rMax(obj,val)
        obj.set_nrn_fld('rMax');
    end
    function set.RMaxType(obj,val)
        obj.set_nrn_fld('RMaxType');
    end
    function set.fano(obj,val)
        obj.set_nrn_fld('fano');
    end
    function set.normType(obj,val)
        obj.set_nrn_fld('normType');
    end
    function set.bNoise(obj,val)
        obj.set_nrn_fld('bNoise');
    end
    function set.var0(obj,val)
        obj.set_nrn_fld('var0');
    end
    function set.rho(obj,val)
        obj.set_nrn_fld('rho')';
    end
    %- GET
    function out=get.bFourierF(obj,val)
        out=obj.get_nrn_fld('bFourierF');
    end
    function out=get.rMax(obj,val)
        out=obj.get_nrn_fld('rMax');
    end
    function out=get.RMaxType(obj,val)
        out=obj.get_nrn_fld('RMaxType');
    end
    function out=get.fano(obj,val)
        out=obj.get_nrn_fld('fano');
    end
    function out=get.normType(obj,val)
        out=obj.get_nrn_fld('normType');
    end
    function out=get.bNoise(obj,val)
        out=obj.get_nrn_fld('bNoise');
    end
    function out=get.var0(obj,val)
        out=obj.get_nrn_fld('var0');
    end
    function out=get.rho(obj,val)
        out=obj.get_nrn_fld('rho')');
    end
    %- READ ONLY
    function out=get.nF(obj)
        out=obj.get_nrn_fld('f');
    end
    function out=get.f(obj)
        out=obj.get_nrn_fld('f');
    end
    function out=get.w(obj)
        out=obj.get_nrn_fld('w');
    end
%-Stim
    %- SET
    function set.sFName(obj,val)
        % NOTE alias
        obj.set_stim_fld('fname',val);
    end
    function set.sAlias(obj,val)
        % NOTE alias
        obj.set_stim_fld('alias',val);
    end
    function set.trainORtestS(obj,val)
        % NOTE alias
        obj.set_stim_fld('trainORtest',val);
    end
    function out=set.measure(obj,val)
        obj.set_stim_fld('measure',val);
    end
    function out=set.units(obj,val)
        obj.set_stim_fld('units',val);
    end
    %- GET
    function out=get.sFName(obj)
        % NOTE alias
        out=obj.get_stim_fld('fname');
    end
    function out=get.sAlias(obj,val)
        % NOTE alias
        out=obj.get_stim_fld('alias');
    end
    function out=get.trainORtestS(obj)
        % NOTE alias
        out=obj.get_stim_fld('trainORtest');
    end
    function out=get.measure(obj)
        out=obj.get_stim_fld('measure');
    end
    function out=get.units(obj)
        out=obj.get_stim_fld('units');
    end
    %- READ ONLY
    function out=get.nStim(obj)
        out=obj.get_stim_fld('nStim');
    end
    function out=get.nPix(obj)
        out=obj.get_stim_fld('nPix');
    end
    function out=get.nPix(obj)
        out=obj.get_stim_fld('PszRC');
    end
    function out=get.nDim(obj)
        out=obj.get_stim_fld('nDim');
    end
    function out=get.nSplit(obj)
        out=obj.get_stim_fld('nSplit');
    end
    function out=get.bBino(obj)
        out=obj.get_stim_fld('bBino');
    end
    function out=get.bBino(obj)
        out=obj.get_stim_fld('X');
    end
%- Obj
    function set.amaType(obj,val)
    end
    function set.meanCon(obj,val)
    end
    function set.fitType(obj,val)
    end
%-TRN
    %- SET
    %- GET
    %- READ ONLY
    function out=get.f0(obj)
        out=obj.get_trn_fld('f0');
    end
    function out=get.w0(obj)
        out=obj.get_trn_fld('w0');
    end
end
end
