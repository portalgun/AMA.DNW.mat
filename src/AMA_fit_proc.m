classdef AMA_fit_proc < handle
properties
end
methods
    function cont(obj)
        obj.train('bContinue',true);
    end
    % F
    function new(obj,varargin)
        obj.fitType='F';
        obj.train(varargin{:});
    end
    function recurse(obj,varargin)
        obj.fitType='rF';
        obj.train('bRecurse',true,varargin{:});
    end
    % W
    function newW(obj,varargin)
        obj.fitType='W';
        obj.train('bRecurse',true,'bRecurseW',false,varargin{:});
    end
    function recurseW(obj,varargin)
        obj.fitType='rW';
        obj.train('bRecurse',true,'bRecurseW',true,varargin{:});
    end
    function recurseAllW(obj,varargin)
        obj.fitType='AW';
        obj.train('nWset',obj.nF,'bRecurse',true,'bRecurseW',true,'bAllW',true,'FFix',1:obj.nF,varargin{:});
    end
    function recurseAllOneW(obj,varargin)
        obj.fitType='OW';
        obj.train('nWset',1,'bRecurse',true,'bRecurseW',true,'bAllOneW',true,'FFix',1:obj.nF,varargin{:});
    end
    % AF
    % AF
    function newAll(obj,varargin)
        obj.fitType='AF';
        obj.train('bRecurse',false,'bAll',true,varargin{:});
    end
    function recurseAll(obj,varargin)
        obj.fitType='rAF';
        obj.train('bRecurse',true,'bAll',true,varargin{:});
    end
    % OF
    function newAllOne(obj,varargin)
        obj.fitType='OF';
        obj.train('bRecurse',true,'bAllOne',true,varargin{:});
    end
    function recurseAllOne(obj,varargin)
        obj.fitType='rOF';
        obj.train('bRecurse',true,'bAllOne',true,varargin{:});
    end
    % FW
    function newFW(obj,varargin)
        obj.fitType='FW';
        obj.train('bRecurse',false,'bAll',false,'bRecurseW',false,'bAllW',false,'MaxIter',10000,'MaxFunEvals',1000000,varargin{:});
    end
    function recurseFW(obj,varargin)
        obj.fitType='rFW';
        obj.train('bRecurse',true,'bAll',false,'bRecurseW',true,'bAllW',false,'MaxIter',10000,'MaxFunEvals',1000000,varargin{:});
    end
    % OFW
    function newAllOneFW(obj,varargin)
        obj.fitType='OFW';
        obj.train('bRecurse',false,'bAllOne',true,'bAllOneW',true,10000,'MaxFunEvals',1000000,varargin{:});
    end
    function recurseAllOneFW(obj,varargin)
        obj.fitType='rOFW';
        obj.train('bRecurse',true,'bAllOne',true,'bAllOneW',true,10000,'MaxFunEvals',1000000,varargin{:});
    end
    % AFW
    function newAllFW(obj,varargin)
        obj.fitType='AFW';
        obj.train('bRecurse',false,'bAll',true,'bRecurseW',true,'bAllW',true,'MaxIter',10000,'MaxFunEvals',1000000,varargin{:});
    end
    function recurseAllFW(obj,varargin)
        obj.fitType='rAFW';
        obj.train('bRecurse',true,'bAll',true,'bRecurseW',true,'bAllW',true,'MaxIter',10000,'MaxFunEvals',1000000,varargin{:});
    end
    function newCtg(obj,varargin)
        obj.fitType='C';
        obj.train('bCtg',true,varargin{:});
    end
%- Group modes
    function newRecHyper(obj)
        obj.new();
        obj.recurse('bSkipPrompt',true);
        obj.recurseAllOne('bSkipPrompt',true);
    end
end
end
