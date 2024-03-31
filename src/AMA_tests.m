classdef AMA_tests < handle
properties
end
methods(Static)
end
methods
    function TEST_COHERE(obj)
        nF=obj.nF;
        f=obj.f;
        inds=1:obj.nStim;
        s=obj.S.s(:,inds);

        figure(1001)
    end
    function TEST_TRAIN(obj)
        Args={...
            'bTEST',true,'';
            'bPlot',false,'';
            ...
            'nRecFit',1,'';
            'bContinue',false,'';
            'nFset',1,'';
            'nWset',1,'';
            'FFix',[],'';
            'WFix',[],'';
            ...
            'bLrnF',true,'';
            'bRecurse',false,'';
            'bAll',false,'';
            'bAllOne',true,'';
            'inds',[],'';
            'rinds',[],'';
            'excl',[],'';
            ...
            'bLrnW',false,'';
            'bRecurseW',false,'';
            'bAllW',false,'';
            'bAllOneW',true,'';
            'indsW',[],'';
            'rindsW',[],'';
            'exclW',[],'';
        };
        Args=Args(:,1:2)';
        obj.train(Args{:});
    end
end
end
