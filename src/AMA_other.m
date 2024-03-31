classdef AMA_other < handle
properties
    other=struct('ica',[],'pca',[],'sparse',[],'factoran',[],'svm',[],'ecoc',[]);
end
methods(Static)
end
methods
    function svm(obj,bFourier)
        if nargin < 3 ||isempty(bFourier)
            bFourier=obj.Stim.bFourier;
        end
        if bFourier
            S=obj.Stim.SF;
        else
            S=obj.Stim.SS;
        end

        %C=obj.Stim.ctgInd;
        C=obj.Stim.XCtg;


        K='polynomial';
        %K='rbf';
        %mdl=fitrsvm(S',C,'KernelFunction',K,'CategoricalPredictors','all');
        mdl=fitrgp(S',C,'OptimizeHyperparameters','auto');

        obj.other.svm.mdl=mdl;

        %f=obj.other.svm.mdl.SupportVectors';
        f=obj.other.svm.mdl.ActiveSetVectors';
        nF=size(f,2);
        %figure(1)
        %obj.plotSVMR();
        n=min([30,nF]);
        obj.plotSVM(n);
    end
    function ecoc(obj,bFourier)
        if nargin < 3 ||isempty(bFourier)
            bFourier=obj.Stim.bFourier;
        end
        if bFourier
            S=obj.Stim.SF;
        else
            S=obj.Stim.SS;
        end

        C=obj.Stim.ctgInd;
        C=obj.Stim.XCtg;
        code='sparserandom';
        %code='denserandom';
        code='onevsall';
        t=templateLinear('Learner','svm','Regularization','ridge','Solver','dual','PassLimit',3,'ObservationsIn','columns');
        mdl = fitcecoc(S',C,'Coding',code,'Learners',t);

        obj.other.ecoc.mdl=mdl;
        nF=numel(obj.other.ecoc.mdl.BinaryLearners);
        f=zeros(obj.Stim.nPix,nF);
        for i = 1:nF
            f(:,i)=obj.other.ecoc.mdl.BinaryLearners{i}.Beta;
        end
        obj.other.ecoc.f=f;
        obj.other.ecoc.bFourier=bFourier;

        %figure(1)
        %obj.plotSVMR();
        %figure(2)
        %n=min([30,nF]);
        %obj.plotSVM(n);
    end
    function pca(obj,bFourier)
        if nargin < 3 ||isempty(bFourier)
            bFourier=obj.Stim.bFourier;
        end
        if bFourier
            S=obj.Stim.SF;
        else
            S=obj.Stim.SS;
        end

        lambda = pca(S');

        obj.other.pca.lambda=lambda;
        obj.other.pca.bFourier=bFourier;
    end
    function factoran(obj,n,bFourier)
        if nargin < 3 ||isempty(bFourier)
            bFourier=obj.Stim.bFourier;
        end
        if bFourier
            S=obj.Stim.SF;
        else
            S=obj.Stim.SS;
        end

        lambda = factoran(S',n);

        obj.other.factoran.lambda=lambda;
        obj.other.factoran.bFourier=bFourier;
    end
    function ICA(obj,n,bFourier)
        if nargin < 3 ||isempty(bFourier)
            bFourier=obj.Stim.bFourier;
        end
        if bFourier
            S=obj.Stim.SF;
        else
            S=obj.Stim.SS;
        end

        Mdl = rica(S',n,'IterationLimit',100);

        obj.other.ica.mdl=Mdl;
        obj.other.ica.bFourier=bFourier;
    end
    function sparse(obj,n,bFourier)
        if nargin < 3 ||isempty(bFourier)
            bFourier=obj.Stim.bFourier;
        end
        if bFourier
            S=obj.Stim.SF;
        else
            S=obj.Stim.SS;
        end

        Mdl = sparsefilt(S,n,'IterationLimit',100);

        obj.other.sparse.mdl=Mdl;
        obj.other.sparse.bFourier=bFourier;
    end
    function plotSparse(obj)
        f=obj.other.sparse.mdl.TransformWeights;
        n=size(f,2);
        obj.Plot.F_(n,false,false,[],f);
    end
    function plotICA(obj,n)
        f=obj.other.ica.mdl.TransformWeights;
        if nargin < 2 || isempty(n)
            n=size(f,2);
        end
        obj.Plot.F_(n,false,false,[],f);
    end
    function plotFactoran(obj)
        f=obj.other.factoran.lambda;
        obj.Plot.F_(n,false,false,[],f);
    end
    function plotPCA(obj)
        f=obj.other.pca.lambda;
        n=size(f,2);
        obj.Plot.F_(n,false,false,[],f);
    end
    function plotSVM(obj,n)
        %f=obj.other.svm.mdl.SupportVectors';
        f=obj.other.svm.mdl.ActiveSetVectors';
        if nargin < 2 || isempty(n)
            n=size(f,2);
        end
        obj.Plot.F_(n,false,false,[],f);
    end
    function plotSVMR(obj,fPairs,ctg2plt)
        if nargin < 2
            fPairs=[];
        end
        if nargin < 3
            ctg2plt=[];
        end
        P=obj.Plot.initRParms();
        %f=obj.other.svm.mdl.SupportVectors';
        f=obj.other.svm.mdl.ActiveSetVectors';
        obj.Plot.R_(fPairs,ctg2plt,[],[],P,f);
    end
    function plotECOC(obj,n)
        f=obj.other.ecoc.f;
        if nargin < 2 || isempty(n)
            n=size(f,2);
        end
        obj.Plot.F_(n,false,false,[],f);
    end
    function plotECOCR(obj,fPairs,ctg2plt)
        if nargin < 2
            fPairs=[];
        end
        if nargin < 3
            ctg2plt=[];
        end
        P=obj.Plot.initRParms();
        f=obj.other.ecoc.f;
        obj.Plot.R_(fPairs,ctg2plt,[],[],P,f);
    end
end
end
