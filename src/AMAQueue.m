% trainData '_' normType '_' fitType
stAlias='DSP21_1';
normType='';
amaType='GSS';
fitType='';
bMeanCon=0;

otherArgs=struct('nF',8,'nFset',1);

lstArgs={'normType'};
lst={...
    %'cov2';
    %'snew';
    %'cov3';

    %'ratio';
    %'sign';  % BAD
    %'sigma'; % BAD
    %'w';
    %'nw';
    %'w02';
    %'nw02';

    %'w2';
    %'nw2';
    %'w22';
    %'nw22';


    %'new1';
    %'w2act3';
    %'w2act4';
    %'new2';
    %'new3';
    %'new4';
    %%%'new70'; % BAD
    %'new80';
    %'new9';

    %'nnew1';
    %'new5';
    %'new7';
    %'nnew2';
    %'nnew3';
    %'nnew4';

    %'nnew7';
    %'nnew8';
    %'nnew9';


    'W1';

    %'mnew1';
    %'mnew5';
    %'mnew6';
    %'mnew2';
    %'mnew3';
    %'mnew4';

    %'mnew7';
    %'mnew8';
    %'mnew9';

    %'mnnew1';
    %'mnnew5';
    %'mnnew6';
    %'mnnew2';
    %'mnnew3';
    %'mnnew4';

    %%'mnnew7'; BAD
    %'mnnew8';
    %'mnnew9';

    %'onew1';
    %'onew5';
    %'onew6';
    % HERE
    %%'onew2'; STUCK
    %'onew3';
    %'onew4';

    %'onew7';
    %'onew8';
    %'onew9';


    %'onnew1';
    %'new6';
    %%'new8'; BAD
    %'onnew2';
    %'onnew3';
    %'onnew4';

    %'onnew7';
    %'onnew8';
    %'onnew9';

    %'snew1';
    %'snew2';
    %'snew3';
    %'snew4';
    %'snew5';
    %'snew6';
    %'snew7';
    %'snew8';
    %'snew9';

    %'snew3';
    %'snew31';
    %'snew32';
    %'snew33';
    %'snew3';
    %'snew34';
    %'snew35';
    %'snew36';
    %'snew37';
    %'snew38';
    %'snew39';
    %'snew310';
    %'snew311';
    %'snew312';
    %%'snew313';
    %%'snew314';

    %'dnew1';
    %'dnew2';
    %'dnew3';
    %'dnew4';
    %'dnew5';
    %'dnew6';

    %W1
    %W2
};


%AMAQ(baseArgs);


% new6 :: 0.75

% base            2.53
% hinge           2.53
%
%                                 n   nn an  a    d  s  F  G
% sum     -      : 2.62          RS                        -
% sum2    -      : 1.90         R2S                        -
% sqr     -      : 1.90           r   2                    -
%
% WEIGHTS
% w       -      : 1.81      w   WR                        -
% nw      -      : 1.45     -w   WR                        -
% w02     -      : 1.76     w2   WR                        -
% nw02    -      : 1.45    -w2   WR                        -
%
% w2      -      : 1.75      w  WR2                        -
% nw2     -      : 1.56     -w  WR2                        -
% w22     -      : 1.82     w2  WR2                        -
% nw22    -      : 1.56    -w2  WR2                        -
%
% ratio   -      : B              r              RS        -
% sign    -      :                r   2S         RS        -
% act     -      : 2.18           r             R2S        -
% act2    -      : 2.12           r   2         R2S        -
% sign2   -      : 2.19           r   2S        R2S        -
% sigma   -      :                m   2S        R2S        -
%
% new1    - 1.92 : 1.63      w    WR           WR2S        -
% w2act3  - 1.76 : 1.19      w    WR  2        WR2S        -
% w2act4  -      : 1.26      w    WR     2     WR2S        -
% new2    - 2.75 : 2.40      w     r           WR2S        -
% new3    - 2.58 : 2.19      w     r  2        WR2S        -
% new4    - 2.60 : 1.91      w     r     2     WR2S        -
% new70   -      : B         w     r     2     WR2S        -
% new80   -      : 1.59      w    WR        2   WRS        -
% new9    -      : B             WR2        2   WRS        -
%
% nnew1   - 2.07 : 1.58     -w    WR           WR2S        -
% new5    -      : 1.13     -w    WR  2        WR2S        -
% new7    -      : 1.58     -w    WR     2     WR2S        -
% nnew2   - 2.60 : 2.41     -w     r           WR2S        -
% nnew3   - 2.57 : 2.47     -w     r  2        WR2S        -
% nnew4   - 2.50 : 2.17     -w     r     2     WR2S        -
% nnew7   -      :          -w   WR2     2     WR2S        -
% nnew8   -      : B        -w    WR        2   WRS        -
% nnew9   -      :          -w   WR2        2   WRS        -

% -1 to 1 weights                                          -
% mnew1   -      : 2.05     w2    WR        A  WR2S        -
% mnew5   -      : 1.46     w2    WR  2     A  WR2S        -
% mnew6   -      : 1.89@4   w2    WR     2  A  WR2S        -
% mnew2   -      : ...      w2     r        A  WR2S        -
% mnew3   -      : 2.17@4   w2     r  2     A  WR2S        -
% mnew4   -      : 1.86@5   w2     r     2  A  WR2S        -
% mnew7   -      : s        w2   WR2        A  WR2S        -
% mnew8   -      :          w2    WR        2   WRS        -
% mnew9   -      :          w2   WR2        2   WRS        -
%
% mnnew1  -      : 1.58    -w2    WR        A  WR2S        -
% mnnew5  -      : 1.13    -w2    WR  2     A  WR2S        -
% mnnew6  -      : 1.58    -w2    WR     2  A  WR2S        -
% mnnew2  -      : 2.41    -w2     r        A  WR2S        -
% mnnew3  -      : 2.47    -w2     r  2     A  WR2S        -
% mnnew4  -      : 2.17    -w2     r     2  A  WR2S        -
% mnnew7  -      :         -w2   WR2        A  WR2S        -
% mnnew8  -      :         -w2    WR        2   WRS        -
% mnnew9  -      : x       -w2   WR2        2   WRS        -
%
% no abs denom
% onew1   -      : 2.3@4    w2    WR           WR2S        -
% onew5   -      : 1.52@6   w2    WR  2        WR2S        -
% onew6   -      : 1.89@4   w2    WR     2     WR2S        -
% onew2   -      : ...      w2     r           WR2S        -
% onew3   -      : nan @?   w2     r  2     a  WR2S        -
% onew4   -      : 1.86@5   w2     r     2     WR2S        -
% onew7   -      :          w2   WR2           WR2S        -
% onew8   -      : B        w2    WR        2   WRS        -
% onew9   -      :          w2   WR2        2   WRS        -
%                                                          -
% onnew1  -      : 1.65    -w2    WR           WR2S        -
% new6    -      : 1.13    -w2    WR  2        WR2S        -
% new8    -      : 1.58    -w2    WR     2     WR2S        -
% onnew2  -      : ...     -w2     r           WR2S        -
% onnew3  -      : 2.47    -w2     r  2     a  WR2S        -
% onnew4  -      : 2.17    -w2     r     2     WR2S        -
% onnew7  -      :         -w2   WR2           WR2S        -
% onnew8  -      :         -w2    WR        2   WRS        -
% onnew9  -      :         -w2   WR2        2   WRS        -
%
% strattle
% snew1   -  2.60 : 1.44         -w2    WR           WR2S  1     -
% snew2   -  2.59 : 1.09         -w2    WR  2        WR2S  1     -
% snew3   -  2.58 : 1.07         -w2    WR     2     WR2S  1     -
% snew4   -  2.71 : 2.42         -w2     r           WR2S  1     -
% snew5   -  2.59 : 2.45         -w2     r  2     a  WR2S  1     -
% snew6   -  2.63 : 2.17         -w2     r     2     WR2S  1     -
% snew7   -  2.58 : 1.77         -w2   WR2           WR2S  1     -
% snew8   -  3.23 : 2.98         -w2    WR        2   WRS  1     -
% snew9   -  3.27 : 3.64         -w2   WR2        2   WRS  1     -
% snew    -      :          -w2   WR2        2   WRS -1     -
% snew31  -      :          -w2   WR2        2   WRS -1     -
% snew32  -      :          -w2   WR2        2   WRS -1     -
% XXX check for abs in numerator
%-----------------------
% new5  n^2
% new6  n^2
% w2act3  n^2
% w2act4 b^2
%
% new1
% mnew1
% nnew1
%
%---
% cov 22, 6 - 0.38 @ 783
% cov 22, 3 - 1.03 @ 783
%
% Full @ 100
% snew3
% snew31 1.42
% snew32 1.76
% snew33 1.42
% snew34 0.62 (0.51)
% snew35 looks same as 34
% snew36 1.70
% snew37 1.47
% snew38 1.90
% snew39 0.73
% snew310 nan
% snew312 0.71 (0.63)
% snew313 0.59
% snew314 0.59
%
% snew313 0.59
% snew314 0.59
% snew34 0.62 (0.51)
%
%- 1
% Full @ 100
% dnew1   1.86
% dnew2   1.42
% dnew3   1.39
% dnew4   2.50
% dnew5   1.32
% snew313 0.31
% snew314
