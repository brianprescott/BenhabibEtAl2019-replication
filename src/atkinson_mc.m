function mom    = atkinson_mc(p,param_fixed)

%% Paramters
% p is a vector for the 12 parameters to be estimated with shape 12 times 1
% first 5 are grid points on r; 6 is mu; 7 is A; 8 to 12 is diagonal
% element of transtional matrix of r

param.beta      = param_fixed.beta; % exogeneously given beta
param.gamma     = param_fixed.gamma;% exogeneously given gamma
param.nu        = exp(p(6));%mu
param.theta     = exp(p(7));%A
param.wage      = 1; 
param.Nr        = 5; % number of grid for r

options.pref        = 'crra';        % [log,crra]
options.spline      = 'cubic';       % [linear,cubic]
options.extrap      = 'no';          % ['yes','no'] (extrap in *solution*)
options.simulate    = 'yes';         % ['yes','no']
options.residplot   = 'no';
options.Cplot       = 'no';
options.nsim        = 100000; % number of people simulated
options.Tsim        = 500; % length of simulation period
 
%% State spaces

% Correlated shock to R
Nr          = param.Nr;
% sigz        = exp(p(1));
rgrid       = p(1:5);
Pr          = diag(p(8:12)); % diagonal element of 

fsolve_opt  = optimset('Display','off');
decayw      = fsolve(@(x) decay_weight(x),.6*ones(5,1),fsolve_opt);
pa1         = decayw(1);
pa2         = decayw(2);
pa3         = decayw(3);
pa4         = decayw(4);
pa5         = decayw(5);


Pr          = [p(8)             (1-p(8))*pa1    (1-p(8))*pa1^2  (1-p(8))*pa1^3  (1-p(8))*pa1^4;
               (1-p(9))*pa2     p(9)            (1-p(9))*pa2    (1-p(9))*pa2^2  (1-p(9))*pa2^3;
               (1-p(10))*pa3^2  (1-p(10))*pa3   p(10)           (1-p(10))*pa3   (1-p(10))*pa3^2;
               (1-p(11))*pa4^3  (1-p(11))*pa4^2 (1-p(11))*pa4   p(11)           (1-p(11))*pa4;
               (1-p(12))/4      (1-p(12))/4     (1-p(12))/4     (1-p(12))/4   p(12)];
            

for ii  = 1:Nr
    Pr(ii,end)   = 1-sum(Pr(ii,1:end-1));
end

%the above is just generating the transitional matrix for r

Prss        = Pr^10000; % steady state distribution of r
Prss        = Prss(1,:)'; % 5 times 1 matrix
Nr          = length(Prss); % 5
% rval 	      = rgrid'*Prss;
% rgrid       = ones(Nr,1)*rval;
Pr          = repmat(Pr,Nr,1); % no idea why this is done
betaER      = param.beta*rgrid'*Prss;
% betaER      = param.beta*rgrid;
fprintf('beta*E[R] = %1.3f\n',betaER);
% switch options.Pr
%     case 'iid'
%         Pr = repmat(Pr(1,:),Nr,1);
%     case 'persistent'
%         Pr = eye(Nr);
%     case 'corr'
%         Pr = Pr;
% end

% Epsilon shock space

% egrid       = [-.898 7.3 14.96 22.51 30.68 39.93 51.41 66.41 87.18 161.21]';
egrid       = linspace(1,10,10)';%10 times 1 vector
% elbar       = min(egrid);
% eubar       = max(egrid);
Ne          = length(egrid);%egrid is just index of quantile of wage

Pe          = [ 0.19635   0.15151   0.13048   0.11062   0.09344   0.07916   0.06791   0.06042   0.05457   0.05557;
                0.17597   0.14984   0.13101   0.11243   0.09772   0.08526   0.07398   0.06502   0.05723   0.05157;
                0.16198   0.14960   0.13101   0.11378   0.10002   0.08872   0.07808   0.06797   0.05943   0.04937;
                0.12051   0.12827   0.12390   0.11603   0.10837   0.09973   0.09169   0.08187   0.07209   0.05757;
                0.09469   0.10846   0.11303   0.11402   0.11130   0.10797   0.10231   0.09493   0.08523   0.06808;
                0.07579   0.08937   0.09937   0.10670   0.11068   0.11182   0.11187   0.10791   0.10119   0.08529;
                0.06061   0.07485   0.08658   0.09872   0.10803   0.11381   0.11701   0.11868   0.11572   0.10594;
                0.04869   0.06267   0.07571   0.09014   0.10348   0.11554   0.12385   0.12911   0.12900   0.12177;
                0.03794   0.05025   0.06331   0.07858   0.09461   0.11027   0.12591   0.13927   0.15065   0.14921;
                0.02748   0.03517   0.04560   0.05897   0.07235   0.08771   0.10740   0.13483   0.17487   0.25563];

for ii = 1:size(Pe,2)
    Pe(:,ii)     = Pe(:,ii)/sum(Pe(:,ii));
end

for ii = 1:size(Pe,1)
    Pe(ii,end)  = 1-sum(Pe(ii,1:end-1));
end
% the above is just transtional matrix for index of wage quantile

                

% Joint shock space / 
P           = kron(Pe,Pr); %Kronecker tensor product; transtional matrix of (50, 50)
% P           = Pe;
rgrid       = kron(ones(Ne,1),rgrid); % r value along the new state vector
egrid       = kron(egrid,ones(Nr,1)); % w index value along the new state vector
xgrid       = [rgrid,egrid];
K           = Ne*Nr;%lenth of new state vector

wgrid 	= [9.759892,19.95226,26.84919,33.05149,39.02348,45.05463,51.4006,59.1581,70.23087,100.3422;
	   11.54737,24.00785,32.58343,40.33186,47.70377,54.84686,63.09628,73.05923,87.21278,138.0883;
	   12.062,25.20177,34.95567,43.94973,52.42389,60.70431,69.4242,80.36736,97.51007,169.4898;
	   12.81191,26.4155,36.46041,45.55047,54.37142,63.08901,72.89075,85.09223,103.5408,182.405;
	   11.73645,24.66041,33.56299,42.22728,51.18447,60.34327,70.62766,82.78095,101.3924,183.3833;
	   8.222097,19.08026,26.77656,34.38757,42.95576,51.90672,61.64647,74.34956,93.41895,180.4287];
% Asset space
% Find Jess' amax
gridnumb    = 100000;%generate amax
ylbar       = param.wage*.001;%lower bound for asset; recall that value of wage is 1
param.ylbar = ylbar;

amin        = ylbar;
astep       = 0.1;
agrid       = ylbar+(0:gridnumb-1)*astep;%0.001+(100,000-1)*0.1

amax        = max(agrid);

Na          = 3000;%grid_number for asset
D           = 1;
curv        = 1;%curvature when generating agrid
agrid       = nodeunif(Na,amin.^curv,amax.^curv).^(1/curv);%generate a grid
%this function comes from compecon; is used to generate a

% Function space
switch options.spline
    case 'linear'
        spliorder = 1;
    case 'cubic'
        spliorder = 3;
end


%collocation function
fspace      = fundef({'spli',agrid,0,spliorder},...
                     {'spli',(1:1:K),0,1});%creat a function space
s           = gridmake(funnode(fspace));%Grid of evaluation nodes

                                       
ns          = length(s);%generate a new state vector incorporating asset; ns = 150100
agrid       = s((1:ns/(Ne*Nr)),1);%set the agrid to be the new grid

% Store the values
param.rgrid = rgrid;
param.egrid = egrid;
param.agrid = agrid;
param.wgrid = wgrid;
param.P     = P;
param.ns    = ns;
param.K     = K;

%% Solve problem

% Starting guess
v       = ones(ns,1);%value on the grid points of new state vector
cold    = funfitxy(fspace,s,[v,v,v,v,v,v]);%coefficients of collocation

% Bellman iterations: collocation method
tolc    = 1e-20;
for iterc   = 1:100;
    [C,v1,v2,v3,v4,v5,v6,C1,C2,C3,C4,C5]    = solve(cold,fspace,s,param,options);
    c       = funfitxy(fspace,s,[v1,v2,v3,v4,v5,v6]);
    dc      = norm(c-cold)/norm(cold);
    cold    = c;
    if dc<tolc,break,end
    fprintf('%i\tdc:%1.3e\n',iterc,dc);
end

%% Fine grid
fprintf('Solving on fine grid\n');
Nafine      = Na*D;%you kidding me, D = 1
agridfine   = nodeunif(Nafine,amin.^curv,amax.^curv).^(1/curv);
fspacefine  = fundef({'spli',agridfine,0,spliorder},{'spli',(1:1:K)',0,1});
sfine       = gridmake(funnode(fspacefine));
agridfine   = sfine(1:fspacefine.n(1),1);
ns          = length(sfine);
param.ns    = ns;

%% Check residuals
if strcmp(options.residplot,'yes');
    [~,v1,v2,v3,v4,v5,v6]     = solve(c,fspace,sfine,param,options);
    resid1      = funbas(fspace,sfine)*c(:,1)-v1;
    resid2      = funbas(fspace,sfine)*c(:,2)-v2;
    resid3      = funbas(fspace,sfine)*c(:,3)-v3;
    resid4      = funbas(fspace,sfine)*c(:,4)-v4;
    resid5      = funbas(fspace,sfine)*c(:,5)-v5;
    resid6      = funbas(fspace,sfine)*c(:,6)-v6;

    figure(1);
    subplot(6,1,1);
    plot(resid1); title('A. resid1: \Phic^1-v_1(c)');grid on;
    xlim([0,length(resid1)]);
    subplot(6,1,2);
    plot(resid2); title('B. resid2; \Phic^2-v_2(c)');grid on;
    xlim([0,length(resid2)]);
    subplot(6,1,3);
    plot(resid3); title('B. resid3; \Phic^3-v_3(c)');grid on;
    xlim([0,length(resid3)]);
    subplot(6,1,4);
    plot(resid4); title('B. resid4; \Phic^4-v_4(c)');grid on;
    xlim([0,length(resid4)]);
    subplot(6,1,5);
    plot(resid5); title('B. resid5; \Phic^5-v_5(c)');grid on;
    xlim([0,length(resid5)]);
    subplot(6,1,6);
    plot(resid6); title('B. resid6; \Phic^6-v_6(c)');grid on;
    xlim([0,length(resid6)]);
end

%% Plot consumption function
if strcmp(options.Cplot,'yes');
    splot       = gridmake(agridfine,1);
    param.ns    = size(splot,1);
    C           = solve(c,fspace,splot,param,options);
    
    figure(2);
%     plot(s(1:length(agrid),1),C(1:length(agrid))./s(1:length(agrid),1));grid on;
%     title('c(a)/a');
%     plot(s(1:length(agrid),1),C(1:length(agrid)));grid on;
%     fprintf('At amax, c(a)/a = %1.3f\n',C(length(agrid))/agrid(end));
    
    plot(agridfine,C./agridfine);grid on;
    title('c(a)/a');
    plot(agridfine,C);grid on;
    fprintf('At amax, c(a)/a = %1.3f\n',C(end)/agridfine(end));
end

%% Stationary distribution
fprintf('Solving stationary distribution\n');

param.ns    = ns;
tic
% C           = solve(c,fspace,sfine,param,options);
% x           = s(:,1)-C;
fspace_erg  = fundef({'spli',agridfine,0,1});
Fprimecell  = cell(K,1);
Pprimecell  = cell(K,1);
amaxvec     = zeros(ns,1) + max(s(:,1));
aminvec     = zeros(ns,1) + min(s(:,1));

Q           = sparse(ns,ns);
% w           = param.wage;
for kk = 1:K
    R               = param.rgrid(kk*ones(ns,1));
    E               = param.egrid(kk*ones(ns,1));
    abar            = zeros(ns,1)+max(s(:,1));
    Aprime_kk       = max(min(C,amaxvec),aminvec);
    Fprime          = funbas(fspace_erg,Aprime_kk);%collocation matrix
    Pprime          = zeros(ns,K);
    Pprime(:,kk)    = P(sfine(:,2),kk);     % Transition only to kk
    % Combine transition matrices
    Q           = Q + dprod(Pprime,Fprime);
end
assert(abs(full(sum(sum(Q)))-ns)<1e-6,'Error: Q');

L           = ones(ns,1)./ns;%initial distribution
itermaxL    = 5000;%max iteration
tolL        = 1e-10;%tolerance
for t = 1:itermaxL
    Lnew        = Q'*L;
    dL          = norm(Lnew-L)/norm(L);
    L           = Lnew;
    if dL<tolL,break,end
    if (mod(t,500)==0)
        fprintf('%i\tdc:%1.3e\n',t,dL);
    end
end
%above finds the stationary distribution

La = kron(ones(1,K),speye(ns/K))*L;
toc


%% Simulation

% Simulation parameters and storage
Ni      = options.nsim;
Nt      = options.Tsim;
dist0   = 'zero'; %['zero'|'SS']

% Fit linear interpolant for C (so can extrapolate)
% [C,~]   = solve(c,fspace,sfine,param,options);
fspaceC = fundef({'spli',agridfine,0,1},{'spli',(1:1:K),0,1});
sC      = gridmake(funnode(fspaceC));
cC      = funfitxy(fspaceC,sC,C);
cC1     = funfitxy(fspaceC,sC,C1);%policy function collocation coefficient
cC2     = funfitxy(fspaceC,sC,C2);
cC3     = funfitxy(fspaceC,sC,C3);
cC4     = funfitxy(fspaceC,sC,C4);
cC5     = funfitxy(fspaceC,sC,C5);

% Initialise simulation from stationary distribution without extrapolation
% As their is an erodicity theorem, this is just a starting guess
switch dist0
    case 'SS'
        Lround      = round(Ni*L);
        si0         = [];
        for i = 1:ns;
            si0 	= [si0;repmat(sfine(i,:),Lround(i),1)];
        end
        Ni          = size(si0,1);
    case 'zero'
        si0         = repmat([ylbar,1],Ni,1);
end
%above just generate a matrix that specifies the initial state of the
%simulated people

asimu       = zeros(Ni,Nt);
a1simu      = zeros(Ni,Nt);
a2simu      = zeros(Ni,Nt);
a3simu      = zeros(Ni,Nt);
a4simu      = zeros(Ni,Nt);
a5simu      = zeros(Ni,Nt);%asset holding at different stage in life
ksimu       = zeros(Ni,Nt);
asimu(:,1)  = si0(:,1);
a1simu(:,1) = si0(:,1);
a2simu(:,1) = si0(:,1);
a3simu(:,1) = si0(:,1);
a4simu(:,1) = si0(:,1);
a5simu(:,1) = si0(:,1);
ksimu(:,1)  = si0(:,2);%state in r and w

tic
rng(111);
Pcumsum     = cumsum(P,2); % cdf
Pcumsum(:,end) = 1;
PS          = rand(Ni,Nt);   % Shock
parfor i = 1:Ni;
%     i
    si = si0(i,:);
    for t = 2:Nt;
%         A           = si(1,1);
        Aprime      = max(funeval(cC,fspaceC,si),0);
        A1          = max(funeval(cC1,fspaceC,si),0);
        A2          = max(funeval(cC2,fspaceC,si),0);
        A3          = max(funeval(cC3,fspaceC,si),0);
        A4          = max(funeval(cC4,fspaceC,si),0);
        A5          = max(funeval(cC5,fspaceC,si),0);%these are next period asset
%         C           = funeval(cC,fspaceC,si);
%         X           = A-C;
        shock       = PS(i,t);
        kprime      = find(Pcumsum(si(1,2),:)>shock,1,'first');
%         Aprime      = C;
        asimu(i,t)  = Aprime;
        a1simu(i,t) = A1;
        a2simu(i,t) = A2;
        a3simu(i,t) = A3;
        a4simu(i,t) = A4;
        a5simu(i,t) = A5;
        ksimu(i,t)  = kprime;
        si          = [Aprime,kprime];
    end
end
toc
% saving('ALL_sim_IID.mat','-v7.3');

%% QUINTILES FROM SIM


%last 50 years
aSSfinal    = vec(asimu(:,end-50:end));
aSS1        = vec(a1simu(:,end-50:end));
aSS2        = vec(a2simu(:,end-50:end));
aSS3        = vec(a3simu(:,end-50:end));
aSS4        = vec(a4simu(:,end-50:end));
aSS5        = vec(a5simu(:,end-50:end));
aSS         = [aSS1;aSS2;aSS3;aSS4;aSS5;aSSfinal];%50*Ni total population in the economy

percvec     = [.2 .4 .6 .8 .9 .95 .99]';
Nperc       = length(percvec);
Q           = quantile(aSS,percvec);%compute quantile
aTOTAL      = sum(aSS);
aSS         = sort(aSS);
aJ          = zeros(Nperc+2,1);
aJ(1)       = 1;
aJ(end)     = length(aSS);%document the number of pp above the quantile
aqTOTAL     = zeros(Nperc+1,1);%document pp with in each quantile

for ii  = 1:Nperc
    aJ(ii+1)        = find(aSS >= Q(ii),1,'first');
    aqTOTAL(ii)     = sum(aSS(aJ(ii):aJ(ii+1)-1));
end
aqTOTAL(end)    = sum(aSS(aJ(end-1):end));

aSHARE      = aqTOTAL/aTOTAL;%only seven due to sum euqalling one

%% Transition matrix(mobility matrix)
%refer to Charles and Hurst (2003)
tgap    = 50;
Pcell   = zeros(tgap,5,5);
for gg = 1:tgap

%     a0      = asimu(:,end-tgap);
%     a1      = asimu(:,end-tgap+1);

%     [~,~,child_resid,~]     = reg(a4simu(:,end-tgap+1),wgrid(4,egrid(ksimu(:,end-tgap+1)))',0.05,1);
%     [~,~,parent_resid,~]    = reg(a4simu(:,end-tgap),wgrid(4,egrid(ksimu(:,end-tgap)))',0.05,1);
 
%     [~,~,child_resid,~]     = reg([a1simu(:,end-tgap+1);a2simu(:,end-tgap+1);a3simu(:,end-tgap+1); ...
%                                    a4simu(:,end-tgap+1);a5simu(:,end-tgap+1);asimu(:,end-tgap+1)], ...
%                                    [ones(size(asimu(:,end-tgap+1))),wgrid(1,egrid(ksimu(:,end-tgap+1)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*2,wgrid(2,egrid(ksimu(:,end-tgap+1)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*3,wgrid(3,egrid(ksimu(:,end-tgap+1)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*4,wgrid(4,egrid(ksimu(:,end-tgap+1)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*5,wgrid(5,egrid(ksimu(:,end-tgap+1)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*6,wgrid(6,egrid(ksimu(:,end-tgap+1)))'],0.05,1);
                               
%     [~,~,parent_resid,~]     = reg([a1simu(:,end-tgap);a2simu(:,end-tgap);a3simu(:,end-tgap); ...
%                                    a4simu(:,end-tgap);a5simu(:,end-tgap);asimu(:,end-tgap)], ...
%                                    [ones(size(asimu(:,end-tgap+1))),wgrid(1,egrid(ksimu(:,end-tgap)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*2,wgrid(2,egrid(ksimu(:,end-tgap)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*3,wgrid(3,egrid(ksimu(:,end-tgap)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*4,wgrid(4,egrid(ksimu(:,end-tgap)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*5,wgrid(5,egrid(ksimu(:,end-tgap)))'; ...
%                                    ones(size(asimu(:,end-tgap+1)))*6,wgrid(6,egrid(ksimu(:,end-tgap)))'],0.05,1);
                               
    [~,~,child_resid,~]     = reg([a1simu(:,end-tgap+1);a2simu(:,end-tgap+1);a3simu(:,end-tgap+1); ...
                                   a4simu(:,end-tgap+1);a5simu(:,end-tgap+1);asimu(:,end-tgap+1)], ...
                                   [ones(size(asimu(:,end-tgap+1))); ...
                                   ones(size(asimu(:,end-tgap+1)))*2; ...
                                   ones(size(asimu(:,end-tgap+1)))*3; ...
                                   ones(size(asimu(:,end-tgap+1)))*4; ...
                                   ones(size(asimu(:,end-tgap+1)))*5; ...
                                   ones(size(asimu(:,end-tgap+1)))*6],0.05,1);
                               
    [~,~,parent_resid,~]     = reg([a1simu(:,end-tgap);a2simu(:,end-tgap);a3simu(:,end-tgap); ...
                                   a4simu(:,end-tgap);a5simu(:,end-tgap);asimu(:,end-tgap)], ...
                                   [ones(size(asimu(:,end-tgap+1))); ...
                                   ones(size(asimu(:,end-tgap+1)))*2; ...
                                   ones(size(asimu(:,end-tgap+1)))*3; ...
                                   ones(size(asimu(:,end-tgap+1)))*4; ...
                                   ones(size(asimu(:,end-tgap+1)))*5; ...
                                   ones(size(asimu(:,end-tgap+1)))*6],0.05,1);
    
    a0      = parent_resid;
    a1      = child_resid;

%     quin0   = prctile(a0,[20,40,60,80,90,95,99]);
%     quin1   = prctile(a1,[20,40,60,80,90,95,99]);

    quin0   = prctile(a0,[20,40,60,80]);
    quin1   = prctile(a1,[20,40,60,80]);
    
    id0     = 0*a0;
    id1     = 0*a1;

    id0(              a0<=quin0(1)) = 1;
    id0(a0>quin0(1) & a0<=quin0(2)) = 2;
    id0(a0>quin0(2) & a0<=quin0(3)) = 3;
    id0(a0>quin0(3) & a0<=quin0(4)) = 4;
%     id0(a0>quin0(4) & a0<=quin0(5)) = 5;
%     id0(a0>quin0(5) & a0<=quin0(6)) = 6;
%     id0(a0>quin0(6) & a0<=quin0(7)) = 7;
%     id0(a0>quin0(7))                = 8;
    id0(a0>quin0(4))                = 5;

    id1(              a1<=quin1(1)) = 1;
    id1(a1>quin1(1) & a1<=quin1(2)) = 2;
    id1(a1>quin1(2) & a1<=quin1(3)) = 3;
    id1(a1>quin1(3) & a1<=quin1(4)) = 4;
%     id1(a1>quin1(4) & a1<=quin1(5)) = 5;
%     id1(a1>quin1(5) & a1<=quin1(6)) = 6;
%     id1(a1>quin1(6) & a1<=quin1(7)) = 7;
%     id1(a1>quin1(7))                = 8;
    id1(a1>quin1(4))                = 5;

%     Pest = zeros(8,8);
    Pest = zeros(5,5);
    for ii = 1:5;
        for jj = 1:5;
            Pest(ii,jj) = sum((id0==ii)&(id1==jj))/sum(id0==ii);
        end
    end
    Pcell(gg,:,:)   = bsxfun(@rdivide,Pest,sum(Pest,2));
end

Pest    = squeeze(sum(Pcell,1))/tgap;
aMOBIL  = diag(Pest);
disp(Pest)

mom     = [aSHARE;aMOBIL];%return moments generated

end
