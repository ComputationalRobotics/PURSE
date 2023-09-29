clc
clear
close all
addpath('lib')
addpath(genpath('p3p_problem'))
addpath(genpath('spotless'))

% Add mosekpath according to its installation path
% restoredefaultpath
% mosekpath = "../mosek";
% addpath(genpath(mosekpath))

rng(0)

%% generate a randome problem
N       = 8;
problem = randPoseProblem(N);

%% create purse
purse = PURSE(problem.centers, problem.radii, problem.kpts3d);

%% check membership of groundtruth, should always be true
inpurse = purse.checkMembership(problem.R,problem.t);
fprintf("gt pose is inside purse: %s \n",mat2str(inpurse))

%% random samples and average pose
[Rset,tset,Ravg,tavg] = purse.ransag(100000);
R_err = getAngularError(Ravg,problem.R);
t_err = getTranslationError(tavg,problem.t);
fprintf("Num of samples: %d.\n",length(Rset));

%% compute error bound: rotation
[~,~,Rinfo] = maxRotationDist(Ravg,purse,2);
R_bound = rad2deg(acos(1 - Rinfo.f_sdp/4));

%% compute error bound: translation
[~,~,tinfo] = maxTranslationDist(tavg,purse,2);
t_bound = sqrt(tinfo.f_sdp);

fprintf("rotation err: %2.4f (deg), translation err: %2.4f (m).\n",R_err,t_err)
fprintf("rotation bound: %2.4f (deg), translation bound: %2.4f (m).\n",R_bound,t_bound)


