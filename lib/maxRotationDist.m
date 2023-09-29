function [R,t,info] = maxRotationDist(R_avg,purse,kappa)
if nargin < 3
    kappa = 2;
end

A = purse.quadCon;
B = purse.linCon;

r_avg = R_avg(:);

d       = 12; 
x       = msspoly('x',d); % symbolic decision variables using SPOTLESS
r       = x(1:9);
t       = x(10:12);
R       = reshape(r,3,3);
c1      = R(:,1);
c2      = R(:,2);
c3      = R(:,3);
f       = - sum((r - r_avg).^2);
h       = [c1'*c1 - 1;
           c2'*c2 - 1;
           c3'*c3 - 1;
           c1'*c2;
           c2'*c3;
           c3'*c1;
           cross(c1,c2) - c3;
           cross(c2,c3) - c1;
           cross(c3,c1) - c2];
       
g       = [10^2 - sum(t.^2)]; % translation is bounded
for i = 1:size(A,3)
    Ai = squeeze(A(:,:,i));
    g  = [g; x' * (-Ai) * x];
end
slack = 1e-3;
for i = 1:size(B,2)
    Bi = squeeze(B(:,i));
    Bi = Bi(:);
    g  = [g; x' * Bi - slack];
end
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
problem.inequality      = g;
[SDP,relaxinfo]         = dense_sdp_relax(problem,kappa);

prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
Xmom = Xopt{1};
[V,~] = sorteig(Xmom);
v = V(:,1) / V(1,1);
R = project2SO3(reshape(v(2:10),3,3));
t = v(11:13);

if purse.checkMembership(R,t)
    f_sdp = -obj(1);
    f_est = norm(r_avg - R(:))^2;
    gap = abs(f_est - f_sdp) / (1 + abs(f_est) + abs(f_sdp));
    info.f_sdp = f_sdp;
    info.f_est = f_est;
    info.gap = gap;
    fprintf('f_sdp: %3.2e, f_est: %3.2e, gap: %3.2e.\n',f_sdp,f_est,gap);
else
    info.f_sdp = -obj(1);
    info.f_est = nan;
    info.gap = nan;
    fprintf('R and t not in PURSE, did not achieve global optimality.\n');
end
end