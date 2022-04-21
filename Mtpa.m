clear; close all; clc;

load FluxMap.mat 

PsiD = @(param,id,iq) param(1) + param(2)*id + param(3)*iq + param(4)*id.^2 + param(5)*id.*iq + param(6)*iq.^2;
PsiQ = @(param,id,iq) (param(1)*id.^2+param(2)*id+param(3)).*(param(4)*tanh(param(5)*iq)+param(6)*iq)+param(7);


dPsiDdid = @(param,id,iq) param(2) + 2*param(4)*id + param(5)*iq;
dPsiDdiq = @(param,id,iq) param(3) + param(5)*id + 2*param(6)*iq;

dPsiQdiq = @(param,id,iq) (param(1)*id.^2+param(2)*id+param(3)).*(param(4)*param(5)*sech(param(5)*iq).^2+param(6));
dPsiQdid = @(param,id,iq) (2*param(1)*id+param(2)).*(param(4)*tanh(param(5)*iq)+param(6)*iq)+param(7);

dPsiD2did2 = @(param,id,iq) 2*param(4);
dPsiD2diq2 = @(param,id,iq) 2*param(6);
dPsiD2didiq = @(param,id,iq) param(5);

dPsiQ2did2 = @(param,id,iq) 2*param(1).*(param(4)*tanh(param(5)*iq)+param(6)*iq)+param(7);
dPsiQ2diq2 = @(param,id,iq) (param(1)*id.^2+param(2)*id+param(3)).*(-2*param(4)*param(5)^2*sech(param(5)*iq).^2.*tanh(param(5)*iq));
dPsiQ2diqid = @(param,id,iq) (2*param(1)*id+param(2)).*(param(4)*param(5)*sech(param(5)*iq).^2+param(6));

figure('Name','Flux Map psid');
surfc(id,iq,psid);
xlabel('id [A]');
ylabel('iq [A]');
zlabel('\psi_d [A]');

figure('Name','Flux Map psiq');
surfc(id,iq,psiq);
xlabel('i_d [A]');
ylabel('i_q [A]');
zlabel('\psi_q [A]');



ERR_D = @(PD) psid(:)-PsiD(PD,id(:),iq(:));
ERR_Q = @(PQ) psiq(:)-PsiQ(PQ,id(:),iq(:));

[PD_opt,RES_fitD,RES_D]  = lsqnonlin(ERR_D,zeros(6,1));
[PQ_opt,RES_fitQ,RES_Q]  = lsqnonlin(ERR_Q,[0;1;0;1;1;0;0],zeros(7,1),Inf*ones(7,1));


figure('Name','Flux Map psid - MODEL');
surfc(id,iq,PsiD(PD_opt,id,iq));
xlabel('id [A]');
ylabel('iq [A]');
zlabel('\psi_d [A]');

figure('Name','Flux Map psiq - MODEL');
surfc(id,iq,PsiQ(PQ_opt,id,iq));
xlabel('i_d [A]');
ylabel('i_q [A]');
zlabel('\psi_q [A]');

% MTPA
% Minimize the current amplitude id² + iq² under the constraint to follow
% the reference torque so that Temg = 3/2 * p * (psid*iq-psiq*id)

p = 8;  % Number of pole pair

Temg = @(id,iq) 3/2*p*(PsiD(PD_opt,id,iq)*iq-PsiQ(PQ_opt,id,iq)*id);
dTemgdid = @(id,iq) 3/2*p*(dPsiDdid(PD_opt,id,iq)*iq-dPsiQdid(PQ_opt,id,iq)*id-PsiQ(PQ_opt,id,iq));
dTemgdiq = @(id,iq) 3/2*p*(dPsiDdiq(PD_opt,id,iq)*iq+PsiD(PD_opt,id,iq)-dPsiQdiq(PQ_opt,id,iq)*id);

dTemg2did2 = @(id,iq) 3/2*p*(dPsiD2did2(PD_opt,id,iq)*iq-dPsiQ2did2(PQ_opt,id,iq)*id-2*dPsiQdid(PQ_opt,id,iq));
dTemg2didiq = @(id,iq) 3/2*p*(dPsiD2didiq(PD_opt,id,iq)*(iq-id)+dPsiDdid(PD_opt,id,iq)-dPsiQdiq(PQ_opt,id,iq));
dTemg2diq2 = @(id,iq) 3/2*p*(dPsiD2diq2(PD_opt,id,iq)*iq+2*dPsiDdiq(PD_opt,id,iq)-dPsiQ2diq2(PQ_opt,id,iq)*id);

% 1 - Determination of the cost function = > Lagrangian
% Objective is minimize current id^2+iq^2
% Under the constraint Temg(id,iq)-Tref) = 0
% The constraint is taken into account with the Lagrange multiplier l1
% which becomes a new auxiliary variable
% It leads to the lagrangian 
% L = id^2+iq^2+l1*(Temg(id,iq)-Tref) to minimize
TREF = 0:5:80;  Ntorque = length(TREF);
ID = 0*TREF;    IQ = 0*TREF;
x0 = [0;0;10];

for ntorque = 1 : Ntorque
    Tref = TREF(ntorque);  % reference torque
     
    L = @(id,iq,l1) id^2+iq^2+l1*(Temg(id,iq)-Tref);
    
    % Minimizing the function involves cancelling its derivatives 
    dLdid = @(id,iq,l1) 2*id+l1*dTemgdid(id,iq);
    dLdiq = @(id,iq,l1) 2*iq+l1*dTemgdiq(id,iq);
    dLdl1 = @(id,iq,l1) Temg(id,iq)-Tref;
    
    % We can employ a Newton with the Jacobian or the derivative of the
    % residual with respect to the variables
    
    dL2did2 = @(id,iq,l1) 2+l1*dTemg2did2(id,iq);
    dL2didiq = @(id,iq,l1) l1*dTemg2didiq(id,iq);
    dL2didl1 = @(id,iq,l1) dTemgdid(id,iq);
    
    
    dL2diq2 = @(id,iq,l1) 2+l1*dTemg2diq2(id,iq);
    dL2diql1 = @(id,iq,l1) dTemgdiq(id,iq);
    
    % Initialization
%     x0 = [-10;    % id
%           20;    % iq
%           10];  % l1
    
    res = [dLdid(x0(1),x0(2),x0(3))
           dLdiq(x0(1),x0(2),x0(3))
           dLdl1(x0(1),x0(2),x0(3))];
    
    jac = [dL2did2(x0(1),x0(2),x0(3))  dL2didiq(x0(1),x0(2),x0(3))  dL2didl1(x0(1),x0(2),x0(3))
           dL2didiq(x0(1),x0(2),x0(3)) dL2diq2(x0(1),x0(2),x0(3))   dL2diql1(x0(1),x0(2),x0(3))
           dL2didl1(x0(1),x0(2),x0(3)) dL2diql1(x0(1),x0(2),x0(3))  0];
    
    err = norm(res);
    
    ite = 0;
    tol = 1e-6;
    itemax = 100;
    relax = 0.5;
    
    while (ite<itemax && err>tol)
        ite = ite +1;
        dx = -jac\res;
        x0 = x0 +relax*dx;
        res = [dLdid(x0(1),x0(2),x0(3))
               dLdiq(x0(1),x0(2),x0(3))
               dLdl1(x0(1),x0(2),x0(3))];
    
        jac = [dL2did2(x0(1),x0(2),x0(3))  dL2didiq(x0(1),x0(2),x0(3))  dL2didl1(x0(1),x0(2),x0(3))
               dL2didiq(x0(1),x0(2),x0(3)) dL2diq2(x0(1),x0(2),x0(3))   dL2diql1(x0(1),x0(2),x0(3))
               dL2didl1(x0(1),x0(2),x0(3)) dL2diql1(x0(1),x0(2),x0(3))  0];
        
        err = norm(res);
    
%         fprintf('ite = %d \t err = %.6g \n',ite,err);
    end
    ID(ntorque) = x0(1);
    IQ(ntorque) = x0(2);
    l1 = x0(3);
    fprintf('\n Temg = %.3g Nm\t Id = %.3g A\t Iq = %.3g A\n',Temg(x0(1),x0(2)),x0(1),x0(2))
end
%     Temg(ID,IQ)

