function [  Fmat, gmat, PROBLEM,GAM0,GAM1,PSI,PPI ] = model_solveR( Theta )
% produce the states space solution of the model give by:
% GAM0*y(t) = GAM1*y(t-1) + C + PSI*z(t) + PPI*eta(t).
% and the solution is solved by sims gensys given by
% beta(t) = Fmat*beta(t-1) + gmat*z(t)

%extract parameters
sigma = Theta(1);
beta = 0.999;  %calibrated
delta = Theta(2);
alpha = Theta(3);
omega = Theta(4);
rho1 = 0;  %restricted
rho2 = Theta(5);
rho3 = Theta(6);
% parameter definitions:
kappa=((1-omega)*(1-(beta*omega)))/(alpha*omega);

%********************************************************************
%*      matrices of canonical system
%********************************************************************/

%* Define equation indices **/
eq_IS   = 1;   %* IS curve**/
eq_PHIL       = 2;   %* Phillips Curve **/
eq_RULE       = 3;   %* Monetary Policy Rule **/
eq_g      = 4;   %* AR for g **/
eq_u   = 5;   %* AR process for u **/
eq_v     = 6;   %* AR process for v **/
eq_Ex     = 7;   %* AR process for v **/
eq_Epi     = 8;   %* AR process for v **/
%* variable indices **/
v_x     = 1;
v_pi         = 2;
v_R          = 3;
v_g           =4;
v_u           =5;
v_v           =6;
v_Ex       = 7;  %* E[x] **/
v_Epi       = 8;  %* E[pi] **/

%* shock indices **/
e_x     = 1;
e_pi = 2;
e_i = 3;
%* expectation error indices **/
n_x      = 1;
n_pi        = 2;

%* summary **/
neq  = 8;  %number of equations
neps = 3;  %number of shocks
neta = 2;  %number of expectational errors

%* initialize matrices **/
GAM0 = zeros(neq,neq);
GAM1 = zeros(neq,neq);
   C = zeros(neq,1);
 PSI = zeros(neq,neps);
 PPI = zeros(neq,neta);

%% equations
% GAM0 y(t) = GAM1 y(t-1) + C + PSI z(t) + PPI eta(t) 
%*********************************************************
%**      1. IS Curve
%*********************************************************/
GAM0(eq_IS,v_x)  =  1;
GAM0(eq_IS,v_Ex)  =  -1;
GAM0(eq_IS,v_R)  =  1/sigma;
GAM0(eq_IS,v_Epi)  =  -1/sigma;
GAM0(eq_IS,v_g)  =  -1;


%*********************************************************
%*      2. Phillips Curve
%**********************************************************/
GAM0(eq_PHIL,v_pi)=1;
GAM0(eq_PHIL,v_Epi)=-beta;
GAM0(eq_PHIL,v_x)=-kappa;
GAM0(eq_PHIL,v_u)=-1;

%*********************************************************
%**      3. Policy Rule
%**********************************************************/
GAM0(eq_RULE,v_R)=1;
GAM0(eq_RULE,v_pi)=-delta;
GAM0(eq_RULE,v_v)=-1;
 
%*********************************************************
%**      Shock process
%**********************************************************/
%* g **/
GAM0(eq_g,v_g) = 1;
GAM1(eq_g,v_g) = rho1;
 PSI(eq_g,e_x) = 1;
%* u **/
GAM0(eq_u,v_u) = 1;
GAM1(eq_u,v_u) = rho2;
 PSI(eq_u,e_pi) = 1;
%* v* **/
GAM0(eq_v,v_v) = 1;
GAM1(eq_v,v_v) = rho3;
 PSI(eq_v,e_i) = 1;

%*********************************************************
%**      Expectation error
%**********************************************************/

% E(x)
GAM0(eq_Ex,v_x)  = 1;
GAM1(eq_Ex,v_Ex) = 1;
 PPI(eq_Ex,n_x)  = 1;
%* E(pi) **/
GAM0(eq_Epi,v_pi)  = 1;
GAM1(eq_Epi,v_Epi) = 1;
 PPI(eq_Epi,n_pi)  = 1;

%*******************************************************************
%**      QZ(generalized Schur) decomposition by GENSYS
%********************************************************************/

PROBLEM=0;	
div = 1;        % criteria for stable roots
% y(t) = T1*y(t-1) + TC + T0*z(t)
% A, B, Q and Z are the QZ decomp matrix
% RC(1) = 1 => existence
% RC(2) = 1 => Unique
[T1,TC,T0,Q,A,B,~,RC,LOOSE] = gensys(GAM0,GAM1,C,PSI,PPI,div);

Fmat = zeros(neq-neta);
Fmat(1:neq-neta,1:neq-neta) = T1(1:neq-neta,1:neq-neta);


 gmat= [T0(1:neq-neta,:)];

if RC(1)~=1 | RC(2)~=1  % non-unique and existence
    PROBLEM=1;
end

