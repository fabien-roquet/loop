% Script showing usage of the loop toolbox
%  1) compute equilibrium loop properties for the test configuration
%  2) analyze linear stablity
%  3) create init_state.txt file to init the loop model
%  4) load model outputs and compare with theoretical values
%
% loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
% F. Roquet 2019
% GNU General Public License

% test configuration
nl=360; Zf=.5; foldtrue=1; R=.1; lambda=0; mu=0; ref_t=5; xi_t=1; F_s=.15; tau=0;

% plot torque curves
w_list = (-1:.01:1)';
torque = w_list*NaN;
for kk=1:length(w_list)
    torque(kk) = loop_equilibrium(w_list(kk),R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s);
end

% determine semi-analytically steady-state
func = @(w)(w - loop_equilibrium(w,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s) - tau);
w1 = fzero(func,-.2); torque1 = w1 - tau;
w2 = fzero(func, .3); torque2 = w2 - tau;
w3 = fzero(func, .1); torque3 = w3 - tau;

% determine the jacobian and its eigenvalues for steady-states 1
jacobian1 = loop_jacobian(w1,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s);
EIG1=eig(jacobian1);lEIG1=EIG1(real(EIG1)>0);
disp(sprintf('STATE 1, w=%4.2f: no eigenvalues with positive real part --> stable state',w1));

% determine the jacobian and its eigenvalues for steady-states 2
jacobian2 = loop_jacobian(w2,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s);
EIG2=eig(jacobian2);lEIG2=EIG2(real(EIG2)>0);
disp(sprintf('STATE 2, w=%4.2f: two conjugate eigenvalues with positive real part %4.2f +/- %4.2fi --> oscillatory state',w2,real(lEIG2(1)),imag(lEIG2(1))));

% determine the jacobian and its eigenvalues for steady-states 3
jacobian3 = loop_jacobian(w3,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s);
EIG3=eig(jacobian3);lEIG3=EIG3(real(EIG3)>0);
disp(sprintf('STATE 3, w=%4.2f: one real positive eigenvalue %4.2f --> unstable state',w3,lEIG3));

% plot torque curves. Each intersection is a steady-state
figure(1),clf,hold on
plot(w_list,torque,w_list,w_list);
plot(w1,torque1,'xr',w2,torque2,'or',w3,torque3,'+r')
legend('buoyancy torque','\tau=0')
xlabel('velocity'),ylabel('torque')
title('Torque curves');

% create init_file.txt with theta and salt distribution of state 1
[torque1,tau1,theta1,salt1] = loop_equilibrium(w1,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s);
loop_init_state('init_state.txt',w1,theta1,salt1);

% read output files of the loop model and compare loop outputs with theory
data=loop_read_out('TEST');
[torque1,tau1,theta1,salt1] = loop_equilibrium(w1,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s);
disp(' ');
disp(sprintf('Difference between loop output and theoretical solution of STATE1 velocity   : %4.2g',data.w(end)-w1));
disp(sprintf('STD between loop output and theoretical solution of STATE1 theta distribution: %4.2g',std(data.theta(:,end)-theta1)));
disp(sprintf('STD between loop output and theoretical solution of STATE1 salt  distribution: %4.2g',std(data.salt(:,end)-salt1)));
figure(2), plot(1:nl,data.theta(:,end),1:nl,data.salt(:,end)),legend('STATE1 theta','STATE1 salt')
