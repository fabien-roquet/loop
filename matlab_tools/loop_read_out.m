function [out_0d,out_1d]=loop_read_out(name_exp)
% loop_init_state(init_state,w,theta,salt)
%
% Read output files created by the loop model.
%
% loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
% F. Roquet 2016
% GNU General Public License

name_out_0d = sprintf('%s_out_0D.txt',name_exp);
name_out_1d = sprintf('%s_out_1D.txt',name_exp);

[niter time w mass]=textread(name_out_0d,'%d%f%f%f','headerlines',1);
out_0d.niter = niter;
out_0d.time  = time;
out_0d.w     = w;
out_0d.mass  = mass;

[niter jk theta salt sigma]=textread(name_out_1d,'%f%d%f%f%f','headerlines',1);
nl=max(jk); nk=length(niter)/nl;
out_1d.niter = reshape(niter,nl,nk);
out_1d.theta = reshape(theta,nl,nk);
out_1d.salt  = reshape(salt,nl,nk);
out_1d.sigma = reshape(sigma,nl,nk);

