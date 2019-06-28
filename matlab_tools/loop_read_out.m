function [data]=loop_read_out(name_exp)
% loop_init_state(init_state,w,theta,salt)
%
% Read output files created by the loop model.
%
% loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
% F. Roquet 2019
% GNU General Public License

name_out = sprintf('%s_output.nc',name_exp);
data = ncload_struct(name_out);

