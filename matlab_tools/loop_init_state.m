function loop_init_state(file_init_state,w,theta,salt)
% loop_init_state(init_state,w,theta,salt)
%
% Create initial file named after init_state to initialize the loop model.
%
% loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
% F. Roquet 2016
% GNU General Public License

nl = length(theta);
fid = fopen(file_init_state,'w');
fprintf(fid,'%g\n',w);
for kk=1:nl,
    fprintf(fid,'%g %g %g\n',kk,theta(kk),salt(kk));
end

