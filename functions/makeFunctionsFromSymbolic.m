function makeFunctionsFromSymbolic(b, b_m, str_psai, str_psai_sym)


% save('data_final','f', 'b')
% fid = fopen('force_field_symbolic.m', 'wt');
% fprintf(fid, 'function [f] = force_field_symbolic(x,y,psai)\n');
% fprintf(fid, '%s', str_psai{:});
% fprintf(fid, 'f=[\n');
% fprintf(fid, '\t%s\n', char(f(1)));
% fprintf(fid, '\t%s\n', char(f(2)));
% fprintf(fid, '\t%s\n', char(f(3)));
% fprintf(fid, '];\nend');
% fclose(fid);

fid = fopen('force_field_symbolic.m', 'wt');
fprintf(fid, 'function [f] = force_field_symbolic(x,y,psai,outputId)\n');
fprintf(fid, 'if nargin > 2\n');
fprintf(fid, '\t%s', str_psai{:});
fprintf(fid, 'else\n');
fprintf(fid, '\t%s', str_psai_sym{:});
fprintf(fid, 'end\n');
fprintf(fid, 'h = 1e-4;\n');
fprintf(fid, 'z = h;\n');
fprintf(fid, 'b1=[\n');
fprintf(fid, '\t%s\n', char(b_m(1)));
fprintf(fid, '\t%s\n', char(b_m(2)));
fprintf(fid, '\t%s\n', char(b_m(3)));
fprintf(fid, '];\n');
fprintf(fid, 'z = -h;\n');
fprintf(fid, 'b2=[\n');
fprintf(fid, '\t%s\n', char(b_m(1)));
fprintf(fid, '\t%s\n', char(b_m(2)));
fprintf(fid, '\t%s\n', char(b_m(3)));
fprintf(fid, '];\n');
fprintf(fid, 'f = (b1 - b2)./2/h;\n');
fprintf(fid, 'if nargin > 3\n');
fprintf(fid, '\tif outputId < 4\n');
fprintf(fid, '\t\tf = f(outputId);\n');
fprintf(fid, '\telse\n');
fprintf(fid, '\t\tf = sqrt(f(1)^2+f(2)^2);\n');
fprintf(fid, '\tend\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end');
fclose(fid);





fid = fopen('magnetic_field_symbolic.m', 'wt');
fprintf(fid, 'function [b] = magnetic_field_symbolic(x,y,psai,z)\n');
fprintf(fid, 'if nargin<4\n');
fprintf(fid, '\tz=0;\n');
fprintf(fid, 'end\n');
fprintf(fid, '%s', str_psai{:});
fprintf(fid, 'b=[\n');
fprintf(fid, '\t%s\n', char(b(1)));
fprintf(fid, '\t%s\n', char(b(2)));
fprintf(fid, '\t%s\n', char(b(3)));
fprintf(fid, '];\nend');
fclose(fid);



end

