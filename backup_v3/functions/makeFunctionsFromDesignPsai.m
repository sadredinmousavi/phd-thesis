function makeFunctionsFromDesignPsai(time, eqP, Psai_t, x_eqPoints, y_eqPoints)

fid = fopen('psaiController.m', 'wt');
fprintf(fid, 'function [Psai,eqPoints,targetInd,real_eqPoint_x,real_eqPoint_y] = psaiController(t)\n');
fprintf(fid, 'time = [');
fprintf(fid, '%f ', time);
fprintf(fid, '];\n');
fprintf(fid, 'temp = time - t;\n');
fprintf(fid, 'for k = 1:length(temp)\n');
fprintf(fid, '\tif temp(k) > 0\n');
fprintf(fid, '\t\tind = k-1;\n');
fprintf(fid, '\t\tbreak;\n');
fprintf(fid, '\telseif k == length(temp)\n');
fprintf(fid, '\t\tind = k;\n');
fprintf(fid, '\tend\n');
fprintf(fid, 'end\n');
for i=1:length(eqP)
    fprintf(fid, 'eqPoints_sequence{%d} = [\n', i);
    for k=2:3
        fprintf(fid, '\t');
        fprintf(fid, '%f ', eqP{i}(k,:));
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n];\n');
end
fprintf(fid, 'Psai_sequence = [\n');
for i=1:size(Psai_t,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', Psai_t(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'real_eqPoints_x_sequence = [\n');
for i=1:size(x_eqPoints,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', x_eqPoints(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'real_eqPoints_y_sequence = [\n');
for i=1:size(y_eqPoints,1)
    fprintf(fid, '\t');
    fprintf(fid, '%f ', y_eqPoints(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n];\n');
fprintf(fid, 'Psai = Psai_sequence(:,ind);\n');
fprintf(fid, 'if nargout > 1\n');
fprintf(fid, '\ttargetInd = ind;\n');
%     fprintf(fid, '\teqPoints = eqPoints(:,ind);\n');
fprintf(fid, '\teqPoints = eqPoints_sequence;\n');
fprintf(fid, '\treal_eqPoint_x = real_eqPoints_x_sequence(:,ind);\n');
fprintf(fid, '\treal_eqPoint_y = real_eqPoints_y_sequence(:,ind);\n');
fprintf(fid, 'end\n');
fprintf(fid, 'end');
fclose(fid);
    
end

