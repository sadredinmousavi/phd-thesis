function [Psai,eqPoints,targetInd,real_eqPoint_x,real_eqPoint_y] = psaiController(t)
time = [0.000000 20.000000 40.000000 70.000000 100.000000 ];
temp = time - t;
for k = 1:length(temp)
	if temp(k) > 0
		ind = k-1;
		break;
	elseif k == length(temp)
		ind = k;
	end
end
eqPoints_sequence{1} = [
	-0.050000 -0.050000 -0.030000 0.000000 -0.000000 
	0.080000 0.070000 0.060000 0.050000 0.030000 

];
eqPoints_sequence{2} = [
	0.050000 0.050000 0.030000 0.000000 0.000000 
	0.080000 0.070000 0.060000 0.050000 0.030000 

];
Psai_sequence = [
	-0.000000 -0.000000 1.543008 -0.000000 -0.034147 
	-0.000000 -0.000000 1.539016 -0.000000 -1.060072 
	-0.000000 -0.000000 1.418614 -0.000000 -0.984612 
	-0.000000 -0.000000 1.542289 -0.000000 -0.842966 
	-0.000000 -0.000000 -1.527162 -0.000000 -0.646063 
	-0.000000 -0.000000 1.350364 -0.000000 -0.007556 
	-0.000000 -0.000000 0.904560 -0.000000 0.136173 
	-0.000000 -0.000000 -1.317409 -0.000000 0.175171 

];
real_eqPoints_x_sequence = [
	-0.050000 -0.050000 -0.030000 0.000000 -0.000000 
	0.050000 0.050000 0.030000 0.000000 0.000000 

];
real_eqPoints_y_sequence = [
	0.080000 0.070000 0.060000 0.050000 0.030000 
	0.080000 0.070000 0.060000 0.050000 0.030000 

];
Psai = Psai_sequence(:,ind);
if nargout > 1
	targetInd = ind;
	eqPoints = eqPoints_sequence;
	real_eqPoint_x = real_eqPoints_x_sequence(:,ind);
	real_eqPoint_y = real_eqPoints_y_sequence(:,ind);
end
end