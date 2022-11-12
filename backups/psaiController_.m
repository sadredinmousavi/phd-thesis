function [Psai,target_x,target_y,targetInd,real_x,real_y] = psaiController(t)
Psai_0 = [
    30
    150
    -90
]*pi/180;
Psai_0 = Psai_0 + [
    0
    0
    -75
]*pi/180;
Psai = Psai_0;
if nargout > 1
	targetInd = ind;
	target_x = target_x_sequence;
	target_y = target_y_sequence;
	real_x = real_x_sequence(:,ind);
	real_y = real_y_sequence(:,ind);
end
end