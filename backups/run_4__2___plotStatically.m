%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 

 
 % run symboic.m and f -> force_field_symbolic.m
 
 % todo MagPosR -> input to symbolic.m ; automatically generate force_field_symbolic.m
clc, close all, clearvars -except cursor_info
load('input_data')
% load('cursor_info')
%%
% num = 20;
num = 0;
maxPsaiD = 360;
for zz = 1:num+1
    % [ Rho Theta Phi Psai]
    counter = ((zz-1)/num)*maxPsaiD*(pi/180);
%     Psai = [
%         counter
%         0
%         counter
%         0
%         ];
    Psai = [
        0
        0
        0
        0
        ]*(pi/180);
%     Psai=0.5646
    printFig('input_data', Psai, [num2str(zz) 'a']);
    
    figure
    surf(x,y,F)
    
end