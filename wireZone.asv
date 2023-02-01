clc,clear;
WirePos1 = [-0.025657; -0.014222; 0]';

alp = 29;
tangent = [cosd(alp), sind(29), 0];

width = 16e-3 / 2;
edge = 4e-3;

leftIn = WirePos1 - tangent * width;
leftOut = leftIn - tangent * edge;

rightIn = WirePos1 + tangent * width;
rightOut = rightIn + tangent * edge;

disp(num2str(leftOut, '%.6f  '));
disp(num2str(leftIn, '%.6f  '));
disp(num2str(rightIn, '%.6f  '));
disp(num2str(rightOut, '%.6f  '));