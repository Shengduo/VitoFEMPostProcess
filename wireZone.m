clc,clear;
WirePos1 = [-0.025657, -0.014222, 0];
VSstart = [0.006354, 0.003522, 0];
alp = 29;
tangent = [cosd(alp), sind(29), 0];
WirePos1_1 = WirePos1 - 0.001 * tangent;

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

% Calculate a rectangle that has linear slip weakening properties
edge_width = 0.2e-3;
rect_center = VSstart + 15e-3 * tangent;
inner_length = 8e-3; % Use to be 12
inner_width = 6e-3;
inner_rectangle = [rect_center(1:2) - inner_length / 3 * tangent(1:2), - inner_width / 2; 
                   rect_center(1:2) - inner_length / 3 * tangent(1:2), + inner_width / 2; 
                   rect_center(1:2) + inner_length / 3 * 2 * tangent(1:2), + inner_width / 2; 
                   rect_center(1:2) + inner_length / 3 * 2 * tangent(1:2), - inner_width / 2];

% Calculate outer triangle
inner_length = inner_length + 2 * edge_width;
inner_width = inner_width + 2 * edge_width;

outer_rectangle = [rect_center(1:2) - inner_length / 3 * tangent(1:2), - inner_width / 2; 
                   rect_center(1:2) - inner_length / 3 * tangent(1:2), + inner_width / 2; 
                   rect_center(1:2) + inner_length / 3 * 2 * tangent(1:2), + inner_width / 2; 
                   rect_center(1:2) + inner_length / 3 * 2 * tangent(1:2), - inner_width / 2];

dlmwrite('./SWPatch.txt', [inner_rectangle; outer_rectangle], 'Delimiter', '\t', 'precision', '%.6f');

