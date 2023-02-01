% Calculate the boxes for material properties
clc,clear;
close all;
alp = 29;
n = [-sind(alp), cosd(alp)];

Q = [cosd(alp), -sind(alp);sind(alp), cosd(alp)];

rec1 = [0.006453, 0.003522; 0.063204, 0.035034];

rec2 = [0.005480, 0.003037; 0.064079, 0.035519];

% Find the points
pts = [rec1; rec2];

pts_up = pts;
pts_dn = pts;

lims = 0.101;
% for i = 1:1:size(pts_up, 1)
%     pts_up(i, :) = (lims - pts(i, 2)) / n(2) * n + pts(i, :);
%     pts_dn(i, :) = (-lims - pts(i, 2)) / n(2) * n + pts(i, :);
% end
% 
% for i = 1:1:size(pts_up, 1)
%     fprintf("    %1.24f", pts_dn(i, :));
%     fprintf("\n");
% end

% d1 = 0.001;
% 
d2 = 0.0012;
% format("%1.6f");
fprintf("    %1.6f", rec1 - n * d2);
% 
% % Material properties
% cs = 1.8e3; 
% cp = 4.0e3;
% G = 6e9;
% 
% rho = G / cs^2;
% K = rho * (3 * cp^2 - 4 * cs^2) / 3;
% la = rho * cp^2 - 2 * G;
% 
% shit = Q * [-0.3, 0.3]';
% % fprintf("    %1.6f", shit);
% 
% 
% Rec_ins = [0.005968,    0.004397; ... 
%            0.006938,    0.002647; ...
%            0.062719,    0.035909; ...
%            0.063689,    0.034159]';
% 
% Rec_ins_transferred = Q' * Rec_ins;
% 
% my = Rec_ins_transferred;
% % for i = 1:1:size(Rec_ins_transferred, 2)
% %     fprintf("    %1.6f", Rec_ins_transferred(:, i));
% %     fprintf("\n");
% % end
% 
% % Shit = [my(1, 1), -0.3; my(1, 1), 0.3; 
% %         my()]
