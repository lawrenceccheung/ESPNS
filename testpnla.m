clear all
addpath('PNLA_MATLAB_OCTAVE');
% polysys{1,1}=[1 -3]
% polysys{2,1}=[1 -1 1 -5]
% polysys{3,1}=[1 -2 7]
% polysys{1,2}=[1 1 0; 0 0 0]
% polysys{2,2}=[2 0 0; 0 0 2; 1 1 0; 0 0 0]
% polysys{3,2}=[3 0 0; 1 1 0; 0 0 0]

% Example in phd_kim
polysys{1,1}=[1 -3];
polysys{2,1}=[1 -1 1 -5];
polysys{3,1}=[1 -2 7];
polysys{1,2}=[1 1 0; 0 0 0];
polysys{2,2}=[2 0 0; 0 0 2; 1 0 1; 0 0 0];
polysys{3,2}=[0 0 3; 1 1 0; 0 0 0];

[root d c ns check cr digits] = sparf(polysys);
disp(root)

% The correct answer should be
%root =
%  -2.0000 + 0.0000i  -1.5000 + 0.0000i  -1.0000 + 0.0000i
%   1.8574 + 0.1762i   1.6008 - 0.1518i   0.5000 - 0.8660i
%   1.8574 - 0.1762i   1.6008 + 0.1518i   0.5000 + 0.8660i
%  -2.3574 + 0.6899i  -1.1722 - 0.3430i   0.5000 - 0.8660i
%  -2.3574 - 0.6899i  -1.1722 + 0.3430i   0.5000 + 0.8660i
%   3.0000 + 0.0000i   1.0000 + 0.0000i  -1.0000 + 0.0000i