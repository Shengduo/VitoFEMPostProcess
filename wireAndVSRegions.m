clc,clear;
W1 = 0.085;
VS_start = W1 + 36.6e-3;
VS_preStart = VS_start - 1e-3;
VS_end = VS_start + 60e-3;
VS_afterEnd = VS_end + 1e-3;
leftBd = [-0.1, -0.05543];
rightBd = - leftBd;
dir = -(leftBd - rightBd);
dir = dir ./ norm(dir, 2);

W1pos = leftBd + W1 * dir;
VS_preStartPos = leftBd + VS_preStart * dir;
VS_startPos = leftBd + VS_start * dir;
VS_afterEndPos = leftBd + VS_afterEnd * dir;
VS_endPos = leftBd + VS_end * dir;

W1LL = [-0.031034, -0.017202];
W1L = [-0.027586, -0.015291];
W1R = [-0.020690, -0.011468];
W1RR = [-0.017241, -0.009557];
W1Length = norm(W1R - W1L, 2);
W1LGap = norm(W1L - W1LL, 2);
W1RGap = norm(W1R - W1RR, 2);

W1Correct = (W1pos - (W1L + W1R) / 2);
W1LL = W1LL + W1Correct;
W1L = W1L + W1Correct;
W1R = W1R + W1Correct;
W1RR = W1RR + W1Correct;