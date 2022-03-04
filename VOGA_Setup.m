%% Run this script every time you pull a new version of VOGA
% Moves a copy of addVOGA.m to the userpath
addVOGA
copyfile('addVOGA.m',userpath)
% Set the new version or make the VOGA_VerInfo.txt if it doesn't exist
current_ver = 'v4.6.1';
VOGA__setVersion(current_ver);
VOGA__saveLastUsedParams;