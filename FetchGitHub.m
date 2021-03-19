%% FetchGitHub.m
%This script consists of a few command line Git commands for
%fetching/pulling/merging from the most recent version of VOGA.

%These first two lines really only needs to be run once but can't hurt
!git add * 
!git remote add upstream https://github.com/ayiotis/VOGA.git
!git fetch upstream
!git merge upstream/master
%Added these commands because of some errors I got with merging
!git commit
!git push