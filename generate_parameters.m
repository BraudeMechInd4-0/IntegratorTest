function [xpdotp, Re, J2,mu, whichconsts, consts, finals, dAT] = generate_parameters()
%% function generate_parameters()
% a function to generate all parameters and constants. 
% these parameters and lines repeat themselves in most my scripts and I am
% tired of copy/pasting them. a function is more elegant

xpdotp = 1440.0 / (2.0*pi); 
%mu = 398600.5;
%Re = 6378.137;
%J2 = 0.00108262998905;
%j2 = J2;
whichconsts = 84;
[tumin, mu, Re, xke, j2, j3, j4, j3oj2] = getgravc(whichconsts);
consts = [tumin, mu, Re, xke, j2, j3, j4, j3oj2];
J2 = j2;

woptions = weboptions;
woptions.Timeout=30;
websave('finals2000A.txt','https://datacenter.iers.org/data/latestVersion/finals.all.iau1980.txt',woptions);
%websave('finals2000A.txt','https://datacenter.iers.org/data/9/finals2000A.all',woptions);
%websave('finals2000A.txt','https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt'); %load the information about dT1 and dAT
%websave('finals2000A.txt','https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt'); %load the information about dT1 and dAT
%datfileh = fopen(fullfile('\', 'finals2000A.txt'));
datfileh = fopen('finals2000A.txt');
fseek(datfileh, 0,'eof');
filelength = ftell(datfileh);
fclose(datfileh);
finals2000A = importfile('finals2000A.txt', 1, filelength);
finals = finals2000A;
save finals finals
dAT = 35;