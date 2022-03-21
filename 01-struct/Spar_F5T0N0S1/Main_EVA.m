%% Initialization
clear all; close all; %clc;
restoredefaultpath;
addpath(genpath('C:/Work/FAST/matlab-toolbox/'))

A=[
    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  1.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  1.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00
    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  1.000E+00  0.000E+00  0.000E+00  0.000E+00
    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  1.000E+00  0.000E+00  0.000E+00
    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  1.000E+00  0.000E+00
    0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  0.000E+00  1.000E+00
    0.000E+00  0.000E+00  3.240E-02 -4.073E+01 -1.793E-01  7.662E-09  0.000E+00  0.000E+00 -1.855E-03  1.863E-03 -2.126E-01  1.008E-08
    0.000E+00  0.000E+00  4.078E+01 -3.242E-02 -1.834E-01 -1.510E-06  0.000E+00  0.000E+00 -1.380E-06  2.088E-02  1.037E-04 -6.681E-07
    0.000E+00  0.000E+00 -4.141E-01  6.515E-04  1.866E-03  1.397E-08  0.000E+00  0.000E+00  8.483E-08  1.722E-03  1.566E-05  5.843E-09
    0.000E+00  0.000E+00 -4.066E-04 -4.139E-01 -1.821E-03 -2.517E-10  0.000E+00  0.000E+00 -1.891E-05 -1.890E-03 -2.176E-03 -9.949E-11
    0.000E+00  0.000E+00  4.057E-02 -2.741E-05  1.780E-04  3.599E-08  0.000E+00  0.000E+00  8.197E-06  2.187E-01  1.890E-03  2.086E-08
    0.000E+00  0.000E+00  4.144E-01 -2.978E-04  1.819E-03 -2.677E-07  0.000E+00  0.000E+00 -1.984E-07  1.868E-04  1.264E-06 -4.701E-07
    ];

mbc = eiganalysis(A);
mbc.NaturalFreqs_Hz
mbc.DampedFreqs_Hz
mbc.DampRatios
