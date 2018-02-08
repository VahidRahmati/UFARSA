% After unzipping the UFARSA.zip into your desired directory, simply run this
% script to add UFARSA to the search path of MATALB.
%
% Author: Vahid Rahmati (December, 2017)

UFARSA_folder = fileparts(mfilename('fullpath'));
addpath(sprintf('%s', UFARSA_folder))
addpath(sprintf('%s%sfunctions', UFARSA_folder, filesep));
addpath(sprintf('%s%sfunctions%sBEADS', UFARSA_folder, filesep, filesep));
addpath(sprintf('%s%sfunctions%sGetSn', UFARSA_folder, filesep, filesep));
addpath(sprintf('%s%sfunctions%sSMOOTHN', UFARSA_folder, filesep, filesep));
savepath