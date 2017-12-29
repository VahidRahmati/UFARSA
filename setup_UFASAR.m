% After unzipping the UFASAR.zip into your desired directory, simply run this
% script to add UFASAR to the search path of MATALB.
%
% Author: Vahid Rahmati (December, 2017)

UFASAR_folder = fileparts(mfilename('fullpath'));
addpath(sprintf('%s', UFASAR_folder))
addpath(sprintf('%s%sfunctions', UFASAR_folder, filesep));
addpath(sprintf('%s%sfunctions%sBEADS', UFASAR_folder, filesep, filesep));
addpath(sprintf('%s%sfunctions%sGetSn', UFASAR_folder, filesep, filesep));
addpath(sprintf('%s%sfunctions%sSMOOTHN', UFASAR_folder, filesep, filesep));
savepath