%% GenerateFigures.m
% Script to generate all synaptic model figures. 
%
%
%     Copyright (C) 2019 Clark Lab, Yale University
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
%% Set MATLAB path

addpath('analyses','models','stimuli','utils');

%% Set configuration
% Note: this function will need to be edited, or appropriate arguments
% provided, to set paths for caching and natural scene data access
% appropriately. 

[ config ] = SetConfiguration();

%% Sinusoid timetraces for each step

ThreeInputModelSinusoidTimetraces();

%% ON and OFF 30 deg/s moving edges

ThreeInputModelMovingEdgeResponses(config);

%% Sinusoid linearity analysis

ThreeInputModelSinusoidLinearityAnalysis(config);

%% Sinusoid bar plots with one TF-SF pairing

ThreeInputModelSinusoidResponses(config, 'numRep', 100);

%% Monocomponent sinusoid frequency sweeps

ThreeInputModelMonocomponentSinusoidPowerMap(config);

%% Oriented sinusoid grating responses

ThreeInputModelOrientedGratingResponses(config);

%% Local flashed ON and OFF bars

ThreeInputModelLocalFlashedBarResponse();

%% Bar pair apparent motion

ThreeInputModelBarPairApparentMotionResponses(config);

%% Flashed apparent motion

ThreeInputModelFlashedApparentMotionResponses();

%% Linear receptive field

ThreeInputModelLinearReceptiveField(config);

%% Ternary correlated noise (correlation interval receptive field)

ThreeInputModelTernaryCorrelatedNoiseResponses(config);

%% Two-point gliders

ThreeInputModelTwoPointGliderResponses(config);

%% Three-point gliders

ThreeInputModelThreePointGliderResponses(config);

%% Random checkerboards

ThreeInputModelRandomCheckerboardResponses(config);

%% Natural scenes

ThreeInputModelNaturalSceneResponses(config, 'numRep', 1e6);

