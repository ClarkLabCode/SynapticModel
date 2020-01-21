%% GenerateFlickerFigures.m
% Script to generate Supplementary Figure 5A from
%     Matulis, C.A., Chen, J., Gonzalez-Suarez, A.D., Behnia, R. and
%     Clark, D.A., 2020. Heterogeneous Temporal Contrast Adaptation in
%     Drosophila Direction-Selective Circuits. Current Biology. DOI:
%     https://doi.org/10.1016/j.cub.2019.11.077
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

[ config ] = SetConfiguration('regenerateData', true);

%% Run simulation with TF = 1 Hz, SF = 1/45 1/deg, and 150 ms time constants

ThreeInputModelSinusoidPlusFlickerResponses(config, 1, 1/45, 'tauLp', 0.15, 'tauHp', 0.15);

%% Run simulation with TF = 1 Hz, SF = 1/30 1/deg, and 150 ms time constants

ThreeInputModelSinusoidPlusFlickerResponses(config, 1, 1/30, 'tauLp', 0.15, 'tauHp', 0.15);

%% Run simulation with TF = 1 Hz, SF = 1/45 1/deg, and 200 ms time constants

ThreeInputModelSinusoidPlusFlickerResponses(config, 1, 1/45, 'tauLp', 0.2, 'tauHp', 0.2);

%% Run simulation with TF = 1 Hz, SF = 1/30 1/deg, and 200 ms time constants

ThreeInputModelSinusoidPlusFlickerResponses(config, 1, 1/30, 'tauLp', 0.2, 'tauHp', 0.2);

