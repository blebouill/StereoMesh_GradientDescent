clc

% Current folder
currentFolder = pwd;

% Go to project directory
cd('C:\Users\blebouill\Documents\Visual Studio 2015\Projects\StereoMesh_GradientDescent\StereoMesh_GradientDescent');

% Compile MEX
mex -O -v -L. -output StereoMesh_GradientDescent MEX_StereoMesh_GradientDescent.cpp

% Go back to current folder
cd(currentFolder);
