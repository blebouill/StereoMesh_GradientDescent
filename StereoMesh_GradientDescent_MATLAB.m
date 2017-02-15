clear all; close all; clc;
dbstop if error;


%%  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  %
%                                  PATH                                   %
%  _____________________________________________________________________  %

main_path = 'C:\Users\blebouill\Documents\MATLAB\StereoMesh_GradientDescent';
addpath(genpath(main_path))
cd(main_path)





%%  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  %
%                                  DATA                                   %
%  _____________________________________________________________________  %

middlebury_dataset_training_path = 'F:\Datasets\Middlebury Stereo v3\trainingF';
addpath(genpath(middlebury_dataset_training_path))
cd(middlebury_dataset_training_path)

% Choix du fichier
listing = dir(middlebury_dataset_training_path);

cd(listing(14).name)    % Recycle

% Chargement des images
img_L = imread('im0.png');
img_R = imread('im1.png');

% figure,
% subplot(1,2,1), imshow(img_L)
% subplot(1,2,2), imshow(img_R)


% Chargement des paramètres
% --------------------- AUTOMATISER ---------------------

K_L = [2987.1 0 1307.34; 0 2987.1 958.889; 0 0 1];
K_R = [2987.1 0 1483.79; 0 2987.1 958.889; 0 0 1];

baseline = 178.232; %(mm)

ndisp = 260;            % Conservative bound on max disparity - CAN be used
disparity_min = 32;     % Minimum disparity - CANNOT be used
disparity_max = 227;    % Maximum disparity - CANNOT be used
dyavg = 0.394;          % Average disparity on y axis
dymax = 1.465;          % Max disparity on y axis





%%  ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  %
%                            PYRAMIDE D'IMAGES                            %
%  _____________________________________________________________________  %

downscaling_factor = 0.8;
min_hw = 10;   % Hauteur (ou largeur) min de l'image du dernier étage de la pyramide

[h, w, ~] = size(img_L);

a = min(h, w);
N_downscale = round( log(min_hw/a) / log(downscaling_factor) );

% Pyramide d'images
img_L_pyramid = cell(1, N_downscale);
img_R_pyramid = cell(1, N_downscale);

img_L_pyramid{1} = img_L;
img_R_pyramid{1} = img_R;

% Matrices de calibrations associées
K_L_pyramid = cell(1, N_downscale);
K_R_pyramid = cell(1, N_downscale);

K_L_pyramid{1} = K_L;
K_R_pyramid{1} = K_R;

for ind_pyramid = 2:N_downscale
    
    % __ Downsampled images (color)
    
    img_L_pyramid{ind_pyramid} = imresize(img_L_pyramid{ind_pyramid-1}, downscaling_factor, 'Antialiasing', true);
    img_R_pyramid{ind_pyramid} = imresize(img_R_pyramid{ind_pyramid-1}, downscaling_factor, 'Antialiasing', true);
    
    
    % __ Intrinsic parameters
    
    K_L_pyramid{ind_pyramid} = [downscaling_factor 0 0 ; 0 downscaling_factor 0 ; 0 0 1] * K_L_pyramid{ind_pyramid-1};
    K_R_pyramid{ind_pyramid} = [downscaling_factor 0 0 ; 0 downscaling_factor 0 ; 0 0 1] * K_R_pyramid{ind_pyramid-1};
    
end

% for ind_pyramid = 1:N_downscale
%     figure, imshow(img_R_pyramid{ind_pyramid})
% end





%%  ____________________________________________________________________  %
%   ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  %
%                    STEREOMESH USING GRADIENT DESCENT                    %
%   ____________________________________________________________________  %
%   ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯  %

StereoMesh_GradientDescent_path = 'C:\Users\blebouill\Documents\Visual Studio 2015\Projects\StereoMesh_GradientDescent\StereoMesh_GradientDescent';
addpath(genpath(StereoMesh_GradientDescent_path))
% cd(StereoMesh_GradientDescent_path)

for ind_pyramid = N_downscale:-1:5
    
    % ===============================================
    % == Récupération des images et des paramètres
    
    img_L = img_L_pyramid{ind_pyramid};
    img_R = img_R_pyramid{ind_pyramid};
    
    K_L = K_L_pyramid{ind_pyramid};
    K_R = K_R_pyramid{ind_pyramid};
    
    [height_img, width_img, ~] = size(img_L);
    
%     figure,
%     subplot(1,2,1), imshow(img_L)
%     subplot(1,2,2), imshow(img_R)

    
    % ===============================================
    % == Segmentation
    
    % __ Paramètres
%     nb_pixel_mean_triangle = 10;
%     nb_triangles = height_img * width_img / nb_pixel_mean_triangle;
    nb_triangles = height_img; 
    lambda = 600;
    energy_shape_definition = 1;
    nb_iter_max = height_img;


    % __ Filtre gaussien
    sigma_filtre = 0;
    
    if (sigma_filtre > 0)
        img_segmentation = imgaussfilt(img_L, sigma_filtre,  'Padding', 'symmetric');
    else
        img_segmentation = img_L;
    end


    % __ Segmentation
    [img_label, v_vertex, v_edge, v_triangle, ind_triangles_using_vertex, ind_vertex_in_triangle, ind_triangles_using_edge, ~] = ...
                        segmentation_triangle_MATLAB(img_segmentation, nb_triangles, lambda, energy_shape_definition, nb_iter_max);
    
    N_V = length(v_vertex(1,:));
    N_E = length(v_edge(1,:));
    N_T = max(img_label(:));
    
    
    
    % ===============================================
    % == StereoMesh_GradientDescent
    
    % __ Paramètres
    delta_init = 10^(2);
    lambda_BREACH = 10^(-3);
    lambda_NORMAL = 0; %10^(-3);
    N_ITERmax = 1000;   
    
    % __ Matrices S
    S = double(reshape([reshape(v_triangle,2, []) ; ones(1,3*N_T)], 3, 3, N_T));
    
    % __ Gradient de l'image
    h_y = fspecial('sobel');
    grad_img_R_x = imfilter(img_R, h_y');
    
%     figure, imshow(abs(double(grad_img_R_x))/255, []), colorbar
    

    % __ Initialisation de D
    
    if (ind_pyramid == N_downscale)
        
        % Initialisation à 0
        D_init = zeros(3*N_T, 1);
        
    else
        
        % Initialisation aux valeurs précédentes
        ind_pix_S = S(2,:,:) + (S(1,:,:)-1).*height_img;
        ind_pix_S = ind_pix_S(:);
        
        D_init = disparity_map_next_level(ind_pix_S);
    end



    % __ Descente de gradient
    
%     compile_MEX_StereoMesh_GradientDescent;

    [final_disparity_map, final_img2_interp, v_energy, D] = ...
                                               StereoMesh_GradientDescent( permute( double(img_L)/255.0, [2 1 3]), ...
                                                                           permute( double(img_R)/255.0, [2 1 3]), ...
                                                                           permute( double(grad_img_R_x)/255.0, [2 1 3]), ...
                                                                           K_L', ...
                                                                           baseline, ...
                                                                           img_label' - 1, ...
                                                                           N_V, ...
                                                                           N_E, ...
                                                                           N_T, ...
                                                                           ind_vertex_in_triangle - 1, ...
                                                                           ind_triangles_using_vertex - 1, ...
                                                                           ind_triangles_using_edge - 1, ...
                                                                           D_init, ...
                                                                           permute( S - [ones(2,3,N_T) ; zeros(1,3,N_T)], [2 1 3]), ...
                                                                           delta_init, ...
                                                                           lambda_BREACH, ...
                                                                           lambda_NORMAL, ...
                                                                           N_ITERmax);
    
                                                                       
                                                                       
    % ===============================================
    % == Affichages
    
%     figure, imshow(img_L), title('Image L')
%     figure, imshow(permute(final_img2_interp, [2 1 3])), title('Image R interpolée')
%     figure, imshow(img_R), title('Image R')
%     figure, imshow(abs(permute(final_img2_interp, [2 1 3]) - double(img_R)/255))
%     figure, imshow(final_disparity_map', []), colormap(jet), colorbar
    
%     figure, plot(v_energy(v_energy~=0)), title(['Evolution de l''energie [ ' num2str(ind_pyramid) ' ]' ])


    % ===============================================
    % == Interpolation de la carte des disparités pour l'étage inférieur
    
    if (ind_pyramid > 1)
        
        % Size next image
        [h_nxt, w_nxt, ~] = size(img_L_pyramid{ind_pyramid-1});
        
        % Upscaling factor
        upscaling_factor = w_nxt / width_img;
        
        % Upscale disparity map
        disparity_map_next_level = imresize(final_disparity_map' .* upscaling_factor, [h_nxt w_nxt], 'Antialiasing', true);
    
        % Show results
%         figure,
%         subplot(1,2,1), imshow(final_disparity_map', []), colormap(jet), colorbar
%         title(['Carte des disparités [ ' num2str(ind_pyramid) ' ]'])
%         subplot(1,2,2), imshow(disparity_map_next_level, []), colormap(jet), colorbar
%         title(['Carte des disparités pour l''initialisation de l''étage inférieur [ ' num2str(ind_pyramid-1) ' ]'])
    
    end
    
end

                                                                       
% ===============================================
% == Affichages
    
figure, imshow(img_L), title('Image L')
figure, imshow(permute(final_img2_interp, [2 1 3])), title('Image R interpolée')
figure, imshow(img_R), title('Image R')
figure, imshow(abs(permute(final_img2_interp, [2 1 3]) - double(img_R)/255))
figure, imshow(final_disparity_map', []), colormap(jet), colorbar