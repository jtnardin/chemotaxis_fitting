%find_best_sim.m written 3-13-17 by JTN to find
%from multiple best-fit simulations, the one that best-fit 
%the data

clear all; clc

pred_ind = num2str(3);

for well = 3:5

        welllet = 'mean';
       
        filename = ['FRET_interp_est_well_' num2str(well) '_'];
        filename2 = ['_fewer_6_pred_' pred_ind '_space'];
%         %load in data
%         load('cell_data_1d_struct_mod.mat')
% 
%         %choose data
%         data = cell_data_1d_mod{5,well-1}';
%         %data grids
%         xdata = linspace(0,1,540);
%         tdata = 0:1/3:1/3*(size(data,1)-1);

        J_final = zeros(4,1);
        q_final = cell(4,1);
%         WLS_weights_final = cell(6,1);
%         weight_matrix_final = cell(6,1);

        for i = 1:4 %do for all models and grid sizes

            
            if exist([filename num2str(i) filename2 '.mat'],'file')
                load([filename num2str(i) filename2 '.mat'])

                best_fit = J_all == min(J_all);

                J_final(i) = J_all(best_fit);
                q_final{i} = q_all{best_fit};
            end
        end

        save([filename filename2 '_final.mat'],...
            'J_final','q_final','F1','F2','m0','m1','dm0','dm1')
        

end