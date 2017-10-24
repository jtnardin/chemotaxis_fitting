% clear all; clc

pred_ind = 3;


for well = 4


    filename = ['FRET_interp_est_well_' num2str(well) '_'];
    filename2 = ['_fewer_6_pred_' num2str(pred_ind) '_space_final'];
    
    for sim = 3

        load([filename filename2 '.mat'])

        save_im = 0;
        plot_res = 0;
        stat_model = 'fewer';


        %options for each type of simulation

         run_sim(sim,well,pred_ind,q_final{sim},save_im,plot_res,stat_model,F1,F2,m0,m1,dm0,dm1)

%          pred_data(sim,well,pred_ind,q_final{sim},save_im,plot_res,stat_model,m0,m1,dm0,dm1)

    end
          
end