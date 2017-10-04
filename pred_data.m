function pred_data(l,m,pred_ind,q_est,save_im,res_plot,stat_model,m0,m1,dm0,dm1)



    il = mod(l,4);
    if il == 0
        il = 4; %to determine grid size.
    end
    
     %load in data
    load('ind_cell_prof_data.mat')
   
    
    
    %choose data
    cell_data = squeeze(ind_cell_data{m-1,2}(pred_ind,:,16:end))';
    FRET_data = squeeze(ind_fret_data{m-1,2}(pred_ind,:,16:end))';
    cell_data_std = ind_cell_data_sv{m-1,2}(:,16:end);
    
   
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    tdata = 5:1/3:1/3*(tndata-1+5*3);
    %make ndgrid for interpolation
    [Xdata,Tdata] = ndgrid(xdata,tdata);
    
    %IC type
    IC_type = 'interp';
        
    
    %%%%% Fit linear interpolant to (t,x) FRET data
    F1 = griddedInterpolant(Xdata,Tdata,FRET_data');
    % the vector x will be constant from now on, so can call interpolated
    % FRET data by: F(x,t*ones(length(x),1))
      
    
%     m0 = 0;
%     m1 = max(max(FRET_data));
    
    %now, we can interpolate the derivative as well. First iterate through
    %each timepoint and smooth and take derivative:
    FRET_data_dx = zeros(size(FRET_data));
    for i = 1:tndata
        FRET_data_dx(i,:) = [diff(smooth(FRET_data(i,:)))' ...
            diff(FRET_data(i,end-1:end))]/(xdata(2)-xdata(1));
    end
    
     F2 = griddedInterpolant(Xdata,Tdata,FRET_data_dx');
    
%     dm0 = min(min(min(FRET_data_dx)));
%     dm1 = max(max(max(FRET_data_dx)));
    


    %%%% Now predict migration data
       
    %grid sizes and models considered
    xnsize = [25 50 100 200];
    
  

    %generate grids for computation
    xn = xnsize(il);
    dt = 1e-3;

    [x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
    x=x';
    tn = length(t);
    dx = x(2)-x(1);
    %interior, boundary points
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

    %initial condition
        switch IC_type
        case 'interp'

                IC = interp1(xdata,cell_data(1,:),x);
                cutoff_x = leading_edge_calc(IC,x,.05,0);
                IC(x>cutoff_x)=0;

        case 'step'

                LE_loc = leading_edge_calc(smooth(cell_data(1,:)),xdata,0.8,0);
                IC = double(x<=LE_loc);
        end


    %boundary conditions
    BC_x_0 = @(t) 1;
    BC_x_1 = @(t) 0;

    %sparse matrix as a function for computation

    A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],...
        [vw.*(-1+sw/2); (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],xn,xn);
    
    A_pos_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],...
        [ve.*(1-1*se/2); ve.*se/2],xn,xn);

    A_pos_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
        (-vw.*sw/2)],total,total);



    A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],...
        [(-vw.*sw/2); (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],xn,xn);

    A_neg_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
        (vw.*sw/2-vw)],xn,xn);

    A_neg_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],...
        [ve.*se/2; (ve-ve.*se/2)],xn,xn);


    tic
                
        toCompare = 6;
        cell_data = cell_data(1:3*toCompare:end,:);
        data_std = cell_data_std(:,1:3*toCompare:end);
        
        tdata = tdata(1:3*toCompare:end);
        
           
        name_save = [stat_model '_' num2str(toCompare)];
        
        [J,OLS_SV,res,model] = MLE_cost_D0_fewer_compare_space(cell_data...
                ,q_est,F1,F2,m0,m1,dm0,dm1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,...
                tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg);
      

    toc
    
    figure
    hold on
    
    colors = distinguishable_colors(length(tdata));
    
    if strcmp(stat_model,'fewer')
        
        for i = 1:length(tdata)
            plot(xdata,model(i,:),'color',colors(i,:))
            plot(xdata(1:3:end),cell_data(i,1:3:end),'.','color',[colors(i,:)],'markersize',15)
            fill([xdata(1:3:end) fliplr(xdata(1:3:end))],[(cell_data(i,1:3:end) +...
                data_std(1:3:end,i)') fliplr(cell_data(i,1:3:end) - data_std(1:3:end,i)')],colors(i,:),'Facealpha',0.1,'edgecolor','none')
        end
        
    else
    
        count = 1;
        
        for i = 1:floor(tndata/4):tndata
            plot(xdata,model(i,:),colors(count))
            plot(xdata,cell_data(i,:),[colors(count) '.'])
            
            count = count + 1;
        end
    end
    
%     legend('1','2','3','4','5','location','northeast')
    
    cell_dens = [1700,2500,3000,4000];

%     h=title(['Predicting cell density, density = ' num2str(cell_dens(m-1)) ' cells/mm$^2, x_n$ = ' num2str(xn)]);
    h=title(['Predicting, density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$']);
    set(h,'interpreter','latex')
    xlabel('Location (x)')
    ylabel('u(t,x)')
    
    if save_im == 1
        exportfig(gcf,['cell_predicting_' num2str(m) '_' num2str(l) name_save '_pred_' num2str(pred_ind) '.eps'],'color','rgb','fontsize',2)
        saveas(gcf,['cell_predicting_' num2str(m) '_' num2str(l) name_save '_pred_' num2str(pred_ind) '.fig'])
    end
    
%     close
    
%     figure
%     hold on
%     
%     plot(tdata,ppval(p,tdata),'b')
%     
%     if plot_data == 1
%         plot(tdata,FRET_data(1:tndata),'b.')
%         
%     end
%     
%     title(['FRET interpolation for cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$'],'interpreter','latex')
%     xlabel('Time (t)')
%     ylabel('FRET')
%     
%     if save_im == 1
%         exportfig(gcf,['FRET_interp_' num2str(m) '_' name save '.eps'],'color','rgb')
%         saveas(gcf,['FRET_interp_' num2str(m) '_' name save '.fig'])
%     end
    
%     close

    if res_plot == 1
       
       [X,T] = meshgrid(xdata,tdata);
       
       
       if strcmp(stat_model,'WLS')
            res_final = (res./weight_matrix)/sqrt(2*WLS_SV);
       elseif strcmp(stat_model,'autor')
           res_final = stat_options.Ainv*res/sqrt(N);
%             res_final = reshape(res_final,)
           figure

           surf(X,T,reshape(res_final,length(tdata),length(xdata)),'edgecolor','none')

           xlabel('x','fontsize',30)
           ylabel('t','fontsize',30)

           set(gca,'fontsize',20)

           title(['WLS residuals, cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$'],'interpreter','latex')
%            axis([0 1 0 tdata(end)])

           caxis([-2 2])
           colorbar
       
       
            if save_im == 1
                exportfig(gcf,['FRET_interp_predict_' num2str(m)   '_' num2str(l) '_' name save '_res.eps'],'color','rgb')
                saveas(gcf,['FRET_interp_predict_' num2str(m) '_' num2str(l) '_' name save '_res.fig'])
            end
    
       elseif strcmp(stat_model,'fewer')
               
           model_vec = model(:);
           
           res_final = res;

            figure
            hold on
            plot(res_final,'b.')
            plot([0 1e6],[0 0],'k')
             axis([0 length(model_vec) -.3 .3])
 
            xlabel('Data Points','fontsize',30)
            ylabel('$y_{ij} - f(t_i,x_j)$','fontsize',30,'interpreter','latex')
%             title(['Residual, cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2,\ x_n$ = ' num2str(xn)],'interpreter','latex')
            title(['Residual, cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$'],'interpreter','latex','fontsize',2)
            
            
            if save_im == 1
                exportfig(gcf,['FRET_interp_predict_' num2str(m)   '_' num2str(l) '_' name_save '_res_pred_' num2str(pred_ind) '.eps'],'color','rgb','fontsize',2)
                saveas(gcf,['FRET_interp_predict_' num2str(m) '_' num2str(l) '_' name_save '_res_pred_' num2str(pred_ind) '.fig'])
            end
            
            figure
            hold on
            plot(model_vec,res_final,'b.')
%             plot([-.1 1.1],[0 0],'k')
%             axis([0 1.01 -.2 .2])
 
            xlabel('$f(t_i,x_j)$','fontsize',30,'interpreter','latex')
            ylabel('$y_{ij} - f(t_i,x_j)$','fontsize',30,'interpreter','latex')
            title(['Predicted residual, cell density = ' num2str(cell_dens(m-1)) ' cells/mm$^2$'],'interpreter','latex')
           
            if save_im == 1
                exportfig(gcf,['FRET_interp_predict_' num2str(m)   '_' num2str(l) '_' name_save '_res_model_pred_' num2str(pred_ind) '.eps'],'color','rgb')
                saveas(gcf,['FRET_interp_predict_' num2str(m) '_' num2str(l) '_' name_save '_res_model_pred_' num2str(pred_ind) '.fig'])
            end
       end
        
    end
    
    %%% plot velocity data from data, predicted
    
       
    %get finer data
    cell_data = squeeze(ind_cell_data{m-1,2}(pred_ind,:,16:3:end))';
    %estimate of how fast population migrated over time
    data_LE=leading_edge_calc(cell_data,xdata,0.5,0);
    
    tdata = 0:size(cell_data,1)-1;
%     
%     %compute v(m)
%     FRET = @(t) ppval(p,t);
%     n=5;
%     msamp = augknt([m0,m1,linspace(m0,m1,n)],2);
%     %create spline functions
%     v_spline = spmak(msamp,q_est');
%     %evaluate V(m(t))
%     Vx_c = @(t) fnval(v_spline,FRET(t));
%     
%     %plot v(m)
%     figure
%     plot(linspace(m0,m1,100),fnval(v_spline,linspace(m0,m1,100)))
%     title(['Estimated v(m) for a density of '  num2str(cell_dens(m-1)) ' cells/mm^2'])
%     xlabel('m')
%     ylabel('speed')
%     
%     if save_im == 1
%         exportfig(gcf,['FRET_fit_vm_' num2str(m) '_' num2str(l) '_' name_save '_pred_' num2str(pred_ind) '.eps'],'color','rgb')
%         saveas(gcf,['FRET_fit_vm_' num2str(m) '_' num2str(l) '_' name_save '_pred_' num2str(pred_ind) '.fig'])
%     end
%     
%     %plot FRET of predicted data set, along with resulting speed function
%     figure
%     hold on
%     plot(tdata,FRET_data(1:3:end)/max(FRET_data(1:3:end)))
%     plot(tdata(1:end-1),Vx_c(tdata(1:end-1))/max(Vx_c(tdata(1:end-1))))
%     legend('FRET ratio','Predicted Speed','location','northwest')
% 
%     xlabel('Time (t)')
%     ylabel('FRET')
%     title(['FRET ratio and predicted v(m) for a density of '  num2str(cell_dens(m-1)) ' cells/mm^2'])
%     
%     
%     if save_im == 1
%         exportfig(gcf,['FRET_pred_vm_' num2str(m) '_' num2str(l) '_' name_save '_pred_' num2str(pred_ind) '.eps'],'color','rgb')
%         saveas(gcf,['FRET_pred_vm_' num2str(m) '_' num2str(l) '_' name_save '_pred_' num2str(pred_ind) '.fig'])
%     end
%     
% %     subplot(2,2,1) %plot data velocity , predicted velocity, and FRET
% %     hold on
% %     plot(tdata(1:end-1),diff(data_LE)/(tdata(4)-tdata(1))/max(diff(data_LE)/(tdata(4)-tdata(1))))
% %     plot(tdata,FRET_data(1:3:end)/max(FRET_data(1:3:end)))
% %     plot(tdata,Vx_c(tdata)/max(Vx_c(tdata)))
% %     
% %     title('normalized Vdata, Vmodel, FRET')
% %     legend('Vdata','Vmodel','FRET','location','southeast')
% %     
% %     %plot v(m)
% %     subplot(2,2,2)
% %     
% %     
% %     subplot(2,2,4)
% %     
% %     plot(tdata,FRET_data(1:3:end))
% %     title('FRET')
    
end