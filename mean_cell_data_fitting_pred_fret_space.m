%mean_cell_data_fitting_pred_fret_space.m written  9-25-17 by JTN to perform
%OLS optimization for the equation u_t + v(t)u_x =0 in a parallelized
%for loop to mean experimental data

%in this script, v(t) = v_1(m) + v_2(dm/dx) to look into chemokinesis and
%chemotaxis, where m(t,x) is interpolated FRET ratio data

%last updated by JTN : 10-4-17 to use csaps for FRET interpolation. Less
%accurate to data, but smoother and less sensitive to the jaggedness.

function mean_cell_data_fitting_pred_fret_space(l,m,pred_ind,stat_model,simnum)

    %for bookkeeping
    well = m;
    
    %interp IC
    IC_type = 'interp';
    
        
    %select which  grid size we're using
    xi = mod(l,4);
    if xi == 0
        xi = 4; %to determine grid size.
    end
    
    %load in data
    load('ind_cell_prof_data.mat')
    
    %helps with FRET interpolation
    dx_large = 10;
    %choose data
    cell_data = avg_cell_data(ind_cell_data{m-1,2}(:,:,16:end),pred_ind)';
    FRET_data = avg_cell_data(ind_fret_data{m-1,2}(:,1:dx_large:end,16:end),pred_ind)';
    
        
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    xdata2 = xdata(1:dx_large:end);
    tdata = 5:1/3:1/3*(tndata-1+5*3);
    %make ndgrid for interpolation
    [Xdata,Tdata] = ndgrid(xdata,tdata);
    [Xdata2,Tdata2] = ndgrid(xdata2,tdata);

    %generate grids for computation, also helpful for interpolation
    xnsize = [25 50 100 200];
    xn = xnsize(xi);
    dt = 1e-3;
    [x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
    tn = length(t);
    dx = x(2)-x(1);
    %interior, boundary points
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);

    
    %find min, max FRET levels behind LE
    [m0,m1] = determine_FRET_behind_LE(xdata,xdata2,cell_data,FRET_data,tdata);
    

    %create time-dependent functions to interpolate FRET data and its
    %gradient. These functions can be called by "F1(x,t*ones(xn,1))" for any
    %t in the range of tdata
    [F1,F2,dm0,dm1] = interp_FRET_gradient(FRET_data,xdata2,Xdata2,Tdata2);
    
    
    %%%% Now fit to migration data. First initialize q and cost vectors
    q_all = cell(simnum,1);
    q0_all = cell(simnum,1);
    J_all = zeros(simnum,1);

   

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

    %load sparse matrices for later computation
    [A_pos,~,~,A_neg,~,~] = aMatrixConstruction(xn);
 
    options = optimset('maxiter',100);
       

    %each row of data matrix corresponds to data at given time point. In this
    %stat model, we only care about rows when t is a multiple of toCompare

    toCompare = 6;
    cell_data = cell_data(1:3*toCompare:end,:);

    tdata = tdata(1:3*toCompare:end);

    name_save = [stat_model '_' num2str(toCompare)];
        
    
    parfor k = 1:simnum

        
            if strcmp(IC_type,'interp') %don't estimate height
                %q = [v_1,v_2,..,v_8]^T;
                q0_all{k} = [.02*rand(4,1);.002*rand(4,1)];
                LB = [zeros(8,1)];
                UB = [inf(8,1)];
            else
                q0_all{k} = [.0075*rand(1,5),1];
                LB = [zeros(1,length(q0_all{k})-1),.7];
                UB = [inf*ones(1,length(q0_all{k})-1),1.2];
            end



            tic

            [q_all{k},J_all(k)] = fmincon(@(q) MLE_cost_D0_fewer_compare_space(cell_data...
                ,q,F1,F2,m0,m1,dm0,dm1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,...
                tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg),q0_all{k}...
                ,[],[],[],[],LB,UB,[],options);

            toc

                       
    end

    %save
    save(['/scratch/summit/jona8898/chem_fitting/FRET_interp_est_well_' ...
        num2str(well) '_' num2str(l) '_' name_save '_pred_'...
        num2str(pred_ind) '_space.mat' ],...
        'q_all','q0_all','J_all','F1','F2','m0','m1','dm0','dm1')

end
