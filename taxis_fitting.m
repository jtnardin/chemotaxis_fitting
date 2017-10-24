%taxis_fitting.m written  10-24-17 by JTN to perform
%OLS optimization for the equation u_t + v(t)u_x =0 in a parallelized
%for loop to mean experimental data

%in this script, v(t) = v_1(m) + v_2(dm/dx) to look into chemokinesis and
%chemotaxis, where m(t,x) is interpolated FRET ratio data

%last updated by JTN : 10-4-17 to use csaps for FRET interpolation. Less
%accurate to data, but smoother and less sensitive to the jaggedness.

%updated by JTN : 10-24-17 to include diffusion, only consider every 5th
%data point

function taxis_fitting(l,m,pred_ind,sim_type,simnum)

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
    
    D = dt/dx^2*sparse(repmat(x_int,1,3),[x_int-1 x_int x_int+1],[ones(1,xn-2) ...
    -2*ones(1,xn-2) ones(1,xn-2)],xn,xn);
    
 
    options = optimset('maxiter',100,'display','iter');
       

    %each row of data matrix corresponds to data at given time point. In this
    %stat model, we only care about rows when t is a multiple of toCompare

    %%% also only looking at every 5th point in space now (10-25)
    
    toCompare = 6;
    
    cell_data = cell_data(1:3*toCompare:end,1:5:end);

    tdata = tdata(1:3*toCompare:end);
    xdata = xdata(1:5:end);

    name_save = [sim_type '_' num2str(toCompare)];
        
    
    for k = 1:simnum

        
         %q = [v_1,v_2,..,v_8]^T;
         q0_all{k} = [6e-3:1e-3:9e-3 1e-4];
         LB = [zeros(5,1)];
         UB = [inf(5,1)];
         

        tic

        [q_all{k},J_all(k)] = fmincon(@(q) MLE_cost_D0_fewer_compare_space(cell_data...
            ,q,F1,F2,m0,m1,dm0,dm1,x,dx,xn,x_int,xbd_0,xbd_1,t,dt,tn,...
            tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg,D,sim_type),q0_all{k}...
            ,[],[],[],[],LB,UB,[],options);

        toc

                       
    end

    %save
    save(['/scratch/summit/jona8898/chem_fitting/FRET_fitting_well_' ...
        num2str(well) '_xn_' num2str(l) '_' name_save '_pred_'...
        num2str(pred_ind) '.mat' ],...
        'q_all','q0_all','J_all','F1','F2','m0','m1','dm0','dm1')

end
