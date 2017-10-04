% plot_sim_v_data written 10-4-17 by JTN to simply take a parameter, q, and
% run the model against some data.

%define parameters for simulation:
    l = 3;
    well = 4;
    pred_ind = 3;
    stat_model = 'fewer';
    IC_type = 'interp';


%load parameter vector
    filename = ['FRET_interp_est_well_' num2str(well) '_'];
    filename2 = ['_fewer_6_pred_' num2str(pred_ind) '_space_final'];
    load([filename filename2 '.mat'])
    q = q_final{l};
    q(1) = q(1)/10000;
    q(3:4) = q(3:4)/4000;
    q(5:8) = zeros(1,4);
%     q(5:8) = q(5:8)/40;


%load in data
    load('ind_cell_prof_data.mat')
    
    %helps with FRET interpolation
    dx_large = 10;
    %choose data
    cell_data = avg_cell_data(ind_cell_data{well-1,2}(:,:,16:end),pred_ind)';
    FRET_data = avg_cell_data(ind_fret_data{well-1,2}(:,1:dx_large:end,16:end),pred_ind)';
    cell_data_std = ind_cell_data_sv{well-1,2}(:,16:end);
    
        
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    xdata2 = xdata(1:dx_large:end);
    tdata = 5:1/3:1/3*(tndata-1+5*3);
    %make ndgrid for interpolation
    [Xdata,Tdata] = ndgrid(xdata,tdata);
    [Xdata2,Tdata2] = ndgrid(xdata2,tdata);
    

%generate grids for computation, also helpful for interpolation
    %select which  grid size we're using
    xi = mod(l,4);
    if xi == 0
        xi = 4; %to determine grid size.
    end
    
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
 
    
   
%run simulation
    [model] = FRET_dep_convection_space(q,F1,F2,m0,m1,dm0,dm1,x,dx,xn,...
        x_int,xbd_0,xbd_1,t,dt,tn,tdata,xdata,IC,IC_type,BC_x_0,BC_x_1,A_pos,A_neg);
    
   
    
%plot simulation against data when t is a multiple of toCompare
    toCompare = 6;
    cell_data2 = cell_data(1:3*toCompare:end,:);
    tdata2 = tdata(1:3*toCompare:end);
    cell_data_std = cell_data_std(:,1:3*toCompare:end);

    [Xdata,Tdata] = meshgrid(xdata,tdata);
    [Xint,Tint] = meshgrid(xdata,tdata2);
    
    
    udata = interp2(Xdata,Tdata,model,Xint,Tint);


    figure
    hold on
    
    colors = distinguishable_colors(length(tdata2));
    
        
    for i = 1:length(tdata2)
        plot(xdata,udata(i,:),'color',colors(i,:))
        plot(xdata(1:3:end),cell_data2(i,1:3:end),'.','color',[colors(i,:)],...
            'markersize',15)
        fill([xdata(1:3:end) fliplr(xdata(1:3:end))],[(cell_data2(i,1:3:end) +...
            cell_data_std(1:3:end,i)') fliplr(cell_data2(i,1:3:end) -...
            cell_data_std(1:3:end,i)')],colors(i,:),'Facealpha',0.1,'edgecolor','none')
    end

    
%Plot video of simulation against speed
    figure
    n=4;
  
    
    %%%% 10-4-17, creating a chemokinesis function v(m) that satisfies v(0)=0.
    %then linear up to m0 . May consider v(m) = 0 for all m < m0 in the
    %future.
    msamp1 = augknt([0,m1,0,linspace(m0,m1,n)],2);
    %create spline functions
    v_spline1 = spmak(msamp1,[0 q(1:4)']);
% %     msamp1 = augknt([m0,m1,linspace(m0,m1,n)],2);
% %     %create spline functions
% %     v_spline1 = spmak(msamp1,q(1:4)');

    %%%%
    
    %for v(dm/dx)
    msamp2 = augknt([dm0,dm1,linspace(dm0,dm1,n)],2);
    %create spline functions
    v_spline2 = spmak(msamp2,q(5:end)');
    
    %interpolated FRET
    FRET = @(t) F1(x,t*ones(xn,1));
    dFRETdx = @(t) F2(x,t*ones(xn,1));
     
    FRET_speed = @(t) fnval(v_spline1,FRET(t));
    chem_speed = @(t) fnval(v_spline2,dFRETdx(t));

    % title_m = ['sim_' num2str(m) '_pred_' num2str(pred_ind) '.avi'];
    % f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    % vid = VideoWriter(title_m); %%title here
    % vid.Quality = 100;
    % vid.FrameRate = 5;
    % open(vid);

    for i = 1:3:length(tdata)
        subplot(2,1,1)
        yyaxis left
         hold off
        plot(xdata,model(i,:))
         hold on
         plot(xdata,cell_data(i,:),'k.')
        title(['t = ' num2str(tdata(i)) ' hours'])
        yyaxis right
        hold off
        plot(x,FRET_speed(tdata(i)))
        hold on
        plot(x,chem_speed(tdata(i)))

        subplot(2,1,2)
        hold off
        plot(x,FRET(tdata(i)))
        hold on
        plot(x,dFRETdx(tdata(i)))

        pause(1)

    %     writeVideo(vid, getframe(f1));
    end

    % close(vid)

