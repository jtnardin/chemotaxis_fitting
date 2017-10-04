%for bookkeeping
m = 5;
pred_ind = 3;

dx_large = 10;

%load in data
load('ind_cell_prof_data.mat')

%choose data
cell_data = avg_cell_data(ind_cell_data{m-1,2}(:,:,16:end),pred_ind)';
% FRET_data = avg_cell_data(ind_fret_data{m-1,2}(:,:,16:end),pred_ind)';
FRET_data = squeeze(ind_fret_data{m-1,2}(pred_ind,1:dx_large:end,16:end))';

%initialize data grids
[tndata,xndata] = size(cell_data);
xdata = linspace(0,1,xndata);
xdata2 = xdata(1:dx_large:end);
tdata = 5:1/3:1/3*(tndata-1+5*3);

[Xdata,Tdata] = ndgrid(xdata,tdata);
[Xdata2,Tdata] = ndgrid(xdata2,tdata);

x = linspace(0,1,100);
x=x';
t = linspace(tdata(1),tdata(end),400);
[X,T] = ndgrid(x,t);

for i = 1:size(FRET_data,1)
    pp = csaps(xdata2,FRET_data(i,:),.9999);
    yyaxis left
    hold off
    fnplt(pp)
    hold on
    plot(xdata2,FRET_data(i,:),'k.')
    yyaxis right
    hold off
    fnplt(fnder(pp))
    axis([0 1 -4e4 1e4])
    pause(.25)
end


