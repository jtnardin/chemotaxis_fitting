function [F1,F2,dm0,dm1] =  interp_FRET_gradient(FRET_data,xdata2,Xdata2,Tdata2)


    %initialize
    FRET_data_interp = zeros(size(FRET_data));
    FRET_grad_data_interp = zeros(size(FRET_data));
    
    for i = 1:size(FRET_data_interp,1)
        %interpolate FRET data ... very smooth, not as accurate using csaps
        pp = csaps(xdata2,FRET_data(i,:),.999999);
        %create grid to use for time-dep sampling
        FRET_data_interp(i,:) = fnval(pp,xdata2);
        %do same for smooth gradient
        FRET_grad_data_interp(i,:) = fnval(fnder(pp),xdata2);
    end
    
    
    %keep interpolated FRET level positive
    FRET_data_interp = max(FRET_data_interp,0);
    
    
    %create time-dep samplers
    F1 = griddedInterpolant(Xdata2,Tdata2,FRET_data_interp');
    F2 = griddedInterpolant(Xdata2,Tdata2,FRET_grad_data_interp');
    
    %find min,max values
    dm0 = min(min(min(FRET_grad_data_interp)));
    dm1 = max(max(max(FRET_grad_data_interp)));
    

end