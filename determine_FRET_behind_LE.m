% determine_FRET_behind_LE.m written 10-4-17 by JTN to determine the
% min and max FRET levels behind the Leading edge 

function [m0,m1] = determine_FRET_behind_LE(xdata,xdata2,cell_data,FRET_data,tdata)

    %initialize
    FRET_val_back = zeros(size(FRET_data));
    
    %figure out min, max FRET values behind the leading edge.
    %this is done at each time step by calculating the leading edge, then
    %only consider x values behind the leading edge, all other set to zero
    for i = 1:length(tdata)
        FRET_val_back(i,:) = (xdata2<=leading_edge_calc(cell_data(i,:),xdata,0.2,0)).*FRET_data(i,:);
    end
    
    
    %now take the min and max over nonzero values
    m0 = min(FRET_val_back(FRET_val_back~=0));
    m1 = max(FRET_val_back(FRET_val_back~=0));
    


end