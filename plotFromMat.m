function isSuccess = plotFromMat(Filename)
% Load from file and plot
    isSuccess=false;
    try
        load(Filename);
        figure(1)
        loglog(sort(savedata(1,:)),sort(savedata(2,:)))
        hold on
        loglog(sort(savedata(1,:)),sort(savedata(3,:)))
        isSuccess=true;
    end
end