function isSuccess = plotFromMat(Filename)
% Load from file and plot
    isSuccess=false;
    try
        load(Filename);
        loglog(sort(savedata(1,:)),sort(savedata(2,:)))
        isSuccess=true;
    end
end