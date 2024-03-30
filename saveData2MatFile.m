function isSuccess = saveData2MatFile(fileTitle,data, loc)
% Save the cvx data
isSuccess=false;
    try
        load(fileTitle);
        shape=size(savedata);
        if loc == -1      
            if shape==[0,0]
                savedata=data;
            else
                savedata(:,shape(2)+1)=data;          
            end
        else
            savedata(:,loc)=data;
        end
        save(fileTitle,'savedata');
        isSuccess=true;
    end
end