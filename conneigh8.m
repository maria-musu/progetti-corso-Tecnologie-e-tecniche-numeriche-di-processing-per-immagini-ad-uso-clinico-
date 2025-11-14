function neighbors8=conneigh8(Nrows,Ncols,ind)
    % in ingresso diamo tutta la matrice tranne una 
    % striscia di bordi
    [x,y] = ind2sub([Nrows Ncols],ind);

%             if x==1 || y==1 || x==Ncols || y==Nrows
            rowSub=[y-1 y y+1 y-1 y+1 y-1 y y+1];
            colSub=[x-1 x-1 x-1 x x x+1 x+1 x+1];
                
                    
%             end

    neighbors8=sub2ind([Nrows Ncols],colSub,rowSub);

end