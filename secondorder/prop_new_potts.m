function Z=prop_new_potts(Z_old,G_start,chess)
N=size(Z_old);
Z=Z_old+randi(G_start-1,N);
chess2=G_start*(Z==G_start);
Z=chess.*(mod(Z,G_start))+chess2;
end


            
        