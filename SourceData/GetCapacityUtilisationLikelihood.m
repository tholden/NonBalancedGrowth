function LogL = GetCapacityUtilisationLikelihood( CUForecast, Model, CUMat, Missing )
    
    CUMat( Missing, 1 ) = CUForecast;
    [ ~, ~, ~, LogL ] = infer( Model, CUMat( :, 1 ), 'X', CUMat( :, 2 : end ) );
    
end
