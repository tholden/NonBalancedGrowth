function Table = DetrendTable( Table, EstimationPeriods )

    VariableNames = Table.Properties.VariableNames;
    
    T = size( Table, 1 );
    
    X = [ ones( T, 1 ), ( 1 : T ).' ];
    
    for i = 1 : length( VariableNames )
        
        VariableName = VariableNames{ i };
        
        Current = Table.( VariableName );
        
        Mean = mean( Current( EstimationPeriods ), 'includenan' );
        
        fprintf( '%s0 = %.64g;\n', VariableName( 5 : end ), exp( Mean ) );
        
        Beta = regress( Current( EstimationPeriods ), X( EstimationPeriods, : ) );
        
        Table.( VariableName ) = Current - X * Beta + Mean;
        
    end

end
