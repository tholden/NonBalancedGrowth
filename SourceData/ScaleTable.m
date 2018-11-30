function Table = ScaleTable( Table, Scale )

    VariableNames = Table.Properties.VariableNames;
    
    Scale = Scale(:) .* ones( length( VariableNames ), 1 );
    
    for i = 1 : length( VariableNames )
        
        Table.( VariableNames{ i } ) = Table.( VariableNames{ i } ) .* Scale( i );
        
    end

end
