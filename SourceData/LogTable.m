function Table = LogTable( Table )

    VariableNames = Table.Properties.VariableNames;
    
    for i = 1 : length( VariableNames )
        
        Table.( VariableNames{ i } ) = log( Table.( VariableNames{ i } ) );
        
    end

end
