function Out = TransposeTable( In )

    Data = num2cell( table2array( In ).', 1 );
    Out = table( Data{:}, 'VariableNames', matlab.lang.makeUniqueStrings( matlab.lang.makeValidName( In.Properties.RowNames, 'ReplacementStyle', 'delete' ) ), 'RowNames', cellfun( @FixRowName, In.Properties.VariableNames, 'UniformOutput', false ) );

end

function Name = FixRowName( Name )

    if length( Name ) >= 2 && Name( 1 ) == 'x' && ~isletter( Name( 2 ) )
        
        Name = Name( 2 : end );
        
    end

end
