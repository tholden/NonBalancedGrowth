function Table = AddYearsToTable( Table )

    RowNames = str2double( Table.Properties.RowNames );

    Table.Year = RowNames;
    
    Table = [ Table( :, end ) Table( :, 1 : ( end - 1 ) ) ];

end
