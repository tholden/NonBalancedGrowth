function SourceData = ObtainSingleBEASpreadsheet( URL, SheetNamePrefix, SheetNameSuffix, Tables )

    FileName = 'BEA.xlsx';

    fprintf( '\nDownloading %s.\n', URL );
    
    websave( FileName, URL, weboptions( 'Timeout', Inf ) );
    
    SourceData = containers.Map;
    
    for i = 1 : length( Tables )
    
        SheetName = [ SheetNamePrefix Tables{i} SheetNameSuffix ];
        
        fprintf( '\nWorking on sheet %s.\n', SheetName );
        
        fprintf( 'Reading scale.\n' );
        
        ScaleText = table2array( readtable( FileName, 'Sheet', SheetName, 'Range', 'A2:A2', 'ReadVariableNames', false, 'ReadRowNames', false ) );
        
        if iscell( ScaleText )
            ScaleText = ScaleText{ 1 };
        end
        
        ScaleText = regexp( regexprep( ScaleText, '\W', ' ' ), '\w+', 'match', 'once' );
        
        switch lower( ScaleText )
            case 'percent'
                Scale = 1e-2;
            case 'index'
                Scale = 1;
            case 'tens'
                Scale = 1e1;
            case 'hundreds'
                Scale = 1e2;
            case 'thousands'
                Scale = 1e3;
            case 'millions'
                Scale = 1e6;
            case 'billions'
                Scale = 1e9;
            case 'trillions'
                Scale = 1e12;
            case ''
                Scale = 1;
                warning( 'No units detected in sheet "%s".', SheetName );
            otherwise
                Scale = 1;
                warning( 'Unrecognised unit scale "%s" in sheet "%s".', ScaleText, SheetName );
        end
        
        fprintf( 'Detecting table size.\n' );
        
        ImportOptions = detectImportOptions( FileName, 'Sheet', SheetName, 'Range', 'D9:IV65536' );
        
        EndRow = regexp( ImportOptions.DataRange, '(?<=\:[A-Za-z]+)\d+', 'match', 'once' );
        EndColumn = regexp( ImportOptions.DataRange, '(?<=\:)[A-Za-z]+', 'match', 'once' );

        ImportOptions.RowNamesRange = [ 'B9:B' EndRow ];
        ImportOptions.VariableNamesRange = [ 'D8:' EndColumn '8' ];

        fprintf( 'Reading table.\n' );
        
        Table = readtable( FileName, ImportOptions, 'ReadVariableNames', true, 'ReadRowNames', true );
        
        fprintf( 'Processing table.\n' );

        SourceData( SheetName ) = ScaleTable( TransposeTable( Table ), Scale );
        
    end
    
    delete( FileName );

end
