FileName = 'G17.xls';

URL = 'https://www.federalreserve.gov/datadownload/Output.aspx?rel=G17&series=3623066af5b1e36b7b505cd16d4c429f&lastobs=&from=01/01/1919&to=12/31/2118&filetype=spreadsheetml&label=include&layout=seriescolumn';

fprintf( '\nDownloading %s.\n', URL );

websave( FileName, URL, weboptions( 'Timeout', Inf ) );

fprintf( 'Detecting table size.\n' );

ImportOptions = detectImportOptions( FileName, 'Range', 'B8:IV65536' );

EndRow = regexp( ImportOptions.DataRange, '(?<=\:[A-Za-z]+)\d+', 'match', 'once' );
EndColumn = regexp( ImportOptions.DataRange, '(?<=\:)[A-Za-z]+', 'match', 'once' );

ImportOptions.RowNamesRange = [ 'A8:A' EndRow ];
ImportOptions.VariableNamesRange = [ 'B2:' EndColumn '2' ];

fprintf( 'Reading table.\n' );

CapacityUtilisation = readtable( FileName, ImportOptions, 'ReadVariableNames', true, 'ReadRowNames', true );

fprintf( 'Processing table.\n' );

CapacityUtilisation.Properties.VariableNames = replace( CapacityUtilisation.Properties.VariableNames, 'capacity', 'Capacity' );
CapacityUtilisation.Properties.VariableNames = replace( CapacityUtilisation.Properties.VariableNames, '_', '' );
CapacityUtilisation.Properties.VariableNames = replace( CapacityUtilisation.Properties.VariableNames, 'SaCAPUTL', '' );

Index1 = find( strcmp( CapacityUtilisation.Properties.VariableNames, 'TotalIndex' ), 1 );
Index2 = find( strcmp( CapacityUtilisation.Properties.VariableNames, 'ManufacturingSIC' ), 1 );
IndexOther = setdiff( 1 : size( CapacityUtilisation, 2 ), [ Index1 Index2 ] );
CapacityUtilisation = [ CapacityUtilisation( :, Index1 ), CapacityUtilisation( :, Index2 ), CapacityUtilisation( :, IndexOther ) ];

CUMat = table2array( CapacityUtilisation ) ./ 100;
CUMat = log( CUMat ./ ( 1 - CUMat ) );

CUMat = flip( CUMat );

CUMat = [ CUMat( :, 1 ), CUMat( :, 2 : end ), CUMat( :, 2 : end ) .* CUMat( :, 2 : end ) ];

Model = regARIMA( 'ARLags', [ 1 12 ], 'MALags', 1, 'Distribution', 't' );
Model = estimate( Model, CUMat( :, 1 ), 'X', CUMat( :, 2 : end ) );

Present = isfinite( CUMat( :, 1 ) );
Missing = ~Present;

CUForecast = forecast( Model, sum( Missing ), 'Y0', CUMat( Present, 1 ), 'X0', CUMat( Present, 2 : end ), 'XF', CUMat( Missing, 2 : end ) );

CUMat = flip( CUMat );
CUForecast = flip( CUForecast );
Missing = flip( Missing );

Model = regARIMA( 'ARLags', [ 1 12 ], 'MALags', 1, 'Distribution', 't' );
Model = estimate( Model, CUMat( :, 1 ), 'X', CUMat( :, 2 : end ) );

CUForecast = fminunc( @( CUF ) - GetCapacityUtilisationLikelihood( CUF, Model, CUMat, Missing ), CUForecast, optimoptions( @fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', Inf, 'MaxIterations', Inf ) );

CUMat( Missing, 1 ) = CUForecast;

CUMat = 1 ./ ( 1 + exp( - CUMat ) );

RowNames = CapacityUtilisation.Properties.RowNames;

RowNames = regexprep( RowNames, '\D.*', '' );

[ RowNames, ~, YearIndex ] = unique( RowNames );

Select = true( size( YearIndex ) );

if sum( YearIndex == 1 ) ~= 12   
    Select( YearIndex == 1 ) = false;
    RowNames = RowNames( 2 : end );
end

if sum( YearIndex == max( YearIndex ) ) ~= 12   
    Select( YearIndex == max( YearIndex ) ) = false;
    RowNames = RowNames( 1 : ( end - 1 ) );
end

CUVec = CUMat( Select, 1 );

CUVec = mean( reshape( CUVec, 12, length( RowNames ) ) );

try
    load SourceData;
catch
    SourceData = containers.Map;
end

delete( FileName );
[ ~, FileName ] = fileparts( FileName );
SourceData( FileName ) = table( CUVec(:), 'VariableNames', { 'CapacityUtilisationTotalIndex' }, 'RowNames', RowNames );

save SourceData SourceData;
