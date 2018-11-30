FileName = 'GS1.xls';

URL = 'https://fred.stlouisfed.org/graph/fredgraph.xls?drp=0&ts=12&tts=12&nt=0&thu=0&trc=0&id=GS1&cosd=1901-01-01&coed=2118-01-01&link_values=false&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Monthly&fam=avg&fgst=lin&line_index=1&transformation=lin&nd=1901-01-01';

fprintf( '\nDownloading %s.\n', URL );

websave( FileName, URL, weboptions( 'Timeout', Inf ) );

fprintf( 'Detecting table size.\n' );

ImportOptions = detectImportOptions( FileName, 'NumHeaderLines', 10 );

fprintf( 'Reading table.\n' );

InterestRates = readtable( FileName, ImportOptions, 'ReadVariableNames', true, 'ReadRowNames', true );

fprintf( 'Processing table.\n' );

RowNames = InterestRates.Properties.RowNames;

RowNames = regexprep( RowNames, '.*\D', '' );

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

RVec = table2array( InterestRates( :, 'GS1' ) );

RVec = RVec( Select );

RVec = mean( reshape( RVec, 12, length( RowNames ) ) );

try
    load SourceData;
catch
    SourceData = containers.Map;
end

delete( FileName );
[ ~, FileName ] = fileparts( FileName );
SourceData( FileName ) = table( 1 + RVec(:) ./ 100, 'VariableNames', { 'InterestRates' }, 'RowNames', RowNames );

save SourceData SourceData;
