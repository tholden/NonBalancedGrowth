load SourceData

DataSheets = keys( SourceData );

HoursSheets = DataSheets( cellfun( @( x ) ~isempty( x ) && ( x == 1 ), strfind( DataSheets, 'T609' ) ) );

Hours = SourceData( HoursSheets{ 1 } );

Hours = Hours( :, 'HoursWorkedByFulltimeAndParttimeEmployees' );

HoursRowNames = Hours.Properties.RowNames;

for i = 2 : length( HoursSheets )

    NewHours = SourceData( HoursSheets{ i } );
    NewHours = NewHours( :, 'HoursWorkedByFulltimeAndParttimeEmployees' );
    
    NewHoursRowNames = NewHours.Properties.RowNames;
    
    Common = intersect( HoursRowNames, NewHoursRowNames );
    
    for j = 1 : length( Common )
        
       assert( all( table2array( Hours( Common{ j }, : ) ) == table2array( NewHours( Common{ j }, : ) ) ) );
       NewHours( Common{ j }, : ) = [];
       
    end
    
    Hours = [ Hours; NewHours ]; %#ok<AGROW>
    
    HoursRowNames = Hours.Properties.RowNames;
    
end

FAName = 'PrivateFixedAssets';

MainTable = SourceData( 'FAAt201-A' );
MainTable = AddYearsToTable( PrefaceVariableNames( MainTable( :, FAName ), 'Nominal' ) );

CurrentTable = SourceData( 'FAAt202-A' );
MainTable = outerjoin( MainTable, AddYearsToTable( PrefaceVariableNames( CurrentTable( :, FAName ), 'Real' ) ), 'MergeKeys', true );

MainTable = outerjoin( MainTable, AddYearsToTable( Hours ), 'MergeKeys', true );

GDPVars = { 'GrossDomesticProduct', 'PersonalConsumptionExpenditures', 'GrossPrivateDomesticInvestment', 'GovernmentConsumptionExpendituresAndGrossInvestment' };

CurrentTable = SourceData( 'T10103-A' );
MainTable = outerjoin( MainTable, AddYearsToTable( PrefaceVariableNames( CurrentTable( :, GDPVars ), 'Real' ) ), 'MergeKeys', true );

CurrentTable = SourceData( 'T10105-A' );
MainTable = outerjoin( MainTable, AddYearsToTable( PrefaceVariableNames( CurrentTable( :, GDPVars ), 'Nominal' ) ), 'MergeKeys', true );

MainTable = AddPricesToTable( MainTable, [ GDPVars, { FAName } ] );

CurrentTable = SourceData( 'GS1' );
MainTable = outerjoin( MainTable, AddYearsToTable( CurrentTable ), 'MergeKeys', true );

CurrentTable = SourceData( 'G17' );
MainTable = outerjoin( MainTable, AddYearsToTable( CurrentTable ), 'MergeKeys', true );

CurrentTable = SourceData( 'T20100-A' );
MainTable = outerjoin( MainTable, AddYearsToTable( CurrentTable( :, { 'CompensationOfEmployees', 'ProprietorsIncomeWithInventoryValuationAndCapitalConsumptionAdj', 'PopulationmidperiodThousands' } ) ), 'MergeKeys', true );
MainTable.Properties.VariableNames{ end } = 'Population';
MainTable.Population = MainTable.Population ./ 1e3;

CurrentTable = SourceData( 'T30100-A' );
MainTable = outerjoin( MainTable, AddYearsToTable( CurrentTable( :, 'CurrentTaxReceipts' ) ), 'MergeKeys', true );

MainTable.GrowthInRealGrossDomesticProduct = MainTable.RealGrossDomesticProduct ./ [ NaN; MainTable.RealGrossDomesticProduct( 1 : ( end - 1 ) ) ];

MainTable.LabourShare = MainTable.CompensationOfEmployees ./ ( MainTable.NominalGrossDomesticProduct - MainTable.ProprietorsIncomeWithInventoryValuationAndCapitalConsumptionAdj );

MainTable.PriceOfInvestmentOverPriceOfCapital = MainTable.PriceGrossPrivateDomesticInvestment ./ MainTable.( [ 'Price' FAName ] );

save MainTable MainTable

OldNames = { [ 'Real' FAName ], 'HoursWorkedByFulltimeAndParttimeEmployees', 'RealPersonalConsumptionExpenditures', 'RealGrossPrivateDomesticInvestment', 'RealGovernmentConsumptionExpendituresAndGrossInvestment', 'NominalGrossDomesticProduct', 'PricePersonalConsumptionExpenditures', 'PriceGrossPrivateDomesticInvestment', 'PriceGovernmentConsumptionExpendituresAndGrossInvestment', 'InterestRates', 'CapacityUtilisationTotalIndex', 'Population', 'CurrentTaxReceipts', 'GrowthInRealGrossDomesticProduct', 'LabourShare', 'PriceOfInvestmentOverPriceOfCapital' };
NewNames = { 'K'              , 'H'                                        , 'C'                                  , 'I'                                 , 'G'                                                      , 'NGDP'                       , 'PC'                                  , 'PI'                                 , 'PG'                                                      , 'R'            , 'U'                            , 'N'         , 'T'                 , 'G_GDP'                           , 'LS'         , 'Omega'                             };
NewScale = [ 1                , 1e-9                                       , 1                                    , 1                                   , 1                                                        , 1e-9                         , 1e-9                                  , 1e-9                                 , 1e-9                                                      , 1              , 1                              , 1e-9        , 1e-9                , 1                                 , 1            , 1                                     ];

LevelTable = ScaleTable( MainTable( :, OldNames ), NewScale );
LevelTable.Properties.VariableNames = NewNames;

save LevelTable LevelTable

FinalTable = LogTable( LevelTable );
FinalTable = PrefaceVariableNames( FinalTable, 'log_' );

save FinalTable FinalTable

SaveTableVariables( FinalTable, 'FinalData' );

EstimationPeriods = find( all( isfinite( table2array( FinalTable ) ), 2 ) );

syms kappa
assume( kappa > 0 & kappa < 1/2 );
Kt0OverU0 = exp( sum( log( LevelTable.K( EstimationPeriods - 1 ) + kappa * LevelTable.Omega( EstimationPeriods ) .* LevelTable.I( EstimationPeriods ) ) ) );
fprintf( '\nKt0OverU0 = ( %s ) ^ ( 1 / %d );\n\n', regexprep( strrep( func2str( matlabFunction( simplify( Kt0OverU0, 'Criterion', 'preferReal', 'Steps', 10 ) ) ), '@(kappa)', '' ), '\.(?=[\*\/\^])', '' ), length( EstimationPeriods ) );

DetrendedTable = DetrendTable( FinalTable, EstimationPeriods );

fprintf( '\n' );

save DetrendedTable DetrendedTable

SaveTableVariables( DetrendedTable, 'DetrendedData' );
