load MainTable

figure( 1 );

subplot( 1, 4, 1 );

plot( MainTable.Year, log( MainTable.PriceGrossPrivateDomesticInvestment ./ MainTable.PricePersonalConsumptionExpenditures ) );

axis square;

title( { 'Logarithm of the price of', 'investment in units of consumption goods.' } );

subplot( 1, 4, 2 );

plot( MainTable.Year, log( 1 ./ MainTable.PriceOfInvestmentOverPriceOfCapital ) );

axis square;

title( { 'Logarithm of the price of', 'capital in units of investment goods.' } );

subplot( 1, 4, 3 );

plot( MainTable.Year, log( MainTable.PriceGrossPrivateDomesticInvestment ./ MainTable.PricePersonalConsumptionExpenditures ./ MainTable.PriceOfInvestmentOverPriceOfCapital ) );

axis square;

title( { 'Logarithm of the price of', 'capital in units of consumption goods.' } );

subplot( 1, 4, 4 );

plot( MainTable.Year, log( MainTable.RealPrivateFixedAssets ./ MainTable.RealGrossDomesticProduct ) );

axis square;

title( 'Logarithm of real capital over real GDP.' );

set( gcf, 'Position', [ 0 0 1920 1200 ] );

savefig RelativePrices
saveas( gcf, 'RelativePrices.emf' );

