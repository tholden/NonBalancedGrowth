function x = lambertwM1( x )

    TooLow  = x <= - 1 / exp( 1 );
    TooHigh = x >= 0;
    Good    = ~( TooLow | TooHigh );
    
    x( Good )    = real( lambertw( -1, x( Good ) ) );
    x( TooLow )  = -1;
    x( TooHigh ) = -Inf;
    
end
