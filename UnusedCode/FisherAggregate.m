function QSum = FisherAggregate( Q, P )

    QL = Q( 1 : ( end - 1 ), : );
    PL = P( 1 : ( end - 1 ), : );
    
    QC = Q( 2 : end, : );
    PC = P( 2 : end, : );
    
    QCPC = sum( QC .* PC, 2 );
    QCPL = sum( QC .* PL, 2 );
    QLPC = sum( QL .* PC, 2 );
    QLPL = sum( QL .* PL, 2 );
    
    QSum = cumprod( [ 1; sqrt( ( QCPC ./ QLPC ) .* ( QCPL ./ QLPL ) ) ] );

end
