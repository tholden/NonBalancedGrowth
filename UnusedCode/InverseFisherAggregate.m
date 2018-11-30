function QE = InverseFisherAggregate( QSum, Q, P, QPE )

    GQSum = QSum( 2 : end ) ./ QSum( 1 : ( end - 1 ) );

    QL = Q( 1 : ( end - 1 ), : );
    PL = P( 1 : ( end - 1 ), : );
    QPEL = QPE( 1 : ( end - 1 ), : );
    
    QC = Q( 2 : end, : );
    PC = P( 2 : end, : );
    QPEC = P( 2 : end, : );
    
    QCPC = sum( QC .* PC, 2 ) + QPEC;
    QCPLb = sum( QC .* PL, 2 );
    % QCPL = QCPLb + QEC .* QPEL ./ QEL;
    QLPCb = sum( QL .* PC, 2 );
    % QLPC = QLPCb + QPEL .* QPEC ./ QEC;
    QLPL = sum( QL .* PL, 2 ) + QPEL;
    
    T1 = QLPL ./ QCPC .* GQSum .* GQSum;
    T2 = QCPLb - T1 .* QLPCb;
    
    QELt = 1;
    
    QE = zeros( length( GQSum ), 1 );
    
    for t = 1 : length( GQSum )
    
        a = QPEL( t ) / QELt;
        b = T2( t );
        c = - T1( t ) * QPEL( t ) * QPEC( t );
        
        QECtP = 0.5 * ( -b + sqrt( b * b - 4 * a * c ) ) / a;
        QECtN = 0.5 * ( -b - sqrt( b * b - 4 * a * c ) ) / a;
        
        if QPE( t ) >= 0
            if QECtP >= 0 && QECtN < 0
                QECt = QECtP;
            elseif QECtP < 0 && QECtN >= 0
                QECt = QECtN;
            else
                error( 'No or multiple solutions with correct sign.' );
            end
        else
            if QECtP >= 0 && QECtN < 0
                QECt = QECtN;
            elseif QECtP < 0 && QECtN >= 0
                QECt = QECtP;
            else
                error( 'No or multiple solutions with correct sign.' );
            end
        end
        
        QE( t ) = QECt;
        
        QELt = QECt;
    
    end
    
    QE = [ 1; QE ];

end
