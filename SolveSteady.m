function Y_ = SolveSteady( hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ )

    if any( ~isfinite( [ hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ ] ) )
        error( 'Non-finite input.' );
    end

    Solution = SolveSteadyInternal( -10, 10, 21, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
    
    Y_ = exp( Solution );
    
end

function Solution = SolveSteadyInternal( LB, UB, NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ )

    if NPoints > 1000
        error( 'Failed to find the steady state.' );
    end
    
    StartPoints = linspace( LB, UB, NPoints );
        
    Resids = GetResids( StartPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );

    if ~any( isfinite( Resids ) )
        
        Solution = SolveSteadyInternal( LB - ( UB - LB ), UB + ( UB - LB ), 3 * NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
        return
        
    end
        
    DSignResids = diff( sign( Resids ) );
    Changes = find( ( DSignResids ~= 0 ) & ( ~isnan( DSignResids ) ) );
    
    if length( Changes ) > 1
        warning( 'Multiple solutions!' );
    end
    
    if isempty( Changes )
        [ ~, BestIndex ] = min( abs( Resids ) );
        if BestIndex == 1
            Solution = SolveSteadyInternal( LB - ( UB - LB ), StartPoints( 2 ), 2 * NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
        elseif BestIndex == length( StartPoints )
            Solution = SolveSteadyInternal( StartPoints( end - 1 ), UB + ( UB - LB ), 2 * NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
        elseif isnan( Resids( BestIndex + 1 ) ) && ( ~isnan( Resids( BestIndex - 1 ) ) )
            Solution = SolveSteadyInternal( StartPoints( BestIndex ), StartPoints( BestIndex + 1 ), NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
        elseif isnan( Resids( BestIndex - 1 ) ) && ( ~isnan( Resids( BestIndex + 1 ) ) )
            Solution = SolveSteadyInternal( StartPoints( BestIndex - 1 ), StartPoints( BestIndex ), NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
        else
            Solution = SolveSteadyInternal( StartPoints( BestIndex - 1 ), StartPoints( BestIndex + 1 ), NPoints, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );
        end
    else
        Change = Changes( 1 );
        Interval = StartPoints( [ Change, Change + 1 ] );
        Solution = fzero( @( In ) GetResids( In, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ ), ...
           Interval, optimset( 'Display', 'off' ) );
    end

end

function Out = GetResids( In, hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, AZ0, AL0, AC0, AG0, AK0, AH0, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ )

    Y_ = exp( In );

    t1 = (nu .* sY);
    t2 = delta_ .^ 2;
    t3 = (t2 .* Xi_);
    t4 = (kappa .* t3);
    t6 = (nu .* tauI_ .* xi_);
    t9 = (nu .* kappa);
    t10 = (xi_ .* t9);
    t14 = (Xi_ .* delta_);
    t15 = (kappa .* t14);
    t18 = (xi_ .* nu);
    t24 = (Xi_ .* kappa);
    t27 = (delta_ .* kappa);
    t35 = 6 .* t10 .* t14 - 2 .* t10 .* t3 - 4 .* t18 .* t14 + 4 .* t6 .* t14 - 6 .* t6 .* t15 - 4 .* t18 .* t24 + t18 .* t3 + 4 .* t6 .* t24 + 2 .* t6 .* t27 - t6 .* t3 + 2 .* t6 .* t4 + 4 .* t4;
    t36 = (Xi_ .* nu);
    t37 = (xi_ .* tauI_);
    t53 = -2 .* t18 .* t27 - 4 .* t37 .* t36 + 4 .* xi_ .* t36 - 4 .* t37 .* t9 - 4 .* Xi_ - 2 .* delta_ + 4 .* t10 + 6 .* t14 - 4 .* t15 + 4 .* t27 - 2 .* t3 + 4;
    t55 = 1 ./ Xi_;
    t61 = -2 + delta_;
    t66 = 1 ./ (2 .* t27 - delta_ + 2);
    t69 = exp(-(t66 ./ t61 ./ (-1 + tauI_) ./ xi_ .* t55 .* (t35 + t53)));
    t79 = lambertw(-1,-t55 .* t66 .* nu .* (2 .* t15 - t14 - 2 .* t24 + 2 .* Xi_ + 2 .* kappa) .* t69);
    t88 = log(t66 .* nu .* t69 ./ (-1 + tauK_) .* xi_ .* t61);
    t90 = t79 .* Xi_;
    t107 = log(t55 .* t66 ./ t79 .* (-(delta_ .* t90) - nu .* t14 - 0.2e1 .* nu .* t24 + (2 .* t9 .* t14) + (2 .* t27 .* t90) + (2 .* t36) + (2 .* t9) + (2 .* t90)));
    t108 = sY .* t107;
    t114 = (AK_ ./ AK0 ./ Kt0_) .^ sY;
    t120 = (aY .* PY_ ./ PI_ .* Omega_) .^ sY;
    t122 = log(t120 .* Y_ .* t114);
    t125 = sY - 0.1e1;
    t127 = 0.1e1 ./ nu;
    t128 = 0.1e1 ./ sY;
    t131 = exp(t128 .* t127 .* t125 .* (-nu .* t108 + nu .* t122 + (t79 .* t1) - t88 .* t1 + t108));
    t134 = Y_ .^ (t128 .* t125);
    t141 = (-0.1e1 ./ (-1 + aY) .* (-t131 .* aY + t134)) .^ (0.1e1 ./ t125 .* sY);
    H_ = H0_ ./ AH_ .* AH0 .* t141;
    U_ = exp(t127 .* t107);
        
    h_ = H_ ./ N_;
    
    PK_ = PI_ ./ Omega_;
    
    SR_ = (-xi_.*Xi_.*PK_.*(-1+tauI_).*(-2+delta_).*log(1-U_.^nu)+2.*Xi_.*(delta_-1).*PK_+2.*PK_)./((2.*(-1+tauK_)).*(-1+tauI_).*(Xi_.*((kappa-1./2).*delta_-kappa+1).*PK_+PK_.*kappa).*U_);
    I_ = ( aY .* PY_ ./ PK_ ) .^ sY .* ( Y_ ) .* ( AK_ ./ AK0 ./ Kt0_ ) .^ ( sY - 1 ) ./ ( U_ .* ( ( 1 - delta_ ./ 2 ) ./ delta_ .* Omega_ + kappa .* Omega_ ) ) ./ SR_ .^ sY;
    K_ = ( 1 - delta_ ./ 2 ) ./ delta_ .* Omega_ .* I_;
    % Kt_ = U_ .* ( K_ + kappa .* Omega_ .* I_ );
    G_ = tauGE_ .* ( 1 - tauE_ ) .* PY_ .* Y_ ./ PG_;
    E_ = tauGE_ .* tauE_ .* PY_ .* Y_ ./ PE_;
    C_ = 1 ./ PC_ .* ( PY_ .* Y_ - ( PI_ .* I_ + PG_ .* G_ + PE_ .* E_ - xi_ .* log( 1 - U_ .^ nu ) .* PK_ .* K_ ) );
    Z_ = ( aZ .* ( AC_ .* C_ .* N0 ./ AC0 ./ N_ ./ C0 ) .^ ( ( sZ - 1 ) ./ sZ ) + ( 1 - aZ ) .* ( AG_ .* G_ .* N0 ./ AG0 ./ N_ ./ G0 ) .^ ( ( sZ - 1 ) ./ sZ ) ) .^ ( sZ ./ ( sZ - 1 ) );
    X_ = ( aX .* ( AZ_ .* Z_ ./ AZ0 ) .^ ( ( sX - 1 ) ./ sX ) + ( 1 - aX ) .* ( AL_ .* ( hBar - h_ ) ./ AL0 ./ ( hBar - h0 ) ) .^ ( ( sX - 1 ) ./ sX ) ) .^ ( sX ./ ( sX - 1 ) );
    lambdaX_ = X_ .^ ( - 1 ./ sV );
    lambdaZ_ = aX .* lambdaX_ .* ( X_ ) .^ ( 1 ./ sX ) .* ( AZ_ .* Z_ ./ AZ0 ) .^ ( ( sX - 1 ) ./ sX ) ./ Z_;
    lambdaB_ = aZ .* ( 1 - tauC_ ) .* lambdaZ_ .* ( Z_ ) .^ ( 1 ./ sZ ) .* ( AC_ .* C_ .* N0 ./ AC0 ./ N_ ./ C0 ) .^ ( ( sZ - 1 ) ./ sZ ) .* N_ ./ C_;
        
    Out = ( 1 - aX ) ./ ( 1 - tauH_ ) .* lambdaX_ ./ lambdaB_ .* ( X_ ) .^ ( 1 ./ sX ) .* ( AL_ .* ( hBar - h_ ) ./ AL0 ./ ( hBar - h0 ) ) .^ ( ( sX - 1 ) ./ sX ) ./ ( hBar - h_ ) ./ ( ( 1 - aY ) .* PY_ ./ PC_ .* ( Y_ ) .^ ( 1 ./ sY ) .* ( AH_ .* H_ ./ AH0 ./ H0_ ) .^ ( ( sY - 1 ) ./ sY ) ./ H_ ) - 1;
    
    Out( imag( Out ) ~= 0 ) = NaN;
    
end
