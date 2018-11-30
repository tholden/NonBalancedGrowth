function delta_ = SolveSteadySimple( kappa, nu, aY, sY, Xi_, Upsilon_, PC_, PI_, PG_, PK_, C_, G_, tauK_, tauI_, tauE_, U_, Kt_ )

    StartPoints = -10 : 0.1 : 0;
    
    Resids = GetResids( StartPoints, kappa, nu, aY, sY, Xi_, Upsilon_, PC_, PI_, PG_, PK_, C_, G_, tauK_, tauI_, tauE_, U_, Kt_ );
    
    if any( isfinite( Resids ) )
        DSignResids = diff( sign( Resids ) );
        Change = find( ( real( DSignResids ) ~= 0 ) & ( imag( DSignResids ) == 0 ), 1 );
        Interval = StartPoints( [ Change, Change + 1 ] );
        
        try
            Solution = fzero( @( log_delta_ ) GetResids( log_delta_, kappa, nu, aY, sY, Xi_, Upsilon_, PC_, PI_, PG_, PK_, C_, G_, tauK_, tauI_, tauE_, U_, Kt_ ), ...
               Interval, optimset( 'Display', 'off' ) );
        catch
        end
    else
    	Solution = -3;
    end
    
    delta_ = exp( Solution );

end

function resid = GetResids( log_delta_, kappa, nu, aY, sY, Xi_, Upsilon_, PC_, PI_, PG_, PK_, C_, G_, tauK_, tauI_, tauE_, U_, Kt_ )

    delta_ = exp( log_delta_ );
    
    t4 = Xi_ .* delta_ - Xi_ + 1;
    t6 = (U_ .^ nu);
    t7 = -1 + t6;
    t9 = (kappa - 0.1e1 ./ 0.2e1) .* delta_;
    t1_ = 0.1e1 + t9;
    t16 = log(-t7);
    t25 = 0.1e1 ./ (-Xi_ .* t7 .* t1_ .* t16 + t6 .* nu .* (Xi_ .* (t9 - kappa + 0.1e1) + kappa));
    t27 = 1 ./ (-1 + tauI_);
    t28 = t25 .* t27;
    t33 = tauE_ - 1;
    t43 = (PG_ .* G_);
    t54 = Kt_ .* PI_;
    t88 = t4 .* t7 .* t1_;
    t9_ = 1 ./ (-2 + delta_);
    t95 = 0.1e1 ./ U_;
    t97 = 0.1e1 - delta_ ./ 0.2e1;
    t98 = 1 ./ delta_;
    t104 = t95 ./ (t97 .* t98 .* Upsilon_ + kappa .* Upsilon_);
    t124 = (aY .* ((PC_ .* C_) + t54 .* t104 + t43 + 0.2e1 .* t88 .* t28 .* t9_ .* t16 .* PK_ .* t97 .* t98 .* Upsilon_ .* Kt_ .* t104 - (t43 .* tauE_ ./ t33)) ./ PK_) .^ sY;
    t125 = 0.1e1 ./ Kt_;
    t127 = t125 .^ (sY - 1);
    t139 = ((1 + Xi_ .* (delta_ - 1)) .* nu .* t6 ./ (-1 + tauK_) .* t95 .* t28) .^ sY;
    resid = log(t124 .* t127 ./ t139 .* t125);
    
end
