@#define ValueFunction          = 0
@#define RandomParamInit        = 1
@#define IRFLength              = 0
@#define SimulationLength       = 50
@#define Estimation             = 0
@#define CMAES                  = 1
@#define Detrend                = 0
@#define LoadCurrentBest        = 1
@#define LoadInitParams         = 0
@#define UseDLL                 = 1
@#define GrowthIterations       = 1
@#define MeasurementError       = 0

@#if Detrend
    @#define GeneralSteadyState = 0
@#else
    @#define GeneralSteadyState = 1
    var GrowthSwitch;
@#endif

parameters N0 Omega0 h0 C0 G0 U0;

N0        = 0.23917944827248993622248462997958995401859283447265625;
Omega0    = 0.1001466637724954811261568465852178633213043212890625;

h0        = 173.26879279666428601558436639606952667236328125 / N0;
C0        = 45.3949137688586432659576530568301677703857421875;
G0        = 62.9526721260935033797068172134459018707275390625;
U0        = 0.80977743586948125464886061308789066970348358154296875;

@#define Kt0OverU0 = "( (kappa*2.131498647350859+2.4468e1)*(kappa*2.673962398868458+3.2258e1)*(kappa*5.979182466439458+3.41e2/5.0)*(kappa*7.756236898312263+1.04706e2)*(kappa*3.160212583655161+3.9674e1)*(kappa*2.544845683940409+2.6374e1)*(kappa*7.145087223483435+9.8858e1)*(kappa*4.50959986066657+6.1493e1)*(kappa*3.397589585666294+7.19e2/2.0e1)*(kappa*1.993307489919355+2.3669e1)*(kappa*3.161770624518119+3.4546e1)*(kappa*5.520698134898874+9.7193e1)*(kappa*7.624412223101486+7.7559e1)*(kappa*4.77162010372731+6.0109e1)*(kappa*6.485020521462299+9.8058e1)*(kappa*2.723003784245451+2.7538e1)*(kappa*7.406129222524263+1.0e2)*(kappa*1.338986494111875+1.8174e1)*(kappa*7.488588083186454+9.3527e1)*(kappa*2.601750763778565+2.8745e1)*(kappa*3.530883307641953+4.076e1)*(kappa*1.591901541282585+2.0298e1)*(kappa*2.697002654734859+2.9839e1)*(kappa*7.657976061879295+9.1034e1)*(kappa*2.251022907153729+2.5357e1)*(kappa*6.506810907636585+7.0198e1)*(kappa*6.883867675272215+9.5686e1)*(kappa*6.827049062024195+8.2652e1)*(kappa*3.872954473757136+4.216e1)*(kappa*5.618079207177916+6.6481e1)*(kappa*7.525436955720581+8.8769e1)*(kappa*4.687534228148148+5.0439e1)*(kappa*4.602798607119817+5.3807e1)*(kappa*3.924605891292202+4.3823e1)*(kappa*6.943592272729071+8.4548e1)*(kappa*6.140209455558063+9.7562e1)*(kappa*3.450979207665422+4.8283e1)*(kappa*1.788618346534654+2.2267e1)*(kappa*5.480295808637803+6.4948e1)*(kappa*1.622538211071107+1.8797e1)*(kappa*7.688403142005601+1.0142e2)*(kappa*8.074384192904029+1.03053e2)*(kappa*3.739216987069166+4.9279e1)*(kappa*7.286878347290541+7.4882e1)*(kappa*6.906911288963955+7.2395e1)*(kappa*4.913411889533934+5.8583e1)*(kappa*7.022442023668898+8.0353e1)*(kappa*1.769401375568129+4.31e2/2.0e1)*(kappa*4.720446396991996+6.2514e1)*(kappa*5.00669557375679+6.3577e1)*(kappa*7.28137654683528+8.6563e1)*(kappa*3.494789104523742+4.5558e1)*(kappa*2.890076041994304+3.3352e1)*(kappa*1.478634284332689+2.0993e1)*(kappa*1.785776138027187+2.2979e1)*(kappa*4.640916680528123+5.2086e1)*(kappa*3.069311968976155+3.751e1)*(kappa*2.852291505425646+3.1017e1)*(kappa*4.691298120348596+5.5455e1)*(kappa*4.747909196429344+5.7018e1)*(kappa*8.016105455973176+1.06298e2)*(kappa*2.695503782484682+3.8777e1)*(kappa*1.645122723004695+1.9553e1)*(kappa*3.87741305128303+4.6927e1) ) ^ ( 1 / 64 )"

@#includepath "DynareTransformationEngine"

@#include "Initialize.mod"

@#if RandomParamInit
    @#define DefaultSteady = "( 0.02 * ( rand - 0.5 ) )"
@#else
    @#define DefaultSteady = "0"
@#endif

@#define EndoVariables = EndoVariables + [ "Y", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "C", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "I", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "G", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "E", "-Inf", "Inf" ]

@#define EndoVariables = EndoVariables + [ "K", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "R", "0", "Inf" ]

@#define EndoVariables = EndoVariables + [ "ht", "0", "1" ]
@#define EndoVariables = EndoVariables + [ "Ut", "0", "1" ]

@#define EndoVariables = EndoVariables + [ "H", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "U", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "LS", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "T", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "G_GDP", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "NGDP", "0", "Inf" ]

@#define EndoVariables = EndoVariables + [ "lambdaK", "0", "Inf" ]

@#if ValueFunction

    @#define EndoVariables = EndoVariables + [ "V", "0", "Inf" ]
    @#define EndoVariables = EndoVariables + [ "Ft", "0", "Inf" ]
    
@#endif

@#define Taxes = [ "GE", "E", "H", "C", "K", "I", "Y" ]

@#for Tax in Taxes

    @#define TaxVar = "tau" + Tax
    @#define EndoVariableName = TaxVar + "t"
    @#define ShockProcessName = "S" + EndoVariableName
    
    @#define EndoVariables = EndoVariables + [ EndoVariableName, "-Inf", "Inf" ]
    @#define ShockProcesses = ShockProcesses + [ ShockProcessName, "-Inf", "Inf", "0", "rho_" + ShockProcessName, "sigma_" + ShockProcessName ]
    
    parameters @{EndoVariableName}_STEADY phi_@{Tax} psi_@{Tax} chi_@{Tax} gamma_@{EndoVariableName} iota0_@{EndoVariableName} iota1_@{EndoVariableName};
    
    @{EndoVariableName}_STEADY = ( 1 + @{DefaultSteady} ) / 2;
    phi_@{Tax} = 0.001;
    psi_@{Tax} = 0.001;
    chi_@{Tax} = 0.001;
    gamma_@{EndoVariableName} = 0.001;
    iota0_@{EndoVariableName} = 0.001;
    iota1_@{EndoVariableName} = 0.001;

    parameters rho_@{ShockProcessName} sigma_@{ShockProcessName};
    
    rho_@{ShockProcessName} = 0.9;
    sigma_@{ShockProcessName} = 0.001;

@#endfor

tauGEt_STEADY = -log( 2 / ( 0.213745098247544 + 1 ) - 1 );

@#define LogI1ShockProcesses = [ "N", "deltat", "xi", "thetaK", "thetaI", "Omega", "AX", "AF", "AZ", "AL", "AC", "AG", "AK", "AH", "PY" ]

@#define PriceVars = [ "C", "I", "G", "E" ]

@#define PriceLogI1ShockProcesses = [ "" ] - [ "" ]

@#for PriceVar in PriceVars

    @#define PriceVarName = "P" + PriceVar
    @#define LogI1ShockProcess = PriceVarName + "t"

    @#define EndoVariables = EndoVariables + [ PriceVarName, "0", "Inf" ]
    @#define PriceLogI1ShockProcesses = PriceLogI1ShockProcesses + [ LogI1ShockProcess ]
    
    parameters @{PriceVarName}_STEADY;

@#endfor

@#define LogI1ShockProcesses = LogI1ShockProcesses + PriceLogI1ShockProcesses

PC_STEADY = 48.0136444414831515814512385986745357513427734375;
PI_STEADY = 14.6635323596550488645107179763726890087127685546875;
PG_STEADY = 11.316876658901950492008836590684950351715087890625;
PE_STEADY = 11.316876658901950492008836590684950351715087890625;

@#for LogI1ShockProcess in LogI1ShockProcesses

    @#define ShockProcessName = "S" + LogI1ShockProcess

    @#define EndoVariables = EndoVariables + [ LogI1ShockProcess, "0", "Inf" ]
    @#define ShockProcesses = ShockProcesses + [ ShockProcessName, "0", "Inf", "1", "rho_" + ShockProcessName, "sigma_" + ShockProcessName ]
    
    parameters @{LogI1ShockProcess}_STEADY G@{LogI1ShockProcess}_STEADY gamma_@{LogI1ShockProcess} iota0_@{LogI1ShockProcess} iota1_@{LogI1ShockProcess};
    
    @{LogI1ShockProcess}_STEADY = exp( @{DefaultSteady} );
    G@{LogI1ShockProcess}_STEADY = exp( @{DefaultSteady} );
    gamma_@{LogI1ShockProcess} = 0.001;
    iota0_@{LogI1ShockProcess} = 0.001;
    iota1_@{LogI1ShockProcess} = 0.001;

    parameters rho_@{ShockProcessName} sigma_@{ShockProcessName};
    
    rho_@{ShockProcessName} = 0.9;
    sigma_@{ShockProcessName} = 0.002;

@#endfor

N_STEADY     = N0;
Omega_STEADY = Omega0;

deltat_STEADY  = 1 / ( 1 / 0.1 - 1 );
xi_STEADY      = 0.0703741487211254;
thetaK_STEADY  = 0.1;
thetaI_STEADY  = 0.05;
PY_STEADY      = 4425.36066522422;

@#define FixedVarianceLogAR1Shocks = [ "A0", "A1" ]

@#for FixedVarianceLogAR1Shock in FixedVarianceLogAR1Shocks

    @#define ShockProcesses = ShockProcesses + [ FixedVarianceLogAR1Shock, "0", "Inf", "1", "rho_" + FixedVarianceLogAR1Shock, "1" ]

    parameters rho_@{FixedVarianceLogAR1Shock};
    
    rho_@{FixedVarianceLogAR1Shock} = 0.9;
    
@#endfor

@#define FixedVarianceAR1Shocks = [ "tau0", "tau1" ]

@#for FixedVarianceAR1Shock in FixedVarianceAR1Shocks

    @#define ShockProcesses = ShockProcesses + [ FixedVarianceAR1Shock, "-Inf", "Inf", "0", "rho_" + FixedVarianceAR1Shock, "1" ]

    parameters rho_@{FixedVarianceAR1Shock};
    
    rho_@{FixedVarianceAR1Shock} = 0.9;
    
@#endfor

@#define FixedVarianceShocks = FixedVarianceLogAR1Shocks + FixedVarianceAR1Shocks

@#define CESVars = [ "V", "X", "Z", "Y" ]

@#for CESVar in CESVars
    
    parameters a@{CESVar} s@{CESVar};
    
    a@{CESVar} = 0.5  * ( 1 + @{DefaultSteady} );
    s@{CESVar} = 0.99 * ( 1 + @{DefaultSteady} );
    
@#endfor

aV = 0.04 * ( 1 + @{DefaultSteady} );

parameters eta hBar kappa nu beta;

eta   = 0.1              * ( 1 + @{DefaultSteady} );
hBar  = 3665.43821935251 * ( 1 + @{DefaultSteady} );
kappa = 0.4              * ( 1 + @{DefaultSteady} );
nu    = 2                * ( 1 + @{DefaultSteady} );

beta  = 1 - aV;

@#define CESProductivityVars = [ "AX", "AF", "AZ", "AL", "AC", "AG", "AK", "AH" ]

@#define EstimationFixedSteady = [ "N", "deltat", "xi", "Omega", "GPY" ] + CESProductivityVars + PriceLogI1ShockProcesses

@#define ThetaVariables = [ "thetaK", "thetaI" ]

@#include "CreateShocks.mod"
@#include "ClassifyDeclare.mod"

parameters Initial_level_tauGEt Initial_level_tauEt Initial_level_tauHt Initial_level_tauCt Initial_level_tauKt Initial_level_tauIt Initial_level_tauYt Initial_log_PCt Initial_log_PIt Initial_log_PGt Initial_log_PEt Initial_log_N Initial_log_deltat Initial_log_xi Initial_log_thetaK Initial_log_thetaI Initial_log_Omega Initial_log_AX Initial_log_AF Initial_log_AH Initial_log_AZ Initial_log_AL Initial_log_AC Initial_log_AG Initial_log_AK Initial_log_PY;

Initial_level_tauGEt = 0.3636409520570230147207269055798;
Initial_level_tauEt = 0.067676425207796431227080802273122;
Initial_level_tauHt = 2.7882543706137341565920451103011;
Initial_level_tauCt = 0.94633636823736699827946949881152;
Initial_level_tauKt = 0.80519050529227786938690769602545;
Initial_level_tauIt = -3.6996936395881019343789830600144;
Initial_level_tauYt = 0;
Initial_log_PCt = -4.2636116095333367326247753226198;
Initial_log_PIt = -5.4497332195831296175470015441533;
Initial_log_PGt = -5.7088017180654677673601327114739;
Initial_log_PEt = -5.7088017180654677673601327114739;
Initial_log_N = -1.430541179190143896704512371798;
Initial_log_deltat = -log( 1 / exp( -3.6198027589576429186024597584037 ) - 1 );
Initial_log_xi = -2.3237213962137608369573626987403;
Initial_log_thetaK = 2.0066560087511180476838035247056;
Initial_log_thetaI = -8.9717174095365983532701648073271;
Initial_log_Omega = -2.3011195297317619257171372737503;
Initial_log_AX = 0;
Initial_log_AF = 0;
Initial_log_AH = 0;
Initial_log_AZ = 0;
Initial_log_AL = 0;
Initial_log_AC = 0;
Initial_log_AG = 0;
Initial_log_AK = 0;
Initial_log_PY = 8.1350968392448024246732529718429;

@#if LoadCurrentBest
    load CurrentBest;
    M_.params = best_M_params;
@#endif

@#if LoadInitParams
    load_params_and_steady_state( 'InitParams.txt' );
@#endif

@#if Estimation

    varobs log_K log_H log_C log_I log_G log_NGDP log_PC log_PI log_PG log_R log_U log_N log_T log_G_GDP log_LS log_Omega;

    @#if !GeneralSteadyState
    
        @#for LogI1ShockProcess in ThetaVariables

            @#define ShockProcessName = "S" + LogI1ShockProcess
            
            gamma_@{LogI1ShockProcess} = 0;
            iota0_@{LogI1ShockProcess} = 0;
            iota1_@{LogI1ShockProcess} = 0;
            rho_@{ShockProcessName}    = 0;
            sigma_@{ShockProcessName}  = 0;
            
        @#endfor
        
    @#endif

    estimated_params;
    
        @#for Tax in Taxes

            @#define TaxVar = "tau" + Tax
            @#define EndoVariableName = TaxVar + "t"
            @#define ShockProcessName = "S" + EndoVariableName
            
            @#if Detrend && ( Tax != "GE" )
                @{EndoVariableName}_STEADY, @{EndoVariableName}_STEADY;
            @#endif
            
            phi_@{Tax}, phi_@{Tax};
            psi_@{Tax}, psi_@{Tax};
            chi_@{Tax}, chi_@{Tax};
            gamma_@{EndoVariableName}, gamma_@{EndoVariableName}, 0, 1;
            iota0_@{EndoVariableName}, iota0_@{EndoVariableName};
            iota1_@{EndoVariableName}, iota1_@{EndoVariableName};
            rho_@{ShockProcessName}, rho_@{ShockProcessName}, 0, 1;
            sigma_@{ShockProcessName}, sigma_@{ShockProcessName}, 0, Inf;
            
        @#endfor
        
        @#for PriceVar in PriceVars

            @#define PriceVarName = "P" + PriceVar
            
            @#if Detrend
                @{PriceVarName}_STEADY, @{PriceVarName}_STEADY, 0, Inf;
            @#endif

        @#endfor
        
        @#for LogI1ShockProcess in LogI1ShockProcesses

            @#define ShockProcessName = "S" + LogI1ShockProcess

            @#if Detrend && ( !( LogI1ShockProcess in EstimationFixedSteady ) )
                @{LogI1ShockProcess}_STEADY, @{LogI1ShockProcess}_STEADY, 0, Inf;
            @#endif
            
            @#if ( !Detrend ) || ( !( LogI1ShockProcess in ThetaVariables ) )
            
                gamma_@{LogI1ShockProcess}, gamma_@{LogI1ShockProcess}, 0, 1;
                iota0_@{LogI1ShockProcess}, iota0_@{LogI1ShockProcess};
                iota1_@{LogI1ShockProcess}, iota1_@{LogI1ShockProcess};
                
                @#if !Detrend
                    G@{LogI1ShockProcess}_STEADY, G@{LogI1ShockProcess}_STEADY;
                @#endif
                
                rho_@{ShockProcessName}, rho_@{ShockProcessName}, 0, 1;
                sigma_@{ShockProcessName}, sigma_@{ShockProcessName}, 0, 1;
            
            @#endif
            
        @#endfor
        
        @#for FixedVarianceShock in FixedVarianceShocks

            rho_@{FixedVarianceShock}, rho_@{FixedVarianceShock}, 0, 1;
            
        @#endfor
        
        @#for CESVar in CESVars
            
            a@{CESVar}, a@{CESVar}, 0, 1;
            s@{CESVar}, s@{CESVar}, 0, Inf;
            
        @#endfor
        
        eta, eta, 0, 1;
        
        @#if !Detrend
            hBar, hBar, 0, Inf;
        @#endif
        
        kappa, kappa, 0, 0.5;
        nu, nu, 0, Inf;
           
        Initial_level_tauGEt, Initial_level_tauGEt, -Inf, Inf;
        Initial_level_tauEt, Initial_level_tauEt, -Inf, Inf;
        Initial_level_tauHt, Initial_level_tauHt, -Inf, Inf;
        Initial_level_tauCt, Initial_level_tauCt, -Inf, Inf;
        Initial_level_tauKt, Initial_level_tauKt, -Inf, Inf;
        Initial_level_tauIt, Initial_level_tauIt, -Inf, Inf;
        Initial_level_tauYt, Initial_level_tauYt, -Inf, Inf;
        Initial_log_AX, Initial_log_AX, -Inf, Inf;
        Initial_log_AF, Initial_log_AF, -Inf, Inf;
        Initial_log_AH, Initial_log_AH, -Inf, Inf;
        Initial_log_PCt, Initial_log_PCt, -Inf, Inf;
        Initial_log_PIt, Initial_log_PIt, -Inf, Inf;
        Initial_log_PGt, Initial_log_PGt, -Inf, Inf;
        Initial_log_PEt, Initial_log_PEt, -Inf, Inf;
        Initial_log_N, Initial_log_N, -Inf, Inf;
        Initial_log_deltat, Initial_log_deltat, -Inf, Inf;
        Initial_log_xi, Initial_log_xi, -Inf, Inf;
        Initial_log_thetaK, Initial_log_thetaK, -Inf, Inf;
        Initial_log_thetaI, Initial_log_thetaI, -Inf, Inf;
        Initial_log_Omega, Initial_log_Omega, -Inf, Inf;
        Initial_log_AZ, Initial_log_AZ, -Inf, Inf;
        Initial_log_AL, Initial_log_AL, -Inf, Inf;
        Initial_log_AC, Initial_log_AC, -Inf, Inf;
        Initial_log_AG, Initial_log_AG, -Inf, Inf;
        Initial_log_AK, Initial_log_AK, -Inf, Inf;
        Initial_log_PY, Initial_log_PY, -Inf, Inf;
        
        @#if MeasurementError
        
            stderr log_K, 100, 0, Inf;
            stderr log_H, 100, 0, Inf;
            stderr log_C, 100, 0, Inf;
            stderr log_I, 100, 0, Inf;
            stderr log_G, 100, 0, Inf;
            stderr log_NGDP, 100, 0, Inf;
            stderr log_PC, 100, 0, Inf;
            stderr log_PI, 100, 0, Inf;
            stderr log_PG, 100, 0, Inf;
            stderr log_R, 100, 0, Inf;
            stderr log_U, 100, 0, Inf;
            stderr log_N, 100, 0, Inf;
            stderr log_T, 100, 0, Inf;
            stderr log_G_GDP, 100, 0, Inf;
            stderr log_LS, 100, 0, Inf;
            stderr log_Omega, 100, 0, Inf;
        
        @#endif
            
    end;
    
    @#if LoadCurrentBest && MeasurementError
        estim_params_.var_endo( 1:16, 2 ) = best_xparam1( 1:16 );
    @#endif
    
@#endif

@#if UseDLL

model( use_dll );

@#else

model;

@#endif

    @#include "InsertNewModelEquations.mod"
    
    @#if Detrend
        #GrowthSwitch = 0;
        #GrowthSwitch_LEAD = 0;
    @#else
        GrowthSwitch = GrowthSwitch(-1);
        #GrowthSwitch_LEAD = GrowthSwitch(+1);
    @#endif
    
    #V0 = 1;
    #X0 = 1;
    #Z0 = 1;
    #Y0 = 1;
    
    @#for CESProductivityVar in CESProductivityVars

        #@{CESProductivityVar}0 = 1;

    @#endfor
    
    #delta_LAG  = 1 / ( 1 + 1 / deltat_LAG );
    #delta      = 1 / ( 1 + 1 / deltat );
    #delta_LEAD = 1 / ( 1 + 1 / deltat_LEAD );

    #Kt0 = U0 * @{Kt0OverU0};

    #H0 = h0 * N0;

    #h_LAG  = ht_LAG * hBar;
    #h      = ht * hBar;
    #h_LEAD = ht_LEAD * hBar;

    @#for LogI1ShockProcess in LogI1ShockProcesses
    
        #G@{LogI1ShockProcess} = G@{LogI1ShockProcess}_STEADY ^ GrowthSwitch * S@{LogI1ShockProcess} / S@{LogI1ShockProcess}_LAG ^ gamma_@{LogI1ShockProcess} * A1 ^ iota1_@{LogI1ShockProcess} * A0 ^ iota0_@{LogI1ShockProcess} / A0_LAG ^ iota0_@{LogI1ShockProcess};
        #G@{LogI1ShockProcess}_LEAD = G@{LogI1ShockProcess}_STEADY ^ GrowthSwitch_LEAD * S@{LogI1ShockProcess}_LEAD / S@{LogI1ShockProcess} ^ gamma_@{LogI1ShockProcess} * A1_LEAD ^ iota1_@{LogI1ShockProcess} * A0_LEAD ^ iota0_@{LogI1ShockProcess} / A0 ^ iota0_@{LogI1ShockProcess};
        
        [ dynamic ] log( @{LogI1ShockProcess} ) = log( @{LogI1ShockProcess}_LAG ) + log( G@{LogI1ShockProcess} );
        [ static  ] log( @{LogI1ShockProcess} ) = log( @{LogI1ShockProcess}_STEADY );

    @#endfor
    
    @#for PriceVar in PriceVars

        @#define PriceVarName = "P" + PriceVar

        log( @{PriceVarName} ) = log( @{PriceVarName}t ) + log( PY );

    @#endfor
    
    @#for Tax in Taxes
        
        @#define TaxVar = "tau" + Tax
        @#define EndoVariableName = TaxVar + "t"
        
        #D@{EndoVariableName} = S@{EndoVariableName} - gamma_@{EndoVariableName} * S@{EndoVariableName}_LAG + iota1_@{EndoVariableName} * tau1 + iota0_@{EndoVariableName} * tau0 - iota0_@{EndoVariableName} * tau0_LAG;
        
        [ dynamic ] @{EndoVariableName} = @{EndoVariableName}_LAG + D@{EndoVariableName};
        [ static  ] @{EndoVariableName} = @{EndoVariableName}_STEADY;

        #@{TaxVar}      = 2 / ( 1 + exp( -( phi_@{Tax} * log( G_GDP      ) + psi_@{Tax} * log( h      / h_LAG ) + chi_@{Tax} * log( Ut      / Ut_LAG ) + @{EndoVariableName} ) ) ) - 1;
        #@{TaxVar}_LEAD = 2 / ( 1 + exp( -( phi_@{Tax} * log( G_GDP_LEAD ) + psi_@{Tax} * log( h_LEAD / h     ) + chi_@{Tax} * log( Ut_LEAD / Ut     ) + @{EndoVariableName} ) ) ) - 1;
        
    @#endfor
    
    #PK = PI / Omega;
    #PK_LEAD = PI_LEAD / Omega_LEAD;
    
    #Kt = Ut * ( K_LAG + kappa * Omega * I );
    #Kt_LEAD = Ut_LEAD * ( K + kappa * Omega_LEAD * I_LEAD );
    
    #Z = Z0 * ( aZ * ( AC * C * N0 / AC0 / N / C0 ) ^ ( ( sZ - 1 ) / sZ ) + ( 1 - aZ ) * ( AG * G * N0 / AG0 / N / G0 ) ^ ( ( sZ - 1 ) / sZ ) ) ^ ( sZ / ( sZ - 1 ) );
    #Z_LEAD = Z0 * ( aZ * ( AC_LEAD * C_LEAD * N0 / AC0 / N_LEAD / C0 ) ^ ( ( sZ - 1 ) / sZ ) + ( 1 - aZ ) * ( AG_LEAD * G_LEAD * N0 / AG0 / N_LEAD / G0 ) ^ ( ( sZ - 1 ) / sZ ) ) ^ ( sZ / ( sZ - 1 ) );

    #X = X0 * ( aX * ( AZ * Z / AZ0 / Z0 ) ^ ( ( sX - 1 ) / sX ) + ( 1 - aX ) * ( AL * max( 1e-8, hBar - h ) / AL0 / ( hBar - h0 ) ) ^ ( ( sX - 1 ) / sX ) ) ^ ( sX / ( sX - 1 ) );
    #X_LEAD = X0 * ( aX * ( AZ_LEAD * Z_LEAD / AZ0 / Z0 ) ^ ( ( sX - 1 ) / sX ) + ( 1 - aX ) * ( AL_LEAD * max( 1e-8, hBar - h_LEAD ) / AL0 / ( hBar - h0 ) ) ^ ( ( sX - 1 ) / sX ) ) ^ ( sX / ( sX - 1 ) );

    #SR = aY * PY / PK * Y0 * ( Y / Y0 ) ^ ( 1 / sY ) * ( AK * Kt / AK0 / Kt0 ) ^ ( ( sY - 1 ) / sY ) / Kt;
    #SR_LEAD = aY * PY_LEAD / PK_LEAD * Y0 * ( Y_LEAD / Y0 ) ^ ( 1 / sY ) * ( AK_LEAD * Kt_LEAD / AK0 / Kt0 ) ^ ( ( sY - 1 ) / sY ) / Kt_LEAD;

    #lambdaX = X ^ ( - 1 / sV );
    #lambdaX_LEAD = X_LEAD ^ ( - 1 / sV );

    #lambdaZ = aX * lambdaX * X0 * ( X / X0 ) ^ ( 1 / sX ) * ( AZ * Z / AZ0 / Z0 ) ^ ( ( sX - 1 ) / sX ) / Z;
    #lambdaZ_LEAD = aX * lambdaX_LEAD * X0 * ( X_LEAD / X0 ) ^ ( 1 / sX ) * ( AZ_LEAD * Z_LEAD / AZ0 / Z0 ) ^ ( ( sX - 1 ) / sX ) / Z_LEAD;

    #lambdaB = aZ * ( 1 - tauC ) * lambdaZ * Z0 * ( Z / Z0 ) ^ ( 1 / sZ ) * ( AC * C * N0 / AC0 / N / C0 ) ^ ( ( sZ - 1 ) / sZ ) * N / C;
    #lambdaB_LEAD = aZ * ( 1 - tauC_LEAD ) * lambdaZ_LEAD * Z0 * ( Z_LEAD / Z0 ) ^ ( 1 / sZ ) * ( AC_LEAD * C_LEAD * N0 / AC0 / N_LEAD / C0 ) ^ ( ( sZ - 1 ) / sZ ) * N_LEAD / C_LEAD;
    
    #Xi_LEAD = beta * AF ^ ( ( sV - 1 ) / sV ) * GN_LEAD ^ ( eta - 1 ) * GAX_LEAD ^ ( ( sV - 1 ) / sV ) * lambdaB_LEAD * PC / lambdaB / PC_LEAD;
    
    #Ht      = h * N;
    #Ht_LEAD = h_LEAD * N_LEAD;
    
    log( Y ) = log( Y0 * ( aY * ( AK * Kt / AK0 / Kt0 ) ^ ( ( sY - 1 ) / sY ) + ( 1 - aY ) * ( AH * Ht / AH0 / H0 ) ^ ( ( sY - 1 ) / sY ) ) ^ ( sY / ( sY - 1 ) ) );
    // Y_LEAD = Y0 * ( aY * ( AK_LEAD * Kt_LEAD / AK0 / Kt0 ) ^ ( ( sY - 1 ) / sY ) + ( 1 - aY ) * ( AH_LEAD * Ht_LEAD / AH0 / H0 ) ^ ( ( sY - 1 ) / sY ) ) ^ ( sY / ( sY - 1 ) );   

    #W = ( 1 - aX ) / ( 1 - tauH ) * lambdaX / lambdaB * X0 * ( X / X0 ) ^ ( 1 / sX ) * ( AL * max( 1e-8, hBar - h ) / AL0 / ( hBar - h0 ) ) ^ ( ( sX - 1 ) / sX ) / max( 1e-8, hBar - h );
    
    log( K ) = log( ( 1 - delta ) * K_LAG + ( 1 - delta / 2 ) * Omega * I );
        
    lambdaK + thetaK * log( K / K_LAG ) * K_LAG / K + Xi_LEAD * PK_LEAD / PK * ( thetaK_LEAD / 2 * log( K_LEAD / K ) ^ 2 + thetaI_LEAD / 2 * log( I_LEAD / I ) ^ 2 - xi_LEAD * log( 1 - Ut_LEAD ^ nu ) ) 
        = Xi_LEAD * PK_LEAD / PK * ( lambdaK_LEAD * ( 1 - delta_LEAD ) + ( 1 - tauK_LEAD ) * SR_LEAD * Ut_LEAD + thetaK_LEAD * log( K_LEAD / K ) );
    
    lambdaK * ( 1 - delta / 2 ) * Omega + ( 1 - tauK ) * SR * Ut * kappa * Omega + Xi_LEAD * PK_LEAD / PK * thetaI_LEAD * log( I_LEAD / I ) * K / I = 1 / ( 1 - tauI ) * Omega + thetaI * log( I / I_LAG ) * K_LAG / I;
    
    1 = R * Xi_LEAD;

    log( ( 1 - tauK ) * SR * Kt * ( 1 - Ut ^ nu ) ) = log( nu * xi * K_LAG * Ut ^ nu );
    
    log( W ) = log( ( 1 - aY ) * PY / PC * Y0 * ( Y / Y0 ) ^ ( 1 / sY ) * ( AH * Ht / AH0 / H0 ) ^ ( ( sY - 1 ) / sY ) / ( Ht ) );
    
    log( G ) = log( ( 1 + tauGE ) / 2 * ( 1 - tauE ) * PY * Y / PG );
    
    E = ( 1 + tauGE ) / 2 * tauE * PY * Y / PE;
    
    log( NGDP ) = log( PC * C + PI * I + PG * G + PE * E );
    
    log( G_GDP ) = log( sqrt( ( PC * C + PI * I + PG * G + PE * E ) / ( PC * C_LAG + PI * I_LAG + PG * G_LAG + PE * E_LAG ) * ( PC_LAG * C + PI_LAG * I + PG_LAG * G + PE_LAG * E ) / ( PC_LAG * C_LAG + PI_LAG * I_LAG + PG_LAG * G_LAG + PE_LAG * E_LAG ) ) );
    
    // GP_GDP = sqrt( ( PC * C + PI * I + PG * G + PE * E ) / ( PC_LAG * C + PI_LAG * I + PG_LAG * G + PE_LAG * E ) * ( PC * C_LAG + PI * I_LAG + PG * G_LAG + PE * E_LAG ) / ( PC_LAG * C_LAG + PI_LAG * I_LAG + PG_LAG * G_LAG + PE_LAG * E_LAG ) );

    log( PY * Y ) = log( PC * C + PI * I + PG * G + PE * E + thetaK / 2 * log( K / K_LAG ) ^ 2 * PK * K_LAG + thetaI / 2 * log( I / I_LAG ) ^ 2 * PK * K_LAG - xi * log( 1 - Ut ^ nu ) * PK * K_LAG ); // = log( PK * SR * Kt + PC * W * Ht )
    
    log( T ) = log( tauC / ( 1 - tauC ) * PC * C + tauI / ( 1 - tauI ) * PI * I + tauY * PY * Y + tauK * PK * SR * Kt + tauH * PC * W * Ht );
    
    // CFC = PK * ( delta * K_LAG + delta / 2 * I * Omega );
    
    log( LS ) = log( PC * W * Ht / NGDP );
    
    U = Ut;
    H = Ht;
    
    @#if ValueFunction
    
        #Ft0 = V0 ^ ( ( sV - 1 ) / sV ) * ( 1 - aV ) / beta;
        Ft = GN_LEAD ^ eta * ( AF * V_LEAD ) ^ ( ( sV - 1 ) / sV );    
        V = V0 * ( aV * ( AX * X / AX0 / X0 ) ^ ( ( sV - 1 ) / sV ) + ( 1 - aV ) * ( Ft / Ft0 ) ) ^ ( sV / ( sV - 1 ) );
    
    @#endif
    
end;

write_latex_dynamic_model;

steady_state_model;

    @#include "InsertNewStartSteadyStateEquations.mod"
    
    @#if !Detrend
        GrowthSwitch = 0;
    @#endif

    beta = 1 - aV;
    
    @#if !Detrend
    
        @#for LogI1ShockProcess in LogI1ShockProcesses
        
            @{LogI1ShockProcess}_STEADY = exp( Initial_log_@{LogI1ShockProcess} );
        
        @#endfor
        
        @#for Tax in Taxes
            
            @#define TaxVar = "tau" + Tax
            @#define EndoVariableName = TaxVar + "t"
            
            @{EndoVariableName}_STEADY = Initial_level_@{EndoVariableName};

        @#endfor
        
    @#endif
    
    @#for PriceVar in PriceVars

        @#define PriceVarName = "P" + PriceVar

        @{PriceVarName}t_STEADY = @{PriceVarName}_STEADY / PY_STEADY;
        @{PriceVarName}_        = @{PriceVarName}_STEADY;

    @#endfor

    @#if !GeneralSteadyState
    
        @#for CESProductivityVar in CESProductivityVars
        
            @{CESProductivityVar}_STEADY = 1;
            
        @#endfor    
    
    @#endif
    
    @#for LogI1ShockProcess in LogI1ShockProcesses

        G@{LogI1ShockProcess}_ = 1; // S@{LogI1ShockProcess}_;
        @{LogI1ShockProcess}_ = @{LogI1ShockProcess}_STEADY;

    @#endfor

    @#for Tax in Taxes
        
        @#define TaxVar = "tau" + Tax
        @#define EndoVariableName = TaxVar + "t"
        
        @{EndoVariableName}_ = @{EndoVariableName}_STEADY;
        @{TaxVar}_ = 2 / ( 1 + exp( - @{EndoVariableName}_ ) ) - 1;

    @#endfor
    
    delta_      = 1 / ( 1 + 1 / deltat_ );
    
    PK_ = PI_ / Omega_;

    Xi_ = beta * AF_ ^ ( ( sV - 1 ) / sV ) * GN_ ^ ( eta - 1 ) * GAX_ ^ ( ( sV - 1 ) / sV );
    R_ = 1 / Xi_;

    H0_ = h0 * N0;
    
    Kt0_ = U0 * @{Kt0OverU0};

    @#if GeneralSteadyState

        Y_ = SolveSteady( hBar, kappa, nu, aX, aZ, aY, sV, sX, sZ, sY, Xi_, N_, delta_, xi_, Omega_, PY_, PC_, PI_, PG_, PE_, AZ_, AL_, AC_, AG_, AK_, AH_, H0_, Kt0_, h0, C0, G0, N0, 1, 1, 1, 1, 1, 1, tauGE_, tauE_, tauH_, tauC_, tauK_, tauI_ );

        t1 = (nu * sY);
        t2 = delta_ ^ 2;
        t3 = (t2 * Xi_);
        t4 = (kappa * t3);
        t6 = (nu * tauI_ * xi_);
        t9 = (nu * kappa);
        t10 = (xi_ * t9);
        t14 = (Xi_ * delta_);
        t15 = (kappa * t14);
        t18 = (xi_ * nu);
        t24 = (Xi_ * kappa);
        t27 = (delta_ * kappa);
        t35 = 6 * t10 * t14 - 2 * t10 * t3 - 4 * t18 * t14 + 4 * t6 * t14 - 6 * t6 * t15 - 4 * t18 * t24 + t18 * t3 + 4 * t6 * t24 + 2 * t6 * t27 - t6 * t3 + 2 * t6 * t4 + 4 * t4;
        t36 = (Xi_ * nu);
        t37 = (xi_ * tauI_);
        t53 = -2 * t18 * t27 - 4 * t37 * t36 + 4 * xi_ * t36 - 4 * t37 * t9 - 4 * Xi_ - 2 * delta_ + 4 * t10 + 6 * t14 - 4 * t15 + 4 * t27 - 2 * t3 + 4;
        t55 = 1 / Xi_;
        t61 = -2 + delta_;
        t66 = 1 / (2 * t27 - delta_ + 2);
        t69 = exp(-(t66 / t61 / (-1 + tauI_) / xi_ * t55 * (t35 + t53)));
        t79 = lambertwM1(-t55 * t66 * nu * (2 * t15 - t14 - 2 * t24 + 2 * Xi_ + 2 * kappa) * t69);
        t88 = log(t66 * nu * t69 / (-1 + tauK_) * xi_ * t61);
        t90 = t79 * Xi_;
        t107 = log(t55 * t66 / t79 * (-(delta_ * t90) - nu * t14 - 0.2e1 * nu * t24 + (2 * t9 * t14) + (2 * t27 * t90) + (2 * t36) + (2 * t9) + (2 * t90)));
        t108 = sY * t107;
        t114 = (AK_ / Kt0_) ^ sY;
        t120 = (aY * PY_ / PI_ * Omega_) ^ sY;
        t122 = log(t120 * Y_ * t114);
        t125 = sY - 0.1e1;
        t127 = 0.1e1 / nu;
        t128 = 0.1e1 / sY;
        t131 = exp(t128 * t127 * t125 * (-nu * t108 + nu * t122 + (t79 * t1) - t88 * t1 + t108));
        t134 = Y_ ^ (t128 * t125);
        t141 = (-0.1e1 / (-1 + aY) * (-t131 * aY + t134)) ^ (0.1e1 / t125 * sY);
        H_ = min( N_ * hBar - 1e-8, H0_ / AH_ * t141 );
        U_ = exp(t127 * t107);

        h_ = H_ / N_;
       
        SR_ = (-xi_*Xi_*PK_*(-1+tauI_)*(-2+delta_)*log(1-U_^nu)+2*Xi_*(delta_-1)*PK_+2*PK_)/((2*(-1+tauK_))*(-1+tauI_)*(Xi_*((kappa-1/2)*delta_-kappa+1)*PK_+PK_*kappa)*U_);
        I_ = ( aY * PY_ / PK_ ) ^ sY * ( Y_ ) * ( AK_ / 1 / Kt0_ ) ^ ( sY - 1 ) / ( U_ * ( ( 1 - delta_ / 2 ) / delta_ * Omega_ + kappa * Omega_ ) ) / SR_ ^ sY;
        K_ = ( 1 - delta_ / 2 ) / delta_ * Omega_ * I_;
        Kt_ = U_ * ( K_ + kappa * Omega_ * I_ );
        G_ = ( 1 + tauGE_ ) / 2 * ( 1 - tauE_ ) * PY_ * Y_ / PG_;
        E_ = ( 1 + tauGE_ ) / 2 * tauE_ * PY_ * Y_ / PE_;
        C_ = max( 1e-8, 1 / PC_ * ( PY_ * Y_ - ( PI_ * I_ + PG_ * G_ + PE_ * E_ - xi_ * log( 1 - U_ ^ nu ) * PK_ * K_ ) ) );
        Z_ = ( aZ * ( AC_ * C_ * N0 / 1 / N_ / C0 ) ^ ( ( sZ - 1 ) / sZ ) + ( 1 - aZ ) * ( AG_ * G_ * N0 / 1 / N_ / G0 ) ^ ( ( sZ - 1 ) / sZ ) ) ^ ( sZ / ( sZ - 1 ) );
        X_ = ( aX * ( AZ_ * Z_ / 1 ) ^ ( ( sX - 1 ) / sX ) + ( 1 - aX ) * ( AL_ * max( 1e-8, hBar - h_ ) / 1 / ( hBar - h0 ) ) ^ ( ( sX - 1 ) / sX ) ) ^ ( sX / ( sX - 1 ) );
        lambdaX_ = X_ ^ ( - 1 / sV );
        lambdaZ_ = aX * lambdaX_ * ( X_ ) ^ ( 1 / sX ) * ( AZ_ * Z_ / 1 ) ^ ( ( sX - 1 ) / sX ) / Z_;
        lambdaB_ = aZ * ( 1 - tauC_ ) * lambdaZ_ * ( Z_ ) ^ ( 1 / sZ ) * ( AC_ * C_ * N0 / 1 / N_ / C0 ) ^ ( ( sZ - 1 ) / sZ ) * N_ / C_;
        
        // To check
        // Y_ = ( aY * ( AK_ * Kt_ / 1 / Kt0_ ) ^ ( ( sY - 1 ) / sY ) + ( 1 - aY ) * ( AH_ * H_ / 1 / H0_ ) ^ ( ( sY - 1 ) / sY ) ) ^ ( sY / ( sY - 1 ) );
        // ( 1 - aY ) * PY_ / PC_ * ( Y_ ) ^ ( 1 / sY ) * ( AH_ * H_ / 1 / H0_ ) ^ ( ( sY - 1 ) / sY ) / H_ = ( 1 - aX ) / ( 1 - tauH_ ) * lambdaX_ / lambdaB_ * ( X_ ) ^ ( 1 / sX ) * ( AL_ * max( 1e-8, hBar - h_ ) / 1 / ( hBar - h0 ) ) ^ ( ( sX - 1 ) / sX ) / max( 1e-8, hBar - h_ );
        // ( 1 - tauK_ ) * SR_ * Kt_ * ( 1 - U_ ^ nu ) = nu * xi_ * K_ * U_ ^ nu;

    @#else
    
        Y_ = 1;
        
        N_ = N0;
        N_STEADY = N0;
        
        Omega_ = Omega0;
        Omega_STEADY = Omega0;

        lambdaX_ = 1;
        lambdaZ_ = aX;
        Z_ = 1;
        X_ = 1;
        H_ = H0_;
        U_ = U0;
        C_ = C0;
        G_ = G0;
        Kt_ = Kt0_;
        h_ = H_ / N_;
        
        delta_ = SolveSteadySimple( kappa, nu, aY, sY, Xi_, Omega_, PC_, PI_, PG_, PK_, C_, G_, tauK_, tauI_, tauE_, U_, Kt_ );
        
        t4 = Xi_ * delta_ - Xi_ + 1;
        t6 = (U_ ^ nu);
        t7 = -1 + t6;
        t9 = (kappa - 0.1e1 / 0.2e1) * delta_;
        t1_ = 0.1e1 + t9;
        t11 = t7 * t1_;
        t16 = log(-t7);
        t25 = 0.1e1 / (-Xi_ * t7 * t1_ * t16 + t6 * nu * (Xi_ * (t9 - kappa + 0.1e1) + kappa));
        t27 = 1 / (-1 + tauI_);
        t28 = t25 * t27;
        t29 = -1 + aY;
        t32 = -1 + tauH_;
        t33 = tauE_ - 1;
        t34 = t32 * t33;
        t35 = -1 + tauC_;
        t42 = C_ * t33 * PC_;
        t43 = (PG_ * G_);
        t44 = t42 - t43;
        t52 = t1_ * Omega_;
        t54 = Kt0_ * PI_;
        hBar = h_ * (-Kt0_ * PK_ * aX * aZ * t4 * t11 * t28 * Omega_ * t29 * t34 * t35 * t16 + ((t29 * t35 * t44 * t32 * aZ + t42) * aX - t42) * U_ * t52 + t54 * aX * aZ * delta_ * t29 * t34 * t35) / (-PK_ * t4 * t11 * t25 * t27 * Kt0_ * Omega_ * t33 * t16 + t44 * U_ * t52 + t54 * delta_ * t33) / t29 / t35 / aZ / t32 / aX;
        t88 = t4 * t7 * t1_;
        t9_ = 1 / (-2 + delta_);
        xi_ = -0.2e1 * t88 * t28 * t9_;

        I_ = Kt_ / U_ / ( ( 1 - delta_ / 2 ) / delta_ * Omega_ + kappa * Omega_  );
        K_ = ( 1 - delta_ / 2 ) / delta_ * Omega_ * I_;
        PY_ = PC_ * C_ + PI_ * I_ + PG_ * G_ - xi_ * log( 1 - U_ ^ nu ) * PK_ * K_ + PG_ * G_ * tauE_ / ( 1 - tauE_ );
        tauGE_ = G_ / ( ( 1 - tauE_ ) * PY_ / PG_ );
        tauGEt_ = -log( 2 / ( tauGE_ + 1 ) - 1 );
        E_ = tauGE_ * tauE_ * PY_ * Y_ / PE_;
        SR_ = (-xi_*Xi_*PK_*(-1+tauI_)*(-2+delta_)*log(1-U_^nu)+2*Xi_*(delta_-1)*PK_+2*PK_)/((2*(-1+tauK_))*(-1+tauI_)*(Xi_*((kappa-1/2)*delta_-kappa+1)*PK_+PK_*kappa)*U_);
        lambdaB_ = aZ * ( 1 - tauC_ ) * lambdaZ_ * ( Z_ ) ^ ( 1 / sZ ) * ( AC_ * C_ * N0 / 1 / N_ / C0 ) ^ ( ( sZ - 1 ) / sZ ) * N_ / C_;

        deltat_ = 1 / ( 1 / delta_ - 1 );

        xi_STEADY = xi_;
        deltat_STEADY = deltat_;
        tauGEt_STEADY = tauGEt_;
        PY_STEADY = PY_;

    @#endif
    
    @#for PriceVar in PriceVars

        @#define PriceVarName = "P" + PriceVar

        @{PriceVarName}t_STEADY = @{PriceVarName}_STEADY / PY_STEADY;
        @{PriceVarName}t_       = @{PriceVarName}t_STEADY;

    @#endfor

    lambdaK_ = Xi_*(-1+xi_*kappa*(-1+tauI_)*log(1-U_^nu))*PK_/(((Xi_*(delta_-1)*PK_+PK_)*kappa-(1/2)*PK_*Xi_*(-2+delta_))*(-1+tauI_));
    W_ = ( 1 - aY ) * PY_ / PC_ * ( Y_ ) ^ ( 1 / sY ) * ( AH_ * H_ / 1 / H0_ ) ^ ( ( sY - 1 ) / sY ) / H_;
    T_ = tauC_ / ( 1 - tauC_ ) * PC_ * C_ + tauI_ / ( 1 - tauI_ ) * PI_ * I_ + tauY_ * PY_ * Y_ + tauK_ * PK_ * SR_ * Kt_ + tauH_ * PC_ * W_ * H_;    
    G_GDP_ = 1;
    NGDP_ = PC_ * C_ + PI_ * I_ + PG_ * G_ + PE_ * E_;
    // CFC_ = PK_ * ( delta_ * K_ + delta_ / 2 * I_ * Omega_ );
    LS_ = PC_ * W_ * H_ / NGDP_;
    
    Ut_ = U_;
    ht_ = h_ / hBar;
    
    @#if ValueFunction
    
        Ft0_ = ( 1 - aV ) / beta;
        V_ = ( aV * ( AX_ * X_ / 1 ) ^ ( ( sV - 1 ) / sV ) / ( 1 - ( 1 - aV ) / Ft0_ * ( GN_ ^ eta * ( AF_ ) ^ ( ( sV - 1 ) / sV ) ) ) ) ^ ( sV / ( sV - 1 ) );
        Ft_ = GN_ ^ eta * ( AF_ * V_ ) ^ ( ( sV - 1 ) / sV );    
    
    @#endif
    
    @#include "InsertNewEndSteadyStateEquations.mod"

end;

shocks;
    @#include "InsertNewShockBlockLines.mod"
end;

@#define OriginalVariableString = "log_C log_I log_G log_Y log_G_GDP log_H log_R log_U log_LS"

@#define VariableString = OriginalVariableString + " log_A0 log_A1 level_tau0 level_tau1"

@#for LogI1ShockProcess in LogI1ShockProcesses

    @#define VariableString = VariableString + " log_S" + LogI1ShockProcess
    
@#endfor

@#for Tax in Taxes
 
    @#define VariableString = VariableString + " level_Stau" + Tax + "t"
    
@#endfor

@#if !Detrend

    options_.non_bgp = 1;
    options_.accurate_nonstationarity = 1;
    options_.accurate_nonstationarity_step_width = 0.1;
    options_.k_order_solver = 0;

    @#if GrowthIterations
        options_.non_bgp_growth_iterations = 20;
        options_.dynatol.f = 1e300;
    @#else
        options_.non_bgp_growth_iterations = 0;
    @#endif

@#endif

options_.qz_criterium = 1 + 1e-4;
options_.endogenous_qz_criterium = 1;

steady;

check;

save_params_and_steady_state( 'InitParams.txt' );

@#if Estimation

    @#if CMAES
    
        @#define OptimizationOptions = "mode_compute = 9, optim = ( 'ResumeRun', 0, 'ResumeFromBest', 1, 'SigmaScale', 1e-2, 'MinSigma', 1e-2, 'MaxFunEvals', 1e12, 'MaxIter', 1e12, 'TolFun', 1e-12, 'TolX', 1e-12, 'UseParallel', 1 )"
        
    @#else
    
        @#define OptimizationOptions = "mode_compute = 1, optim = ( 'Display', 'iter', 'MaxFunEvals', 1e12, 'MaxIter', 1e12, 'TolFun', 1e-12, 'TolX', 1e-12, 'UseParallel', 1 )"

    @#endif

    @#if Detrend

        @#define DataFile = "DetrendedData.mat"
        
    @#else
    
        @#define DataFile = "FinalData.mat"
        options_.extended_kalman_filter = 0;
    
    @#endif
    
    options_.add_empty_presamples = 10;

    estimation( datafile = 'SourceData/@{DataFile}', plot_priors = 0, graph_format = none, lik_init = 6, mh_replic = 0, mh_nblocks = 0, mode_check, @{OptimizationOptions}, smoother, forecast = 100, order = 1, use_univariate_filters_if_singularity_is_detected = 0, keep_kalman_algo_if_singularity_is_detected ) @{VariableString};

@#endif

@#if IRFLength > 0

    stoch_simul( order = 1, periods = 0, irf = @{IRFLength}, nocorr, nodecomposition, nofunctions, nomoments, graph_format = none ) @{VariableString};
    
@#endif

@#if SimulationLength > 0

    M_.Sigma_e = zeros( size( M_.Sigma_e ) );

    stoch_simul( order = 1, periods = @{SimulationLength}, drop = 0, irf = 0, nocorr, nodecomposition, nofunctions, nomoments, graph_format = none ) @{OriginalVariableString};
    
    [ ~, VarIndices ] = ismember( cellstr( var_list_ ), cellstr( M_.endo_names ) );
    figure;
    for i = 1 : length( VarIndices )
        subplot( 3, 3, i );
        plot( oo_.endo_simul( VarIndices( i ), : ).' );
        title( var_list_( i, : ) );
    end
    
@#endif
