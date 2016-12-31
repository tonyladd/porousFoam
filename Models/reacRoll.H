//  Dissolution-precipitation: minimal model //
//  M_S + X <--> A_S            R1
//  M_P + X  --> A_P            R2
//  Ordering: M_S, M_P, X, A_S


    double stoiTable[2][4] =      // Stoichiometric numbers
        {{-1.0, 0.0,-1.0, 1.0},   // R1
         { 0.0,-1.0,-1.0, 0.0}};  // R2

    double kS  = 1.0;             // Rate constants
    double kP  = 10.0;
    double Keq = 1.0;

//  Rates

    RList[0]()  = k0*kS*(CList[0]()-CList[1]()/Keq);
    RList[1]()  = k0*kP*CList[0]();

//  Derivatives for each case

    switch (Case){
        case 0:
            dRList[0]() == k0*scalar(0);
            dRList[1]() == k0*scalar(0);
            break;
        case 1:
            dRList[0]() == k0*scalar(0);
            dRList[1]() == k0*scalar(0);
            break;
        case 2:
            dRList[0]() == k0*kS;
            dRList[1]() == k0*kP;
            break;
        case 3:
            dRList[0]() ==-k0*kS/Keq;
            dRList[1]() == k0*scalar(0);
            break;
    }

