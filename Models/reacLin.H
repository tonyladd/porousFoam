//  Linear reaction kinetics  //


    double stoiTable[1][2] = {{-1.0,-1.0}};  // Stoichiometric coeffs

    RList[0]()  == k0*CList[0]();
    dRList[0]() == k0;

