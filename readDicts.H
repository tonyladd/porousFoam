//  Read additional dictionaries

    Info << "***************************************************" << nl;
    Info << "Reading porousFoam dictionary" << nl;
    IOdictionary porousFoam
    (
        IOobject
        (
            "porousFoamDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
                                                       
    int nF;
    if( !porousFoam.readIfPresent<int>("nF", nF) ){
    SeriousErrorIn("main")
    << "There is no nF parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    int nC;
    if( !porousFoam.readIfPresent<int>("nC", nC) ){
    SeriousErrorIn("main")
    << "There is no nC parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    int nR;
    if( !porousFoam.readIfPresent<int>("nR", nR) ){
    SeriousErrorIn("main")
    << "There is no nR parameter in porousFoam dictionary"
    << exit(FatalError);
    }

    bool QCON;
    if( !porousFoam.readIfPresent<bool>("QCON", QCON) ){
    SeriousErrorIn("main")
    << "There is no QCON parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    bool DARCY;
    if( !porousFoam.readIfPresent<bool>("DARCY", DARCY) ){
    SeriousErrorIn("main")
    << "There is no DARCY parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    bool DEBUG;
    if( !porousFoam.readIfPresent<bool>("DEBUG", DEBUG) ){
    SeriousErrorIn("main")
    << "There is no DEBUG parameter in porousFoam dictionary"
    << exit(FatalError);
    }

    double Kmax;
    if( !porousFoam.readIfPresent<double>("Kmax", Kmax) ){
    SeriousErrorIn("main")
    << "There is no Kmax parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    double Qmax;
    if( !porousFoam.readIfPresent<double>("Qmax", Qmax) ){
    SeriousErrorIn("main")
    << "There is no Qmax parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    double Pmin;
    if( !porousFoam.readIfPresent<double>("Pmin", Pmin) ){
    SeriousErrorIn("main")
    << "There is no Pmin parameter in porousFoam dictionary"
    << exit(FatalError);
    }

    Info << "PorousFoamDict:              "           << nl;
    Info << "# of minerals :              " << nF     << nl;
    Info << "# of concentrations :        " << nC     << nl;
    Info << "# of reactions :             " << nR     << nl;
    Info << "Constant Flow Rate:          " << QCON   << nl;
    Info << "Darcy Flow:                  " << DARCY  << nl;
    Info << "Debugging:                   " << DEBUG  << nl;
    Info << "Max Perm (Darcy flow):       " << Kmax   << nl;
    Info << "Max Flow Rate (QCON false):  " << Qmax   << nl;
    Info << "Min Porosity:                " << Pmin   << nl << nl;

    Info << "Reading fvSolution controls dictionary\n" << endl;
    const dictionary& controls = mesh.solutionDict().subDict("concentration");

    int nCycle;
    if( !controls.readIfPresent<int>("nCycle", nCycle) ){
    SeriousErrorIn("main")
    << "There is no nCycle parameter in porousFoam dictionary"
    << exit(FatalError);
    }
    double tolerance;
    if( !controls.readIfPresent<double>("tolerance", tolerance) ){
    SeriousErrorIn("main")
    << "There is no tolerance parameter in porousFoam dictionary"
    << exit(FatalError);
    }

    Info << "fvSolution.concentrationDict:"              << nl;
    Info << "# cycles:                    " << nCycle    << nl;
    Info << "tolerance:                   " << tolerance << nl << nl;

    Info << "Reading transportProperties\n" << nl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar K0 (transportProperties.lookup("K0"));
    dimensionedScalar D0 (transportProperties.lookup("D0"));
    dimensionedScalar k0 (transportProperties.lookup("k0"));
    scalar Da=1;
    if( !transportProperties.readIfPresent<scalar>("Da", Da) ){
        Info << "No Da in transportProperties : set Da = 1" << endl;
    }


    Info << "K0 = K0/L^2:  " << K0 << endl;
    Info << "D0 = U0*L/Pe: " << D0 << endl;
    Info << "k0 = Da*U0/L: " << k0 << endl;
    Info << "Da: " << Da << endl;

    Info << "***************************************************" << nl
         << endl;

