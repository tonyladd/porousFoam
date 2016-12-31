    List< autoPtr<volScalarField> > FTmpList(nF);

    forAll(FTmpList, n)
    {
        FTmpList[n].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "FTmp" + name(n),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                FList[n]()
            ) 
        );
    }

