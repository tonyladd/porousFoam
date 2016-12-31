dimensionedScalar Uscale("Uscale", dimless, scalar(0));

Info << "Create velocity scale" << endl;
{
    volScalarField Ptmp = P;        // Save current P

    if (runTime.value() == 0)
    {
        #include "porosity.H"       // Calculate P(t=0)
    }
    else
    {
        volScalarField P0           // Read P(t=0)
        (
            IOobject
            (
                "P",
                "0",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        P = P0;
    }


//  Create Uscale 

    #include "permModel.H"
    #include "solveU.H"


//  Restore current mineral fraction and velocity fields

    if (runTime.value() != 0)
    {
        P = Ptmp;                   // Restore P
    }
}
