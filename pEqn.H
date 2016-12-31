//  SIMPLE pressure solver

{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::interpolate(HbyA) & mesh.Sf()
    );

    adjustPhi(phiHbyA, U, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    // Relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
}
