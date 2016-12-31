//  Calculate porosity field (P)

{
    P = 1;
    forAll(FList, n)
    {
        volScalarField& F = FList[n]();
        P = P - F;
    } 
    P.correctBoundaryConditions();
}
    
