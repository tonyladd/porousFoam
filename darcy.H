//  Darcy flow solver

{
    forAll(Kinv, i)
    {
        if (Kinv[i] < 1.0/Kmax)  Kinv[i] = 1.0/Kmax;
    }
    Kinv.correctBoundaryConditions();

    int counter=0;
    while(counter < nCycle)
    {
        counter++;
        double residual = solve
        ( fvm::laplacian(1./Kinv,p) ).initialResidual();

        if (DEBUG)
        Info << "Iteration "    << counter 
             << ": residual = " << residual << endl;

        if (residual < tolerance)
        {
            Info << "Darcy solver: Converged in "
                 << counter << " steps" << endl;
            break;                                  // Done
        }                
        else if (counter == nCycle)
        {
            Info << "Darcy: FAILED to converge in "
                 << nCycle << " steps" << nl << endl;
            exit(1);                                // No convergence 
        }
    }
}
