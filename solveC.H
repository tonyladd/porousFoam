//  Solve the convection-diffusion equation (steady) for C

{
    Info << "Start C solve" << endl;

    int counter = 0;
    while(counter < nCycle)
    {
        counter++;
        double residual = 0;
        forAll(CList, n)
        {   
            label Case = n + nF;
            volScalarField& C = CList[n]();
            #include "dispModel.H"
            #include "kinetics.H"
            fvScalarMatrix CEqn
            (
                fvm::div(phi,C)
              - fvm::laplacian(D,C)
              - fvm::Sp(dR,C) - (R-dR*C)
            );
            CEqn.relax();
            double residualC = solve(CEqn).initialResidual();
            residual = std::max(residual, residualC);
        }

        if (DEBUG)
            Info << "Iteration " << counter
                 << ": residual = " << residual << endl;

        if (residual < tolerance)
        {
            Info << "Convection-diffusion: Converged in "
                 << counter << " steps" << endl;
            break;                                  // Done
        }
        else if (counter == nCycle)
        {
            Info << nl 
            << "Convection-diffusion: FAILED to converge in "
            << nCycle << " steps" << nl << endl;
            exit(1);                                // No convegence 
        }
    }

    Info << "Concentration solver done:" << nl
         << "ExecutionTime = " << runTime.elapsedCpuTime()
         << " s" << ": "
         << "ClockTime = " << runTime.elapsedClockTime() 
         << " s" << endl;
}
