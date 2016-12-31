//  Darcy or Darcy-Brinkman flow solver

    U   == U*Uscale;                // Undo scaling
    phi == phi*Uscale;              // Improves solver convergence

    if (DARCY)                      //  Darcy-flow solver
    {
        Info << "Start Darcy solver" << endl;

        #include "darcy.H"

        phi == -fvc::interpolate(1/Kinv)*fvc::snGrad(p)*mag(mesh.Sf());
        U   == -(1/Kinv)*fvc::grad(p);

    }
    else
    {
        Info << "Start Brinkman solver" << endl;

        while (simple.loop())
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }
    }

    if (DEBUG)
    {
        #include "continuityErrs.H"
    }

    #include "avgU.H"               // Average flow velocity (Ux)

    U   == U/Uscale;                // Normalize U
    phi == phi/Uscale;              // Normalize surface volume flux

    Info << "Flow solver done: Uscale = " << Uscale.value() << nl
         << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
    	 << "    ClockTime = " << runTime.elapsedClockTime() << " s"
    	 << endl;
