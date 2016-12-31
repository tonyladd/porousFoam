//  Momentum predictor: Brinkmann

    fvVectorMatrix UEqn
    (
        fvm::Sp(Kinv,U)
      - fvm::laplacian(U)
    );

    UEqn.relax();
    solve(UEqn == -fvc::grad(p));

