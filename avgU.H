//  Calculate average flow velocity <Ux>

{
    const label patchID = mesh.boundaryMesh().findPatchID("outlet");
    dimensionedScalar Area = gSum(mesh.magSf().boundaryField()[patchID]);
    dimensionedScalar Q = gSum(phi.boundaryField()[patchID]);
    dimensionedScalar avgU = Q/Area;

//  Calculate velocity scaling factor (Uscale)

    if (Uscale.value() == 0)
        Uscale.value() =  avgU.value();     // Initial value of <U>

    if (QCON)
        Uscale.value() = avgU.value();
    else if (avgU.value() > Qmax*Uscale.value())
        Uscale.value() = avgU.value()/Qmax;
}
