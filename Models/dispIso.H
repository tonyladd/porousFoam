//  Isotropic dispersion  (with porosity) //
    
{
    symmTensor T(1,0,0,1,0,1);
    D = D0*P*T;
    D.correctBoundaryConditions();
}
