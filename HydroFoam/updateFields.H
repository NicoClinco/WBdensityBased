{
    // Compute U
    U.internalField() = rhoU.internalField()/rho.internalField();
    U.correctBoundaryConditions();

     
    // Internal Energy
    h.internalField() = (rhoE.internalField()/rho.internalField() - 0.5*magSqr(U.internalField()));

    h.correctBoundaryConditions();

    T = h/Cv;
    T.correctBoundaryConditions();
    forAll(T,cI)
    {
      thetaP[cI] = T[cI] + g_*mesh.C()[cI].y()/Cp.value()-theta0;
    }
    thetaP.correctBoundaryConditions();
    
    /*
    p = rho*(R)*T;
    p.correctBoundaryConditions();
    */
    
    p = (gamma -1.0)*h*rho;
    CorrectHydroBCs(p,rho);
    // rhoE.boundaryField() = Cv*T.boundaryField();
    
    
}
