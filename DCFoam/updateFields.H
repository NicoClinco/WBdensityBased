{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

      
    // Internal Energy
    h = (rhoE/rho - 0.5*magSqr(U));

    h.correctBoundaryConditions();

    T = h/Cv;
    T.correctBoundaryConditions();
    forAll(T,cI)
    {
      thetaP[cI] = T[cI] + g_*mesh.C()[cI].y()/Cp.value()-theta0;
    }
    thetaP.correctBoundaryConditions();
    

    p = (gamma -1.0)*h*rho;
    rho.correctBoundaryConditions();
    
}
