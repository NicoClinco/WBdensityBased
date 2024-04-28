{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

      
    // Internal Energy
    h = (rhoE/rho - 0.5*magSqr(U));

    h.correctBoundaryConditions();

    T = h/Cv;
    T.correctBoundaryConditions();
	//volScalarField dT = T - T0; // Differences in the temperature.
	
    forAll(T,cI)
    {
	  //thetaP[cI] = dT[cI];
      thetaP[cI] = T[cI] + g_*mesh.C()[cI].y()/Cp.value()-theta0;
    }
    thetaP.correctBoundaryConditions();
    
    /*
    p = rho*(R)*T;
    p.correctBoundaryConditions();
    */
    p = (gamma -1.0)*h*rho;
    rho.correctBoundaryConditions();
    // rhoE.boundaryField() = Cv*T.boundaryField();
    
}
