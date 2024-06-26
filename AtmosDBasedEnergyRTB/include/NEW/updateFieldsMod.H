{
    // Compute U
    U = rhoU/rho;
    U.correctBoundaryConditions();

    // Calculate enthalpy from rhoE
    //const volScalarField Cp = thermo->Cp();
    //const volScalarField Cv = thermo->Cv();
    
   
    h = Cp/Cv*(rhoE/rho - 0.5*magSqr(U));
    T = h/Cp;
//     h = rhoE/rho - 0.5*magSqr(U) + p/rho;
    // Alternative formulation, Felipe Alves Portela TU Delft
    /*
    h = Cp/Cv*(rhoE/rho - 0.5*magSqr(U));
    h.correctBoundaryConditions();
    */

    // Bound enthalpy
    /*
    dimensionedScalar CpMin = min(Cp);
    dimensionedScalar CpMax = max(Cp);

    dimensionedScalar hMin = CpMin*TMin;
    dimensionedScalar hMax = CpMax*TMax;
    
    boundMinMax(h, hMin, hMax);
    */

    // Bound density
    // boundMinMax(rho, rhoMin, rhoMax);
   
    p = rho*(R)*T;
    p.correctBoundaryConditions();

    // Bound pressure
    // boundMinMax(p, pMin, pMax);
}
