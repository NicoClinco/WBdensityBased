scalar theta0(300.0);
scalar Pg(100000.0);
scalar g_(9.81);
scalar CpOR = 3.5087;
scalar CvOR = 2.5087;
scalar CpTheta = Cp.value()*theta0;
scalar RTheta = R.value()*theta0;
scalar PgoRTheta= Pg/RTheta;
vectorField VV = mesh.C();

// Initialize the hydrostatic field:
forAll(T0,cI)
{
  // Density
  scalar Par = 1.0 - (g_*VV[cI].y())/CpTheta;
  scalar PgrRTheta = Pg/RTheta;
  rho0[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  rho[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  // Temperature 
 T0[cI] = theta0 - (g_*VV[cI].y())/Cp.value();
 T[cI] = theta0 - (g_*VV[cI].y())/Cp.value(); 
  // The pressure is set into the class numericFlux,
  // calculated from pressure and temperature
}
T.correctBoundaryConditions();
rho.correctBoundaryConditions();

Info << "Done" << endl;
