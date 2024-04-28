scalar theta0(300.0);
scalar Pg(100000.0);
scalar g_(9.81);
//scalar CpOR = 3.5087;
scalar CvOR = 2.5087;
scalar CpTheta = Cp.value()*theta0;
scalar RTheta = R.value()*theta0;
//scalar PgoRTheta= Pg/RTheta;
vectorField VV = mesh.C();
//const surfaceVectorField& faceCentre = mesh.Cf();

// Initialize the hydrostatic field:
forAll(T,cI)
{
  // Density
  scalar Par = 1.0 - (g_*VV[cI].y())/CpTheta;
  scalar PgrRTheta = Pg/RTheta;
  // rho0[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  rho[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  // Temperature 
 //T0[cI] = theta0 - (g_*VV[cI].y())/Cp.value();
 T[cI] = theta0 - (g_*VV[cI].y())/Cp.value(); 
}


scalar r0 = 250.0;

scalar pi = 3.14159;
forAll(T,cI)
{
  vector pos_ = mesh.C()[cI];
  scalar r = Foam::sqrt(Foam::sqr((pos_.x())/r0) + Foam::sqr((pos_.y()-350.0)/r0));
  if(r<=1.00)
  {
  // Density 
  
  scalar thetacI = 0.25*(1.0+Foam::cos(pi*r))+theta0;
  scalar Par = 1.0 - (g_*pos_.y())/(Cp.value()*thetacI);
  scalar PgrRTheta = Pg/(R.value()*thetacI);

  rho[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  // Temperature 
  T[cI] = thetacI-g_*pos_.y()/Cp.value();
  thetaP[cI] = thetacI-theta0;
  }
  // The pressure is set into the class numericFlux  
}

T.correctBoundaryConditions();
thetaP.correctBoundaryConditions();

Info << "Done" << endl;

/*
 Initialize the 
 */
