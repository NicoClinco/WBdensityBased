scalar theta0(300.0);
scalar Pg(100000.0);
scalar g_(9.81);
scalar CpOR = Cp.value()/R.value();
scalar CvOR = Cv.value()/R.value();
scalar CpTheta = Cp.value()*theta0;
scalar RTheta = R.value()*theta0;
scalar PgoRTheta= Pg/RTheta;
vectorField VV = mesh.C();
//const surfaceVectorField& faceCentre = mesh.Cf();

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
  // calculated from the density and temperature
}

scalar rx = 4000.0;
scalar ry = 2000.0;
scalar pi = 3.14159;
forAll(T,cI)
{
  vector pos_ = mesh.C()[cI];
  scalar r = Foam::sqrt(Foam::sqr(pos_.x()/rx) + Foam::sqr((pos_.y()-3000)/ry));
  if(r<=1.001)
  {
  // Density 
  scalar thetacI = theta0-7.5*(1.0+Foam::cos(pi*r));
  scalar Par = 1.0 - (g_*pos_.y())/(Cp.value()*thetacI);
  scalar PgrRTheta = Pg/(R.value()*thetacI);

  rho[cI] = PgrRTheta * Foam::pow(Par,CvOR);
  // Temperature 
  T[cI] = thetacI-g_*pos_.y()/Cp.value();
  }
  // The pressure is set into the class numericFlux,
  // calculated from pressure and temperature
}

//Loop trough the left boundary:
/*
label leftBoundary = mesh.boundaryMesh().findPatchID("left");
fvPatchScalarField& tP = T.boundaryField()[leftBoundary];
const fvsPatchVectorField& pfC = faceCentre.boundaryField()[leftBoundary];

forAll(tP, facei)
{
  //Temperature
  vector pos_ = pfC[facei];
  scalar r = Foam::sqrt(Foam::sqr(pos_.x()/rx) + Foam::sqr((pos_.y()-3000)/ry));
  if(r<=1.001)
  {
    scalar thetacI = theta0-7.5*(1.0+Foam::cos(pi*r));
    T[facei] = thetacI-g_*pos_.y()/Cp.value();
  }
}
*/

Info << "Done" << endl;
