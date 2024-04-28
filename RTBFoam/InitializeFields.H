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
/*
scalar rx = 4000;
scalar ry = 2000;
scalar pi = 3.14;
forAll(T,cI)
{
  vector pos_ = mesh.C()[cI];
  scalar r = Foam::sqrt(Foam::sqr(pos_.x()/rx) + Foam::sqr((pos_.y()-3000)/ry));
  if(r<=1.0)
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
*/
Info << "Done" << endl;
