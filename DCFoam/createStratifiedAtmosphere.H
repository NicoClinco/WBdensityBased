Info << "Initilializing the Background State \n"<< endl;

using namespace GEA::PerfectGasConstant;

doubleScalar Oogamma1 = 1/(gamma-1);
doubleScalar Ookappa      = 1/kappa; //g/(g-1)
doubleScalar Oogamma   = 1/gamma;
doubleScalar rhoHydro0 = pHydro0/(R.value()*theta0);
scalar g_(9.81000);
vectorField VV = mesh.C();

// Brunt-Vaisaala:
scalar N0=0.01;
scalar factor = kappa*g_*g_/(R.value()*theta0*N0*N0);
scalar expfactor = -N0*N0/(g_);


forAll(p0,cI)
{
  scalar posZ = VV[cI].y();
  doubleScalar thetaz = theta0*Foam::exp(-expfactor*posZ);
  p0[cI] = pHydro0*Foam::exp(-g_*posZ/(R.value()*theta0));
  doubleScalar pOpg = p0[cI]/pHydro0;
  rho[cI] = Foam::pow(pOpg,Oogamma)*(pHydro0/(R.value()*thetaz));
  T[cI] = p0[cI]/(R.value()*rho[cI]);
}


Info << "Done" << endl;
