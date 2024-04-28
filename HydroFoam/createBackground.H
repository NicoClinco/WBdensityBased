Info << "Initilializing the Background State \n"<< endl;

using namespace GEA::PerfectGasConstant;

doubleScalar Oogamma1 = 1/(gamma-1);
doubleScalar Ookappa      = 1/kappa;
doubleScalar Oogamma   = 1/gamma;
doubleScalar rhoHydro0 = pHydro0/(R.value()*theta0);
scalar g_(9.81000);
vectorField VV = mesh.C();


forAll(p0,cI)
{
  scalar posZ = VV[cI].y();
  doubleScalar inner = (pHydro0 - kappa*g_*rhoHydro0*posZ);
  doubleScalar first =Foam::pow(pHydro0,-Oogamma1);
  p0[cI]    = (first*Foam::pow(inner,Ookappa));
  doubleScalar pOpg = p0[cI]/pHydro0;
  rho[cI] = rhoHydro0*Foam::pow(pOpg,Oogamma);
  T[cI] = p0[cI]/(R.value()*rho[cI]);
}

Info << "Done" << endl;
