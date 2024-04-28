namespace GEA
{
  namespace PerfectGasConstant
  {
    doubleScalar  g = 9.810;
    doubleScalar R = 287.00;
    doubleScalar cp = 1007.00;
    doubleScalar cv = 720.00;
    doubleScalar theta0 = 300.00;
    doubleScalar gamma = 1.39861111;
    doubleScalar kappa = (gamma-1)/gamma;
    doubleScalar pHydro0 = 100000.0;
    doubleScalar Rtheta0 = R*theta0;
    doubleScalar RoCv      = R/cv;
    // Used for c0:
    doubleScalar c0 = (Foam::pow(R,gamma)/Foam::pow(pHydro0,RoCv));

    // Dimensioned-Variables
    
    dimensionedScalar gDim("gDim",dimAcceleration,g);
    dimensionedScalar RDim("RDim",dimGasConstant,R);
    dimensionedScalar CvDim("CvDim",dimGasConstant,cv);
    dimensionedScalar CpDim("CvDim",dimGasConstant,cp);
  }
}
