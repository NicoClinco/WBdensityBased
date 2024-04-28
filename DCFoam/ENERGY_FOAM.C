/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    NicoFoam

Description
    Density-based compressible explicit time-marching for
    atmospheric flows

Author
    Nicola Clinco

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "hllcFlux.H"
#include "roeFlux.H"
#include "rusanovFlux.H"
#include "betaFlux.H"
#include "AusmUPFlux.H"
#include "MDLimiter.H"
#include "firstOrderLimiter.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"
#include "numericFlux.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createTimeControls.H"

//ofstream FilePostMomentum;
//FilePostMomentum.open("FilePostMomentum.txt",std::ios::app);

ofstream FilePostVelocity;
FilePostVelocity.open("FilePostVel.txt",std::ios::app);
FilePostVelocity << "t" << "\t" << "Uymax" << "\t" << "Source-flux" << "\n";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // Runge-Kutta coefficient
    scalarList beta(4);
    beta[0] = 0.1100;
    beta[1] = 0.2766;
    beta[2] = 0.5000;
    beta[3] = 1.0000;
    dimensionedVector gDim("g",dimAcceleration,vector(0,9.81,0));
    // Switch off solver messages
    lduMatrix::debug = 0;

    while (runTime.run())
    {
#       include "readTimeControls.H"
   // #       include "readFieldBounds.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "\n Time = " << runTime.value() << endl;

        // Low storage Runge-Kutta time integration
        forAll (beta, i)
        {
            // Solve the approximate Riemann problem for this time step
	    // Evaluate the source term

	    dbnsFlux.EvalMomSource();
	    dbnsFlux.EvalEnergySource();
	    dbnsFlux.EvaluateLimiter();

            dbnsFlux.computeFlux();

	    tmp<volScalarField> tEnerSource = -rho*(U&gDim);
	    tmp<volScalarField> tEn = Cv*T;
	    tEn().boundaryField() = Cv.value()*T.boundaryField();
	    tmp<volScalarField> tLapEn = fvc::laplacian(mu,tEn()+0.5*magSqr(U));
            
	     

	    FilePostVelocity << runTime.value() << "\t" << Foam::max(mag(U.component(0))).value() << "\t" << 
	      Foam::max(mag(U.component(1))).value() <<"\t" <<  "\n";

            // Time integration
	     
            solve
            (
              1.0/beta[i]*fvm::ddt(rho)
              + fvc::div(dbnsFlux.rhoFlux())
            );
	    
            solve
            (
              1.0/beta[i]*fvm::ddt(rhoU)
              + fvc::div(dbnsFlux.rhoUFlux()) == dbnsFlux.MomentumSource() + fvc::laplacian(mu,U) 
            );
	     
            solve
            (
              1.0/beta[i]*fvm::ddt(rhoE)
              + fvc::div(dbnsFlux.rhoEFlux())  == tLapEn() + tEnerSource()
            );

            #include "updateFields.H"
	   
        }
	
        runTime.write();
	
        Info<< "    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
