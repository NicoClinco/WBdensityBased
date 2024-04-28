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

Class
    numericFlux

Description
    Single level numeric flux class for density-based solvers

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

SourceFiles
    numericFlux.H
    numericFlux.C

\*---------------------------------------------------------------------------*/

#ifndef numericFlux_H
#define numericFlux_H

#include "numericFluxBase.H"
#include "basicThermo.H"
#include "fvc.H"
#include "IntegrateUtilitiesNew.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class numericFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Flux, class Limiter>
class numericFlux
:
    public numericFluxBase<Flux>
{
    // Private data

        // Reference to primitive fields
            
	    //- density
            const volScalarField& rho_;  
            //- Velocity
            const volVectorField& U_;

            //- Static temperature
            const volScalarField& T_;

	    basicThermo& thermo_;

	    // - pressure
            volScalarField p_; // total pressure
	
            surfaceScalarField rhoFlux_;
	    
	    surfaceVectorField rhoUFlux_;

	    surfaceScalarField rhoEFlux_;

	    // MomentumSource
	    volVectorField MomentumSource_;
            // EnergySource
            volScalarField EnergySource_;
	    
	    // Limiters
	    volVectorField Dp_;

	    volVectorField Drho_; 

	    
	    //- Hydrostatic quantities at the cell center

	     autoPtr<volScalarField> prho0_;  
	    
	     autoPtr<volScalarField> pp0_; 
	   
	     autoPtr<volScalarField> pT0_;  


    // Private Member Functions

        //- Disallow default bitwise copy construct
        numericFlux(const numericFlux&);

        //- Disallow default bitwise assignment
        void operator=(const numericFlux&);


        //- Return internal field of mass flux
        const scalarField& rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }

        //- Return access to internal field of mass flux
        scalarField& rhoFluxI()
        {
            return rhoFlux_.internalField();
        }

        //- Return internal field of momentum flux
        const vectorField& rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }

        //- Return access to internal field of momentum flux
        vectorField& rhoUFluxI()
        {
            return rhoUFlux_.internalField();
        }

         //- Return access to internal field of energy flux
        const scalarField& rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }

        //- Return access to internal field of energy flux
        scalarField& rhoEFluxI()
        {
            return rhoEFlux_.internalField();
        }


public:

    //- Runtime type information
    TypeName("numericFlux");

    // Constructors

        //- Construct from components
        

	numericFlux
	(
	   const volScalarField& rho,
	   const volVectorField& U,
	   const volScalarField& T,
	   basicThermo& thermo
	);


    //- Destructor
    virtual ~numericFlux()
    {}

    void InitializeHydroStaticField
    (
     const volScalarField& rho0,
     const volScalarField& p0,
     const volScalarField& T0
    )
    {
	prho0_ = autoPtr<volScalarField> 
		(
		  new volScalarField("rho0",rho0)
		);
	pp0_ = autoPtr<volScalarField>
		(
		 new volScalarField("p0",p0)
		);
	pT0_ = autoPtr<volScalarField>
		(
		 new volScalarField("T0",T0)
		);
	// Initialize the prho0_, pp0_,pT0_;
    }
    

        // Return fluxes

            //- Return density flux
            virtual const surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            //- Return velocity flux
            virtual const surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }

            //- Return energy flux
            virtual const surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }
	    
	    const volScalarField& p() const
	    {
		return p_;
	    }
	    volScalarField& p()
	    {
	        return p_; 
	    }

	    const volVectorField& MomentumSource() const
	    {
		return MomentumSource_;
	    }
	    
	    volVectorField& MomentumSourceRef()
            {
            	return MomentumSource_;
            }

	    const volScalarField& EnergySource() const
	    {
		return EnergySource_;
	    }

	    volScalarField& EnergySource()
	    {
		return EnergySource_;
	    }



       // Return residuals

            //- Return density equation residual
            virtual tmp<scalarField> rhoResidual() const
            {
                return fvc::div(rhoFlux_)().internalField();
            }

            //- Return momentum equation flux
            virtual tmp<vectorField> rhoUResidual() const
            {
                return fvc::div(rhoUFlux_)().internalField();
            }

            //- Return energy equation flux
            virtual tmp<scalarField> rhoEResidual() const
            {
                return fvc::div(rhoEFlux_)().internalField();
            }
	

        // Update fluxes based on current state

            //- Compute flux
            virtual void computeFlux();

	    //- Compute Momentum Source
	    void EvalMomSource();
	    //- Compute the energy Source
	    void EvalEnergySource();
            //- Evaluate the limiter for pressure and density
	    void EvaluateLimiter();
	    //- Evaluate the limiter for the velocity
	    //void EvaluateLimiterV();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "numericFlux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
