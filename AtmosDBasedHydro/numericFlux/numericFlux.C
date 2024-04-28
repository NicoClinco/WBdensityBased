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

\*---------------------------------------------------------------------------*/

#include "numericFlux.H"
#include "MDLimiter.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Flux, class Limiter>
Foam::numericFlux<Flux, Limiter>::numericFlux
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    basicThermo& thermo
)
:
    numericFluxBase<Flux>(rho.mesh()),
    rho_(rho),
    U_(U),
    T_(T),
    thermo_(thermo),
    p_
    (
      IOobject
       (
	 "p",
	 this->mesh().time().timeName(),
         this->mesh(),
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
       ),
      ( rho_*(GEA::PerfectGasConstant::RDim)*T_)
    ),
    rhoFlux_
    (
        "phi",    
        (linearInterpolate(rho_*U_) & this->mesh().Sf())
    ),
    rhoUFlux_
    (
      "rhoUFlux",
      rhoFlux_*linearInterpolate(U_)
    ),
    rhoEFlux_
    (
        "rhoEFlux",  
        rhoFlux_*linearInterpolate(GEA::PerfectGasConstant::CvDim*T_ + 0.5*magSqr(U_))
    ),
   MomentumSource_
   (
       IOobject
       (
         "MomSource",
         this->mesh().time().timeName(),
         this->mesh(),
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
       ),  
       rho_*dimensionedVector("g",dimAcceleration,vector(0,-9.81,0))*scalar(0.0),
       Foam::zeroGradientFvPatchField<Foam::vector>::typeName
   ),
   EnergySource_
   (
      "Esource",
      rho_*U_ & (dimensionedVector("g",dimAcceleration,vector(0,-9.81,0)))
   ),
   Dp_
   (
    "Dp",
    p_*dimensionedVector("z",dimensionSet(0,-1,0,0,0,0,0),vector(0.0,0.0,0.0))
   ),
   Drho_
   (
    "Drho",
    rho_*dimensionedVector("z",dimensionSet(0,-1,0,0,0,0,0),vector(0.0,0.0,0.0))
   )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const unallocLabelList& owner = this->mesh().owner();
    const unallocLabelList& neighbour = this->mesh().neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = this->mesh().Sf();
    const surfaceScalarField& magSf = this->mesh().magSf();

    const volVectorField& cellCentre = this->mesh().C();
    const surfaceVectorField& faceCentre = this->mesh().Cf();
   

    const doubleScalar Cv = GEA::PerfectGasConstant::cv;
    const doubleScalar R  = GEA::PerfectGasConstant::R;

	/*
	const volVectorField& Drho = this->Drho_;
    const volVectorField& Dp   = this->Dp_;
	*/
	

    // Needed for limiter:
    const tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();
    // Standard Limiter for velocity:
    MDLimiter<vector, Limiter> vectorULimiter
    (
     this->U_,
     gradU
    );
    const volVectorField& ULimiter = vectorULimiter.phiLimiter();
    
    forAll (owner, faceI)
    {
        const label own = owner[faceI];          // own.
        const label nei = neighbour[faceI];      // nei.
	
		// Info << own << " "<< nei << "\n";
        const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];

		/*
		scalar deltapRight   = Dp[nei] & (deltaRRight);
		scalar deltapLeft    = Dp[own] & (deltaRLeft);
		scalar deltarhoRight = Drho[nei] & (deltaRRight);
		scalar deltarhoLeft  = Drho[own] & (deltaRLeft);
		*/
		
        
		doubleScalar rhoFNe = 0.0;
		doubleScalar rhoFOwn = 0.0;
		doubleScalar pFNe = 0.0;
		doubleScalar pFOwn = 0.0;
	
		// Owner to face
        GEA::HydrostaticReconstruction_
		  (
		   cellCentre[own],
		   faceCentre[faceI],
		   p_[own],
		   rho_[own],
		   pFOwn,
		   rhoFOwn
		   );

        // Neighbour to face
		GEA::HydrostaticReconstruction_
		  (
		   cellCentre[nei],
		   faceCentre[faceI],
		   p_[nei],
		   rho_[nei],
		   pFNe,
		   rhoFNe
		   );

	doubleScalar TFown = pFOwn/(R*rhoFOwn);
    doubleScalar TFnei = pFNe/(R*rhoFNe);

	// Evaluate the standard flux.

        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            pFOwn,
            pFNe,
            U_[own]+ cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
            U_[nei]+ cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
            TFown,
            TFnei,
            R,
            R,
            Cv,
            Cv,
            Sf[faceI],
            magSf[faceI]
        );
    } // end internal loop

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        
        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryField()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryField()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryField()[patchi];

        // Patch fields
        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const vectorField& pU = U_.boundaryField()[patchi];
        // const scalarField& pT = T_.boundaryField()[patchi];
	const labelUList& pFaceCells = p_.mesh().boundary()[patchi].faceCells();

       
        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
    	const volVectorField& cellCentre = this->mesh().C();
	const fvsPatchVectorField& pfC = faceCentre.boundaryField()[patchi];

        forAll (pp, facei)
        {
           doubleScalar rhoFNe = 0.0;
           doubleScalar rhoFOwn = 0.0;
           doubleScalar pFNe = 0.0;
           doubleScalar pFOwn = 0.0;

          GEA::HydrostaticReconstruction_
          (
           cellCentre[pFaceCells[facei]],
           pfC[facei],
           p_[pFaceCells[facei]],
	   rho_[pFaceCells[facei]],
           pFOwn,
           rhoFOwn
          );

          doubleScalar TFown = pFOwn/(R*rhoFOwn);
   
          Flux::evaluateFlux
          (
             pRhoFlux[facei],
             pRhoUFlux[facei],
             pRhoEFlux[facei],
             pFOwn,
             pFOwn,
             pU[facei],
             pU[facei],
             TFown,
             TFown,
             R,
             R,
             Cv,
             Cv,
             pSf[facei],
             pMagSf[facei]
          );
        }// end loop specific boundary
   }// end loop boundaries
}// end method

template<class Flux, class Limiter>
void Foam::numericFlux<Flux, Limiter>::EvalMomSource()
{
  //  Info << "Evaluating the momentum source\n";
  volVectorField& Mom = this->MomentumSource_;
  Mom *= scalar(0.0);
  // At every timestep we freeze the background-state
  const volScalarField& p0 = this->p_;
  const volScalarField& rho0 = this->rho_;
  GEA::EvalBalancedSource
  (
    Mom,
    p0,
    rho0
  );

}

template<class Flux,class Limiter>
void Foam::numericFlux<Flux, Limiter>::EvalEnergySource()
{
  volScalarField& Esource = this->EnergySource_;
  dimensionedVector g("g",dimAcceleration,vector(0,-9.81,0));
  Esource.internalField() = (rho_*(g & U_));
  // Esource.correctBoundaryConditions(); // ZeroGradient;
}


template<class Flux,class Limiter>
void Foam::numericFlux<Flux, Limiter>::EvaluateLimiter()
{
  const volScalarField& rho = this->rho_;
  const volScalarField&  p = this->p_;
  volVectorField& Drho = this->Drho_;
  volVectorField& Dp   = this->Dp_;
  
  Drho*=scalar(0.0);
  Dp*=scalar(0.0);
 /*  
  GEA::ConstructLimiter
  (
    Drho,
    Dp,
    rho,
    p
  );
 */
}

// ************************************************************************* //
