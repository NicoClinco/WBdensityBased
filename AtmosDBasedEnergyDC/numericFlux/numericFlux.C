

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
      ( rho_*(dimensionedScalar("R",dimGasConstant,287.0))*T_)
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
        rhoFlux_*linearInterpolate(dimensionedScalar("Cv",dimGasConstant,720.0)*T_ + 0.5*magSqr(U_))
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
    
    const volVectorField& Drho = this->Drho_;
    const volVectorField& Dp   = this->Dp_;
    // const volTensorField& DV   = this->DV_;

    const scalar Cv = 720.0;
    const scalar R  = 287.0;

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

        const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
        const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];
        
	scalar deltapRight   = Dp[nei] & (deltaRRight);
	scalar deltapLeft    = Dp[own] & (deltaRLeft);
	scalar deltarhoRight = Drho[nei] & (deltaRRight);
	scalar deltarhoLeft  = Drho[own] & (deltaRLeft);

	scalar rhoFNe = 0.0;
	scalar rhoFOwn = 0.0;
	scalar pFNe = 0.0;
	scalar pFOwn = 0.0;
	
	// Construct Local equilibrum for the face:
	
	GEA::ConstructLocalEquilibrum
	(
	  faceCentre[faceI].y(),
	  cellCentre[nei].y(),
	  cellCentre[own].y(),
	  rho_[own], 
	  rho_[nei],
	  p_[own],
	  p_[nei],
	  rhoFNe,
	  rhoFOwn,
	  pFNe,
	  pFOwn
	);
	
	 pFOwn += deltapLeft;
	 pFNe  += deltapRight;
	 rhoFOwn += deltarhoLeft;
         rhoFNe  += deltarhoRight;
	scalar TFown = pFOwn/(R*rhoFOwn);
	scalar TFnei = pFNe/(R*rhoFNe);

	// Evaluate the standard flux.

        Flux::evaluateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            pFOwn,
            pFNe,
            U_[own] + cmptMultiply(ULimiter[own], (deltaRLeft & gradU[own])),
            U_[nei] + cmptMultiply(ULimiter[nei], (deltaRRight & gradU[nei])),
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
        //const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

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
           scalar rhoFNe = 0.0;
           scalar rhoFOwn = 0.0;
           scalar pFNe = 0.0;
           scalar pFOwn = 0.0;

           // Construct Local equilibrum for the face:

           GEA::ConstructLocalEquilibrum
           (
            pfC[facei].y(),
            cellCentre[pFaceCells[facei]].y(),
            cellCentre[pFaceCells[facei]].y(),
            rho_[pFaceCells[facei]],
            rho_[pFaceCells[facei]],
            p_[pFaceCells[facei]],
            p_[pFaceCells[facei]],
            rhoFNe,
            rhoFOwn,
            pFNe,
            pFOwn
           );

          scalar TFown = pFOwn/(R*rhoFOwn);
          scalar TFnei = pFNe/(R*rhoFNe);
          

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
             TFnei,
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
  GEA::ConstructLimiter
  (
    Drho,
    Dp,
    rho,
    p
  );

}

// ************************************************************************* //
