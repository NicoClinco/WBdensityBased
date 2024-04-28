// Correct the bcs for the pressure and density based on a hydrostatic
// reconstruction:
void CorrectHydroBCs(volScalarField& _p_, volScalarField& _rho_)
{
   const volVectorField& cellCentre = _p_.mesh().C();
   const surfaceVectorField& faceCentre = _p_.mesh().Cf();
   const surfaceVectorField& Sf = _p_.mesh().Sf();
   const surfaceScalarField& magSf = _p_.mesh().magSf();
   
  forAll(_p_.boundaryField(),patchi)
  {
    // Corresponding GeometricalPfields
    const labelUList& pFaceCells = _p_.mesh().boundary()[patchi].faceCells();
    const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
    const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
    const fvsPatchVectorField& pfC = faceCentre.boundaryField()[patchi];
    
    
    // Corresponding PhysicalPfields
    fvPatchScalarField& pp = _p_.boundaryField()[patchi];
    fvPatchScalarField& prho = _rho_.boundaryField()[patchi];

    forAll(pp,facei)
    {
      // Reconstruct the pressure and density:
          GEA::HydrostaticReconstruction_
          (
           cellCentre[pFaceCells[facei]],
           pfC[facei],
           _p_[pFaceCells[facei]],
	   _rho_[pFaceCells[facei]],
           pp[facei],
           prho[facei]
          );
    }
    
   }
}
