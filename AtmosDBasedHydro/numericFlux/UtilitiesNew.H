namespace GEA
{

/* -------------------end Local Limiters ---------------- */

// Based on the papaer of Kapelli :
// Extrapolate the density and the pressure 
// WITHOUT assuming \theta = 300K  
inline void HydrostaticReconstruction_
(
 const vector& posI,
 const vector& posFi,
 const doubleScalar& p0,
 const doubleScalar& rho0,
 doubleScalar& pf_,
 doubleScalar& rhof_
)
{
   doubleScalar gamma = GEA::PerfectGasConstant::gamma;
   scalar gamma1 = gamma - 1.0;
   scalar gamma2 = 1/gamma1;
   scalar gamma3 = gamma1/gamma;
   scalar g = 9.81;

   doubleScalar& pf = pf_;
   doubleScalar& rhof = rhof_;
   doubleScalar K = p0/(Foam::pow(rho0,gamma));
   doubleScalar rho0g = Foam::pow(rho0,gamma1);
   rhof = Foam::pow((rho0g + g*(gamma3/K)*(posI.y()-posFi.y())),gamma2);
   pf   = K * (Foam::pow(rhof,gamma));
}

  
void EvalBalancedSource
(
   volVectorField& MomentumSource,
   const volScalarField& p0,
   const volScalarField& rho0
)
{
    const fvMesh& mesh = p0.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceVectorField& faceCentre = mesh.Cf();
    const volVectorField& cellCentre = mesh.C();
    const volScalarField& p0_ = p0;
    const volScalarField& rho0_ = rho0;
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    
    forAll(owner, facei)
    {
      vector posFace = faceCentre[facei];
      vector posOwn = cellCentre[owner[facei]];
      vector posNe   = cellCentre[neighbour[facei]]; 
      vector faceVal = Sf[facei];

      //Info << "owner: "<< owner << "\n";
           
      // OwnerToFaceReconstruction
      doubleScalar p0Own = p0_[owner[facei]];
      doubleScalar rho0Own = rho0_[owner[facei]];
      doubleScalar pFown = 0.0;
      doubleScalar rhoFown = 0.0; 
      HydrostaticReconstruction_(posOwn,posFace,p0Own,rho0Own,pFown,rhoFown);

      //NeighbourToFaceReconstruction
      doubleScalar p0Ne = p0_[neighbour[facei]];
      doubleScalar rho0Ne = rho0_[neighbour[facei]];
      doubleScalar pFNe = 0.0;
      doubleScalar rhoFNe = 0.0;
      HydrostaticReconstruction_(posNe,posFace,p0Ne,rho0Ne,pFNe,rhoFNe);

      MomentumSource[owner[facei]]+= (pFown*faceVal);
      MomentumSource[neighbour[facei]]-= (pFNe*faceVal); 
    }//end internal for

   forAll(mesh.boundary(),patchi)
   {
      const unallocLabelList& pFaceCells =
          mesh.boundary()[patchi].faceCells();

      const fvsPatchField<vector>& pssf = Sf.boundaryField()[patchi];
      const fvsPatchField<vector>& pfp = faceCentre.boundaryField()[patchi];
      const fvPatchField<scalar>& pp = p0_.boundaryField()[patchi];

      forAll(mesh.boundary()[patchi], facei){
	   vector faceVal = pssf[facei];
	   
           vector posFace = pfp[facei];
	   vector posOwn = cellCentre[pFaceCells[facei]]; 
           
           doubleScalar p0Own = p0_[pFaceCells[facei]];
	   doubleScalar rho0Own = rho0_[pFaceCells[facei]];
	   doubleScalar pFown = 0.0;
	   doubleScalar rhoFown = 0.0;
	   HydrostaticReconstruction_(posOwn,posFace,p0Own,rho0Own,pFown,rhoFown); 
	   
           MomentumSource[pFaceCells[facei]]+= (pFown*faceVal);
       }// end specific boundary
    }// end boundaries
  MomentumSource.internalField()/=mesh.V();
  MomentumSource.correctBoundaryConditions();

}//end reconstruction

}//end namespace GEA
