#include "fvCFD.H"
namespace GEA
{ 
  
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

    scalar gamma = 1.3986;
    scalar gamma1 = gamma - 1.0;
    scalar gamma2 = 1/gamma1;
    scalar gamma3 = gamma1/gamma;
    scalar g = 9.81;


    forAll(owner, facei)
    {
	scalar yFace = faceCentre[facei].y();  // y-position of the face
	vector faceVal = Sf[facei];          // vector
       
	//OWNER
	scalar yOwnCenter = cellCentre[owner[facei]].y();
	scalar rhoOwn = rho0_[owner[facei]];
	scalar pOwn = p0_[owner[facei]];
	scalar Kown = pOwn/(Foam::pow(rhoOwn,gamma));
  	scalar rho0Owng = Foam::pow(rhoOwn,gamma1);
  	scalar rhoOwnF_ = Foam::pow((rho0Owng + g*(gamma3/Kown)*(yOwnCenter-yFace)),gamma2);
  	scalar pOwnF_   = Kown * (Foam::pow(rhoOwnF_,gamma));
	
	//NEIGHBOUR
	scalar yNeCenter = cellCentre[neighbour[facei]].y();
	scalar rhoNe = rho0_[neighbour[facei]];
	scalar pNe = p0_[neighbour[facei]];
  	scalar KNe = pNe/(Foam::pow(rhoNe,gamma));
  	scalar rho0Neg = Foam::pow(rhoNe,gamma1);
  	scalar rhoNeF_ =Foam::pow((rho0Neg + g*(gamma3/KNe)*(yNeCenter-yFace)),gamma2);
 	scalar pNeF_   = KNe * Foam::pow(rhoNeF_,gamma);

	MomentumSource[owner[facei]]+= (pOwnF_*faceVal);
	MomentumSource[neighbour[facei]]-= (pNeF_*faceVal);


   }// end internal Loop
   
   forAll(mesh.boundary(),patchi)
   {
      const unallocLabelList& pFaceCells =
          mesh.boundary()[patchi].faceCells();
    
      const fvsPatchField<vector>& pssf = Sf.boundaryField()[patchi];
      const fvsPatchField<vector>& pfp = faceCentre.boundaryField()[patchi];
      forAll(mesh.boundary()[patchi], facei){
	   scalar yFace = pfp[facei].y();
	   scalar yOwnCenter  = cellCentre[pFaceCells[facei]].y();
	   scalar rhoOwn      = rho0_[pFaceCells[facei]];
	   scalar pOwn        = p0_[pFaceCells[facei]];
	   scalar Kown = pOwn/(Foam::pow(rhoOwn,gamma));
	   scalar rho0Owng = Foam::pow(rhoOwn,gamma1);
	   scalar rhoOwnF_ = Foam::pow((rho0Owng + g*(gamma3/Kown)*(yOwnCenter-yFace)),gamma2);
           scalar pOwnF_   = Kown * (Foam::pow(rhoOwnF_,gamma));

	   vector faceVal = pssf[facei]; 
	   MomentumSource[pFaceCells[facei]]+= (pOwnF_*faceVal);
       }// end specific boundary
    }// end boundaries
  MomentumSource.internalField()/=mesh.V();
  MomentumSource.correctBoundaryConditions(); 
  }// End function
  
   void ConstructLocalEquilibrum
   ( 
     scalar yF,
     scalar yNeC,
     scalar yownC,
     scalar rho0Own,
     scalar rho0Ne,
     scalar p0Own,
     scalar p0Ne,
     scalar& rhoNeF,
     scalar& rhoOwnF,
     scalar& pNeF,
     scalar& pOwnF
   )
   {
   scalar rho0Own_  = rho0Own;
   scalar rho0Ne_ = rho0Ne;
   scalar p0Own_ = p0Own;
   scalar p0Ne_ = p0Ne;
   scalar& rhoOwnF_ = rhoOwnF;
   scalar& rhoNeF_  = rhoNeF;
   scalar& pOwnF_ = pOwnF;
   scalar& pNeF_ = pNeF;
   scalar yFace = yF;
   scalar yNeCenter = yNeC;
   scalar yOwnCenter = yownC;
   scalar gamma = 1.3986;
   scalar gamma1 = gamma - 1.0;
   scalar gamma2 = 1/gamma1;
   scalar gamma3 = gamma1/gamma;
   scalar g = 9.81;

   // OWNER //
   scalar Kown = p0Own_/(Foam::pow(rho0Own_,gamma));
   scalar rho0Owng = Foam::pow(rho0Own_,gamma1);
   rhoOwnF_ = Foam::pow((rho0Owng + g*(gamma3/Kown)*(yOwnCenter-yFace)),gamma2);
   pOwnF_   = Kown * (Foam::pow(rhoOwnF_,gamma));
   // END

   // NEIGHBOUR // 
   scalar KNe = p0Ne_/(Foam::pow(rho0Ne_,gamma));
   scalar rho0Neg = Foam::pow(rho0Ne_,gamma1);
   rhoNeF_ =Foam::pow((rho0Neg + g*(gamma3/KNe)*(yNeCenter-yFace)),gamma2);
   pNeF_   = KNe *(Foam::pow(rhoNeF_,gamma));
   }
/* ***************************************************** */ 
  void ExtrapolateEquilibrum
  (
   const scalar& p0,
   const scalar& rho0,
   scalar y0,
   scalar y,
   scalar& pE_,
   scalar& rhoE_
  )
  {
   // Recostruction from "0" to a the generic state pE,rhoE in y
   // p/rho^gamma = cost.
   
   scalar gamma = 1.3986;
   scalar gamma1 = gamma - 1.0;
   scalar gamma2 = 1/gamma1;
   scalar gamma3 = gamma1/gamma;
   scalar g = 9.81;
   
   scalar& pE = pE_;
   scalar& rhoE = rhoE_;
   scalar K = p0/(Foam::pow(rho0,gamma));
  
   scalar rho0g = Foam::pow(rho0,gamma1);
                                             
   rhoE = Foam::pow((rho0g + g*(gamma3/K)*(y0-y)),gamma2);
   pE   = K * (Foam::pow(rhoE,gamma));
  }// end ExtrapolateEquilibrum
/* ******************************************************* */

// For vectors
void CalculateLocalLimiter
 (
  vector rI,
  vector rIp,
  vector rIm,
  const vector& pertI,
  const vector& pertIp,
  const vector& pertIm,
  vector& Dp_
 )
 {
   vector& Dp = Dp_;
   scalar deltarIm = Foam::mag(rI - rIm);
   scalar deltarIp = Foam::mag(rIp - rI);
   vector delta_m = pertI - pertIp;
   vector delta_p = pertIp - pertI;

   vector delta_m_dr = delta_m/deltarIm;
   vector delta_p_dr = delta_p/deltarIp;
   vector Abs_dmdr = Foam::cmptMag(delta_m_dr);
   vector Abs_dpdr = Foam::cmptMag(delta_p_dr);

   vector signA(0.0,0.0,0.0);
   vector signB(0.0,0.0,0.0);
   // check if a component is zero.
   for(int i=0; i<3;i++){
     if(Abs_dmdr[i]<1e-13){signA[i] = +1;}
     else{signA[i] = Foam::sign(delta_m_dr[i]);}

     if(Abs_dpdr[i]<1e-13){signB[i] = +1;}
     else{signB[i] = Foam::sign(delta_p_dr[i]);}
   }

   vector minium = Foam::min(Abs_dmdr,Abs_dpdr);
   Dp = 0.5*Foam::cmptMultiply((signA + signB),Foam::min(2*minium,Foam::cmptMag(0.5*(delta_m_dr+delta_p_dr))));

}// end function



/* ******************************************************* */

  // For scalars
  void CalculateLocalLimiter
  (
   vector rI,
   vector rIp,
   vector rIm,
   const scalar& q_pertI,
   const scalar& q_pertIp,
   const scalar& q_pertIm,
   scalar& Dq_
  )
  {
    scalar& Dq = Dq_;
    scalar deltarIm = Foam::mag(rI - rIm);
    scalar deltarIp = Foam::mag(rIp - rI);
    scalar delta_q_m  = q_pertI - q_pertIm; //qp(xi)-qp(xi-1)
    scalar delta_q_p  = q_pertIp - q_pertI; //qp(xi+1)-qp(xi)
    
    scalar dq_m_dr = delta_q_m/deltarIm;
    scalar dq_p_dr = delta_q_p/deltarIp;
    
   
    scalar signA = Foam::sign(delta_q_m);
    scalar signB = Foam::sign(delta_q_p);
    
    // Generalized minmod : 
    //Dq = 0.5*(signA + signB)*Foam::min(Foam::mag(dq_m_dr),Foam::mag(dq_p_dr));
    
    scalar minium = Foam::min(Foam::mag(dq_m_dr),Foam::mag(dq_p_dr));
    Dq = 0.5*(signA + signB)*Foam::min(2*minium,Foam::mag(0.5*(dq_m_dr+dq_p_dr)));
     
  }


/* *********************************************************** */  
 
  // ConstructLimiter for the density and pressure:
  void ConstructLimiter
  (
   volVectorField& Drho_,
   volVectorField& Dp_,
   const volScalarField& rho_,
   const volScalarField& p_
  )
  {
   const fvMesh& mesh = p_.mesh();
   
   const volScalarField& rho(rho_);
   const volScalarField& p(p_);

   volVectorField& Drho(Drho_);
   volVectorField& Dp(Dp_);
     
   const volVectorField& cellCentre = mesh.C();
   scalar Lx = 25600.0; // Density-current
   scalar Ly = 6400.0;  //
   scalar LxLy = (Lx*Ly);
   label Ncells = cellCentre.size();
   
   scalar delta = Foam::sqrt(LxLy/(Ncells));
   scalar Nx_ = Lx/delta;
   scalar Ny_ = Ly/delta;
   int Nx = int(Nx_);
   int Ny = int(Ny_);

   // internal loop cycle for reconstruction 
   
   for(int Ix=1; Ix<Nx-1;Ix++)
   {
    for(int Iy=1; Iy<Ny-1;Iy++)
    {
      label cI = Ix + Iy*Nx ; // Enumeration of OF 
      label cInord = cI + Nx; 
      label cIeast = cI + 1;
      label cIwest = cI - 1;
      label cIsouth = cI - Nx; 

      scalar p_0 = p[cI]; 
      scalar rho_0 = rho[cI];
      scalar p_n = p[cInord];      scalar p0_n   = 0.0;  
      scalar rho_n = rho[cInord];  scalar rho0_n = 0.0;
      scalar p_e = p[cIeast];      scalar p0_e   = 0.0;
      scalar rho_e = rho[cIeast];  scalar rho0_e = 0.0; 
      scalar p_w = p[cIwest];      scalar p0_w   = 0.0;
      scalar rho_w = rho[cIwest];  scalar rho0_w = 0.0;
      scalar p_s = p[cIsouth];     scalar p0_s   = 0.0;
      scalar rho_s = rho[cIsouth]; scalar rho0_s = 0.0;
      scalar yC = cellCentre[cI].y();
        
      scalar pPertcI = 0.0;
      scalar rhoPertcI = 0.0;
      

        //------- Nord ------ //
	
        GEA::ExtrapolateEquilibrum
	(
	 p_0,
	 rho_0,
         yC,
         cellCentre[cInord].y(),
         p0_n,
	 rho0_n
	);

       // Step 2: Evaluate the perturbation in Nord
       scalar pPertNord = p_n - p0_n;
       scalar rhoPertNord = rho_n - rho0_n;  
       // ----------------- //

      //------- South------ //

        GEA::ExtrapolateEquilibrum
        (
         p_0,
         rho_0,
         yC,
         cellCentre[cIsouth].y(),
         p0_s,
         rho0_s
        );
	 
       scalar pPertSouth = p_s - p0_s;
       scalar rhoPertSouth = rho_s - rho0_s;  
      // ----------------- //

      // ----- East ------ //
      
        GEA::ExtrapolateEquilibrum
        (
         p_0,
         rho_0,
         yC,
         cellCentre[cIeast].y(),
         p0_e,
         rho0_e
        );

       scalar pPertEast = p_e - p0_e;
       scalar rhoPertEast = rho_e - rho0_e; 
      // ------------------ //

      // ----- West ------- //

        GEA::ExtrapolateEquilibrum
        (
         p_0,
         rho_0,
         yC,
         cellCentre[cIwest].y(),
         p0_w,
         rho0_w
        );

       scalar pPertWest = p_w - p0_w;
       scalar rhoPertWest = rho_w - rho0_w; 
     // ------------------- //
     
    // Evaluate Limiters:
    
    scalar Dpy = 0.0; 
    // y-LIMITER-pressure
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cInord],
     cellCentre[cIsouth],
     pPertcI,
     pPertNord,
     pPertSouth,
     Dpy
    );
   

    scalar Drhoy = 0.0;
    // y-LIMITER-density
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cInord],
     cellCentre[cIsouth],
     rhoPertcI,
     rhoPertNord,
     rhoPertSouth,
     Drhoy
    );
    
    scalar Dpx = 0.0;
    // x-LIMITER-pressure
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cIeast],
     cellCentre[cIwest],
     pPertcI,
     pPertEast,
     pPertWest,
     Dpx
    );

    
    scalar Drhox = 0.0;
    // x-LIMITER-pressure
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cIeast],
     cellCentre[cIwest],
     rhoPertcI,
     rhoPertEast,
     rhoPertWest,
     Drhox
    );

    Drho[cI].component(0) = Drhox;
    Dp[cI].component(0)   = Dpx;
    Dp[cI].component(1)   = Dpy;
    Drho[cI].component(1) = Drhoy;
    
    }//end for x
  }// end for y

  // Loop for the boundary : 
  int Iy = 0;
  for(int Ix=1; Ix<Nx-1;Ix++)
  {
   label cI = Ix;
   label cIeast = cI + 1;
   label cIwest = cI - 1;
   scalar p_0 = p[cI];
   scalar rho_0 = rho[cI];
   scalar p_e = p[cIeast];      scalar p0_e   = 0.0;
   scalar rho_e = rho[cIeast];  scalar rho0_e = 0.0;
   scalar p_w = p[cIwest];      scalar p0_w   = 0.0;
   scalar rho_w = rho[cIwest];  scalar rho0_w = 0.0;
   scalar yC = cellCentre[cI].y();
   scalar pPertcI = 0.0;
   scalar rhoPertcI = 0.0;
   
   GEA::ExtrapolateEquilibrum
   (
         p_0,
         rho_0,
         yC,
         cellCentre[cIeast].y(),
         p0_e,
         rho0_e
    );

    scalar pPertEast = p_e - p0_e;
    scalar rhoPertEast = rho_e - rho0_e;
   // ------------------ //

   // ----- West ------- //

   GEA::ExtrapolateEquilibrum
   (
         p_0,
         rho_0,
         yC,
         cellCentre[cIwest].y(),
         p0_w,
         rho0_w
    );

    scalar pPertWest = p_w - p0_w;
    scalar rhoPertWest = rho_w - rho0_w;
    

    scalar Dpx = 0.0;
    // x-LIMITER-pressure
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cIeast],
     cellCentre[cIwest],
     pPertcI,
     pPertEast,
     pPertWest,
     Dpx
    );


    scalar Drhox = 0.0;
    // x-LIMITER-pressure
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cIeast],
     cellCentre[cIwest],
     rhoPertcI,
     rhoPertEast,
     rhoPertWest,
     Drhox
    );

    Drho[cI].component(0) = Drhox;
    Dp[cI].component(0)   = Dpx;
 
  }// end bottom

  int Ix = 0;
  for(int Iy=1; Iy<Ny-1;Iy++)
  {
   label cI = Ix + Iy*Nx ;
   label cInord = cI + Nx; 
   label cIsouth = cI - Nx;
   scalar p_0 = p[cI];
   scalar rho_0 = rho[cI];
   scalar p_n = p[cInord];      scalar p0_n   = 0.0;
   scalar rho_n = rho[cInord];  scalar rho0_n = 0.0;
   //scalar p_e = p[cIeast];      scalar p0_e   = 0.0;
   //scalar rho_e = rho[cIeast];  scalar rho0_e = 0.0;
   //scalar p_w = p[cIwest];      scalar p0_w   = 0.0;
   //scalar rho_w = rho[cIwest];  scalar rho0_w = 0.0;
   scalar p_s = p[cIsouth];     scalar p0_s   = 0.0;
   scalar rho_s = rho[cIsouth]; scalar rho0_s = 0.0;
   scalar yC = cellCentre[cI].y();

   scalar pPertcI = 0.0;
   scalar rhoPertcI = 0.0;

   GEA::ExtrapolateEquilibrum
   (
    p_0,
    rho_0,
    yC,
    cellCentre[cInord].y(),
    p0_n,
    rho0_n
   );
   scalar pPertNord = p_n - p0_n;
   scalar rhoPertNord = rho_n - rho0_n;
  // ----------------- //

   //------- South------ //
   GEA::ExtrapolateEquilibrum
   (
     p_0,
     rho_0,
     yC,
     cellCentre[cIsouth].y(),
     p0_s,
     rho0_s
   );

   scalar pPertSouth = p_s - p0_s;
   scalar rhoPertSouth = rho_s - rho0_s;

   scalar Dpy = 0.0;
    // y-LIMITER-pressure
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cInord],
     cellCentre[cIsouth],
     pPertcI,
     pPertNord,
     pPertSouth,
     Dpy
    );

    scalar Drhoy = 0.0;
    // y-LIMITER-density
    GEA::CalculateLocalLimiter
    (
     cellCentre[cI],
     cellCentre[cInord],
     cellCentre[cIsouth],
     rhoPertcI,
     rhoPertNord,
     rhoPertSouth,
     Drhoy
    );

   Dp[cI].component(1)   = Dpy;
   Drho[cI].component(1) = Drhoy;

  }// end loop first boundary

 }//end ConstructLocalLimiter

void ConstructVLimiter
(
 volTensorField& DV,
 const volVectorField& V
)
{
   const fvMesh& mesh = V.mesh();
   const volVectorField& cellCentre = mesh.C();
   scalar Lx = 25600.0;
   scalar Ly = 6400.0;
   scalar LxLy = (Lx*Ly);
   label Ncells = cellCentre.size();

   //loop 
   scalar delta = Foam::sqrt(LxLy/(Ncells));
   scalar Nx_ = Lx/delta;
   scalar Ny_ = Ly/delta;
   int Nx = int(Nx_);
   int Ny = int(Ny_);

   for(int Ix=1; Ix<Nx-1;Ix++)
   {
    for(int Iy=1; Iy<Ny-1;Iy++)
    {
      label cI = Ix + Iy*Nx ; // Enumeration of OF 
      label cInord = cI + Nx;
      label cIeast = cI + 1;
      label cIwest = cI - 1;
      label cIsouth = cI - Nx;
      scalar yC = cellCentre[cI].y();

      vector u_i= V[cI];
      vector u_e= V[cIeast];
      vector u_w= V[cIwest];
      vector u_n= V[cInord];
      vector u_s= V[cIsouth];

      // Horizontal
      vector Du_x(0.0,0.0,0.0);
      CalculateLocalLimiter
      (
       cellCentre[cI],
       cellCentre[cIeast],
       cellCentre[cIwest],
       u_i,
       u_e,
       u_w,
       Du_x
      );
      DV[cI].xx() = Du_x.component(0); // dudx
      DV[cI].yx() = Du_x.component(1); //dvdx

      // Vertical 
      vector Du_y(0.0,0.0,0.0);
      CalculateLocalLimiter
      (
       cellCentre[cI],
       cellCentre[cInord],
       cellCentre[cIsouth],
       u_i,
       u_n,
       u_s,
       Du_y
      );
      DV[cI].xy() = Du_y.component(0);
      DV[cI].yy() = Du_y.component(1);
     }
   }
}// end Limiter for the velocity 

}// end namespace
