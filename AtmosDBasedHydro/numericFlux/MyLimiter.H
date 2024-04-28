using namespace Foam;
template<class Type>
void MonotonizedCentralLimiter
(
 Type& dL_,
 const Type& a,
 const Type& b
)
{
 Type& dL = dL_; // Reference
 for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
 {
   typename Type::cmptType aCmpt = a.component(cmpt);
   typename Type::cmptType bCmpt = b.component(cmpt);
   typename Type::cmptType signApsignB = Foam::sign(aCmpt)+Foam::sign(bCmpt);
   typename Type::cmptType minAbsaAbsb = Foam::min(Foam::mag(aCmpt),Foam::mag(bCmpt));
   dL.component(cmpt) = 0.5*(signApsignB)*Foam::min(2*minAbsaAbsb,0.5*Foam::mag(aCmpt+bCmpt));
 }//end for
}

// Template specialization
template <>
void MonotonizedCentralLimiter<scalar>
(
 scalar& dL_,
 const scalar& a,
 const scalar& b
)
{
  scalar& dL = dL_;

  dL = 0.5*(Foam::sign(a)+Foam::sign(b))*
     Foam::min(2*Foam::min(Foam::mag(a),Foam::mag(b)),0.5*Foam::mag(a+b));
}//end MonotizedCentralLimiter


template<class Type>
void MinModLimiter
(
 Type& dL_,
 const Type& a,
 const Type& b
)
{
 Type& dL = dL_;
 for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
 {
   typename Type::cmptType aCmpt = a.component(cmpt);
   typename Type::cmptType bCmpt = b.component(cmpt);
   dL.component(cmpt) = 0.5*(Foam::sign(aCmpt)+Foam::sign(bCmpt))*Foam::min(Foam::mag(aCmpt),Foam::mag(bCmpt));
 }//end for
}//end MinMod

template <>
void MinModLimiter<scalar>
(
 scalar& dL_,
 const scalar& a,
 const scalar& b
)
{
  scalar& dL = dL_;
  dL = 0.5*(Foam::sign(a)+Foam::sign(b))*Foam::min(Foam::mag(a),Foam::mag(b));
}
