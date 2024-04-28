/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "hllcAusmFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcAusmFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf
) const
{
// Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = (RLeft + CvLeft)/CvLeft;
    const scalar kappaRight = (RRight + CvRight)/CvRight;

    // Compute conservative variables assuming perfect gas law

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*(CvLeft*TLeft+0.5*magSqr(ULeft));
    const scalar rhoERight = rhoRight*(CvRight*TRight+0.5*magSqr(URight));

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (ULeft & normalVector);
    const scalar qRight = (URight & normalVector);

    // Speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft =
        Foam::sqrt(max(0.0,kappaLeft * pLeft/rhoLeft));

    const scalar aRight =
        Foam::sqrt(max(0.0,kappaRight * pRight/rhoRight));

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft-aLeft, qRight-aRight);
    const scalar SRight = max(qLeft+aLeft, qRight+aRight);

    const scalar SStar = (rhoRight*qRight*(SRight-qRight)
      - rhoLeft*qLeft*(SLeft - qLeft) + pLeft - pRight )/
        stabilise((rhoRight*(SRight-qRight)-rhoLeft*(SLeft-qLeft)),VSMALL);

    scalar m_dot = 0.0;

    if (pos(SLeft))
    {
        m_dot = rhoLeft*qLeft;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);
        m_dot = rhoLeft*qLeft + SLeft*(rhoLeft*omegaLeft*(SLeft-qLeft)-rhoLeft);
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);
        m_dot = rhoRight*qRight + SRight*(rhoRight*omegaRight*(SRight-qRight)-rhoRight);
    }
    else if (neg(SRight))
    {
        m_dot = rhoRight*qRight;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }


   const scalar Ku    = 0.3;

    const scalar aTilde = 0.5*(aLeft+aRight);
    const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));

    const scalar MaInf = 0.01;
    const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf)));
    const scalar MaZero    = Foam::sqrt(sqrMaZero);


    const scalar fa = MaZero*(2.0-MaZero);
    const scalar alpha = 3.0/16.0*(-4.0+5.0*sqr(fa));

    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);


    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);

    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
        (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) 
         -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
        (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)
         +16.0*alpha*MaRelRight*Ma2PlusRight)));

    const scalar pU = -Ku*P5alphaPlusLeft*P5alphaMinusRight*(rhoLeft+rhoRight)
                                 *(fa*aTilde)*(qRight-qLeft);  
 
    scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight + pU;


    if(m_dot>0)
    {
        rhoFlux = m_dot * magSf;
        rhoUFlux = (m_dot * ULeft + pTilde * normalVector) *magSf;
        rhoEFlux = (m_dot * HLeft ) *magSf;
    }
    else
    {
        rhoFlux = m_dot * magSf;
        rhoUFlux = (m_dot * URight + pTilde * normalVector) *magSf;
        rhoEFlux = (m_dot * HRight ) *magSf;
    }
  
}


// ************************************************************************* //
void Foam::hllcAusmFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const scalar& RLeft,
    const scalar& RRight,
    const scalar& CvLeft,
    const scalar& CvRight,
    const vector& Sf,
    const scalar& magSf,
    const vector& Cf
) const
{
// Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // Ratio of specific heat capacities
    const scalar kappaLeft = (RLeft + CvLeft)/CvLeft;
    const scalar kappaRight = (RRight + CvRight)/CvRight;

    // Compute conservative variables assuming perfect gas law

    // Density
    const scalar rhoLeft = pLeft/(RLeft*TLeft);
    const scalar rhoRight = pRight/(RRight*TRight);

    // DensityVelocity
    const vector rhoULeft = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // DensityTotalEnergy
    const scalar rhoELeft = rhoLeft*(CvLeft*TLeft+0.5*magSqr(ULeft));
    const scalar rhoERight = rhoRight*(CvRight*TRight+0.5*magSqr(URight));

    // Compute left and right total enthalpies:
    const scalar HLeft = (rhoELeft + pLeft)/rhoLeft;
    const scalar HRight = (rhoERight + pRight)/rhoRight;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (ULeft & normalVector);
    const scalar qRight = (URight & normalVector);

    // Speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft =
        Foam::sqrt(max(0.0,kappaLeft * pLeft/rhoLeft));

    const scalar aRight =
        Foam::sqrt(max(0.0,kappaRight * pRight/rhoRight));

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft-aLeft, qRight-aRight);
    const scalar SRight = max(qLeft+aLeft, qRight+aRight);

    const scalar SStar = (rhoRight*qRight*(SRight-qRight)
      - rhoLeft*qLeft*(SLeft - qLeft) + pLeft - pRight )/
        stabilise((rhoRight*(SRight-qRight)-rhoLeft*(SLeft-qLeft)),VSMALL);

    scalar m_dot = 0.0;

    if (pos(SLeft))
    {
        m_dot = rhoLeft*qLeft;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);
        m_dot = rhoLeft*qLeft + SLeft*(rhoLeft*omegaLeft*(SLeft-qLeft)-rhoLeft);
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);
        m_dot = rhoRight*qRight + SRight*(rhoRight*omegaRight*(SRight-qRight)-rhoRight);
    }
    else if (neg(SRight))
    {
        m_dot = rhoRight*qRight;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }


    const scalar Ku    = 0.75;

    const scalar aTilde = 0.5*(aLeft+aRight);
    const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));

    const scalar MaInf = 0.01;
    const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf)));
    const scalar MaZero    = Foam::sqrt(sqrMaZero);


    const scalar fa = MaZero*(2.0-MaZero);
    const scalar alpha = 3.0/16.0*(-4.0+5.0*sqr(fa));

    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);


    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);

    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
        (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) 
         -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
        (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)
         +16.0*alpha*MaRelRight*Ma2PlusRight)));

    const scalar pU = -Ku*P5alphaPlusLeft*P5alphaMinusRight*(rhoLeft+rhoRight)
                                 *(fa*aTilde)*(qRight-qLeft);  
 
    scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight + pU;


    if(m_dot>0)
    {
        rhoFlux = m_dot * magSf;
        rhoUFlux = (m_dot * ULeft + pTilde * normalVector) *magSf;
        rhoEFlux = (m_dot * HLeft ) *magSf;
    }
    else
    {
        rhoFlux = m_dot * magSf;
        rhoUFlux = (m_dot * URight + pTilde * normalVector) *magSf;
        rhoEFlux = (m_dot * HRight ) *magSf;
    }

    scalar g = 9.81;
    scalar PhiF = g * Cf.component(1);
    rhoEFlux+= ( m_dot * magSf * PhiF );

}
