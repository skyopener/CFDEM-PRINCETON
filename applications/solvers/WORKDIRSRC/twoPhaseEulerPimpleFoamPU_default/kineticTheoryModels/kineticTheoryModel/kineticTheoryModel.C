/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kineticTheoryModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const Foam::phaseModel& phasea,
    const Foam::volVectorField& Ub,
    const Foam::volScalarField& alpha,
    const Foam::dragModel& draga
)
:
    phasea_(phasea),
    Ua_(phasea.U()),
    Ub_(Ub),
    alpha_(alpha),
    phia_(phasea.phi()),
    draga_(draga),

    rhoa_(phasea.rho()),
    da_(phasea.d()),
    nua_(phasea.nu()),

    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            Ua_.time().constant(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    kineticTheory_(kineticTheoryProperties_.lookup("kineticTheory")),
    equilibrium_(kineticTheoryProperties_.lookup("equilibrium")),

    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    conductivityModel_
    (
        conductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    radialModel_
    (
        radialModel::New
        (
            kineticTheoryProperties_
        )
    ),
    granularPressureModel_
    (
        granularPressureModel::New
        (
            kineticTheoryProperties_
        )
    ),
    frictionalStressModel_
    (
        frictionalStressModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_(kineticTheoryProperties_.lookup("e")),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    alphaMinFriction_(kineticTheoryProperties_.lookup("alphaMinFriction")),
    Fr_(kineticTheoryProperties_.lookup("Fr")),
    eta_(kineticTheoryProperties_.lookup("eta")),
    p_(kineticTheoryProperties_.lookup("p")),
    phi_(dimensionedScalar(kineticTheoryProperties_.lookup("phi"))*M_PI/180.0),
    Theta_
    (
        IOobject
        (
            "Theta",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh()
    ),
    mua_
    (
        IOobject
        (
            "mua",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    pa_
    (
        IOobject
        (
            "pa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    ppMagf_
    (
        IOobject
        (
            "ppMagf",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    gs0_
    (
        IOobject
        (
            "gs0",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
    ),
    gs0Prime_
    (
        IOobject
        (
            "gs0prime",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),
//_AO_09/01/2014 
// Shear stress ratio eta (here, we call it as upsilon)    
    upsilon_
    (
        IOobject
        (
            "upsilon",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ), 
    
//_AO_09/01/2014 
// Shear stress tau_ (here, we call it as upsilon)    
    tau_
    (
        IOobject
        (
            "tau",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedSymmTensor("zero", dimensionSet(1, -1, -2, 0, 0), symmTensor(0,0,0,0,0,0))
    ),
//_AO_09/01/2014 
// Effective restitution coefficient (function of e and mu (friction coefficient) )
    eEff_(0),
// Friction coefficient 
    muFric_(0),    	    
// Dilute inertial regime 0 < alpha < alphaf (p. 2) alphaf ~ 0.49
    alphaf_(0),
// Dense intertial regime
    alphac_(0),     
// Yield stress ratio
    upsilons_(0), 
// Modified kinetic theory by Chialvo-Sundaresan on/off     
    mofidiedKineticTheoryPU_(kineticTheoryProperties_.lookup("modifiedKineticTheoryPU"))        
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> Foam::kineticTheoryModel::ppMagfCalc(const volScalarField& alphaUpdate)
{
  const scalar constSMALL = 1.e-06;

  if(mofidiedKineticTheoryPU_)
  {
  // Read from dictionary
   alphaf_ = readScalar(kineticTheoryProperties_.lookup("alphaDiluteInertialUpperLimit"));
   alphac_ = readScalar(kineticTheoryProperties_.lookup("alphaCritical"));
  }else
  {
   alphaf_ = alphaMax_.value();
   alphac_ = alphaMax_.value();
  }
  volScalarField alpha = alphaUpdate;

  //gs0_ = radialModel_->g0(min(alpha, alphaMinFriction_), alphaMax_); 
  //gs0Prime_ = radialModel_->g0prime(min(alpha, alphaMinFriction_), alphaMax_);

  // 
  gs0_ = radialModel_->g0jamming
  (
      Ua_.mesh(),
      min(max(alpha, scalar(constSMALL)), alphac_ - 0.01),
      alphaMax_,
      alphaf_,
      alphac_
  );

  gs0Prime_ = radialModel_->g0jammingPrime
  (
      Ua_.mesh(),
      min(max(alpha, scalar(constSMALL)), alphac_ - 0.01),
      alphaMax_,
      alphaf_,
      alphac_
  );	

  volScalarField ppMagfFriction = frictionalStressModel_->frictionalPressurePrime
  (
      alpha, 
      alphaMinFriction_, 
      alphaMax_,
      Fr_,
      eta_,
      p_
  );

  // Computing ppMagf
  ppMagf_ = Theta_*granularPressureModel_->granularPressureCoeffPrime
  (
      alpha, 
      gs0_, 
      gs0Prime_, 
      rhoa_, 
      e_
  );  

  // NOTE: this might not be appropriate if J&J model is used (verify)
  forAll(alpha, cellI)
  {
      if(alpha[cellI] >= alphaMinFriction_.value())
      {
	  ppMagf_[cellI] = ppMagfFriction[cellI];
      }
  }

  ppMagf_.correctBoundaryConditions();

  Info<< "kinTheory: min(ppMagf) = " << min(ppMagf_).value()
      << ", max(ppMagf) = "          << max(ppMagf_).value() << endl;

  return ppMagf_;
}

void Foam::kineticTheoryModel::solve(const volTensorField& gradUat)
{
 if(kineticTheory_)
 {
     //if (!kineticTheory_)
     //{
     //    return;
     //}

     //const scalar sqrtPi = sqrt(mathematicalConstant::pi);

     const scalar sqrtPi = sqrt(constant::mathematical::pi);

     surfaceScalarField phi = 1.5*rhoa_*phia_*fvc::interpolate(alpha_);

     volTensorField dU = gradUat.T();//fvc::grad(Ua_);
     volSymmTensorField D = symm(dU);

     // NB, drag = K*alpha*beta,
     // (the alpha and beta has been extracted from the drag function for
     // numerical reasons)
     volScalarField Ur = mag(Ua_ - Ub_);
     volScalarField betaPrim = alpha_*(1.0 - alpha_)*draga_.K(Ur);

     // Calculating the radial distribution function (solid volume fraction is
     //  limited close to the packing limit, but this needs improvements)
     //  The solution is higly unstable close to the packing limit.
     gs0_ = radialModel_->g0
     (
         min(max(alpha_, 1e-6), alphaMax_ - 0.01),
         alphaMax_
     );

     // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
     volScalarField PsCoeff = granularPressureModel_->granularPressureCoeff
     (
         alpha_,
         gs0_,
         rhoa_,
         e_
     );

     // 'thermal' conductivity (Table 3.3, p. 49)
     kappa_ = conductivityModel_->kappa(alpha_, Theta_, gs0_, rhoa_, da_, e_);

     // particle viscosity (Table 3.2, p.47)
     mua_ = viscosityModel_->mua(alpha_, Theta_, gs0_, rhoa_, da_, e_);

     dimensionedScalar Tsmall
     (
         "small",
         dimensionSet(0 , 2 ,-2 ,0 , 0, 0, 0),
         1.0e-6
     );

     dimensionedScalar TsmallSqrt = sqrt(Tsmall);
     volScalarField ThetaSqrt = sqrt(Theta_);

     // dissipation (Eq. 3.24, p.50)
     volScalarField gammaCoeff =
         12.0*(1.0 - sqr(e_))*sqr(alpha_)*rhoa_*gs0_*(1.0/da_)*ThetaSqrt/sqrtPi;

     // Eq. 3.25, p. 50 Js = J1 - J2
     volScalarField J1 = 3.0*betaPrim;
     volScalarField J2 =
         0.25*sqr(betaPrim)*da_*sqr(Ur)
	/(max(alpha_, 1e-6)*rhoa_*sqrtPi*(ThetaSqrt + TsmallSqrt));

     // bulk viscosity  p. 45 (Lun et al. 1984).
     lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0+e_)*ThetaSqrt/sqrtPi;

     // stress tensor, Definitions, Table 3.1, p. 43
     volSymmTensorField tau = 2.0*mua_*D + (lambda_ - (2.0/3.0)*mua_)*tr(D)*I;

     if (!equilibrium_)
     {
         // construct the granular temperature equation (Eq. 3.20, p. 44)
         // NB. note that there are two typos in Eq. 3.20
         // no grad infront of Ps
         // wrong sign infront of laplacian
         fvScalarMatrix ThetaEqn
         (
             fvm::ddt(1.5*alpha_*rhoa_, Theta_)
           + fvm::div(phi, Theta_, "div(phi,Theta)")
          ==
             fvm::SuSp(-((PsCoeff*I) && dU), Theta_)
           + (tau && dU)
           + fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
           + fvm::Sp(-gammaCoeff, Theta_)
           + fvm::Sp(-J1, Theta_)
           + fvm::Sp(J2/(Theta_ + Tsmall), Theta_)
         );

         ThetaEqn.relax();
         ThetaEqn.solve();
     }
     else
     {
         // equilibrium => dissipation == production
         // Eq. 4.14, p.82
         volScalarField K1 = 2.0*(1.0 + e_)*rhoa_*gs0_;
         volScalarField K3 = 0.5*da_*rhoa_*
             (
                 (sqrtPi/(3.0*(3.0-e_)))
        	*(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha_*gs0_)
        	+1.6*alpha_*gs0_*(1.0 + e_)/sqrtPi
             );

         volScalarField K2 =
             4.0*da_*rhoa_*(1.0 + e_)*alpha_*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0;

         volScalarField K4 = 12.0*(1.0 - sqr(e_))*rhoa_*gs0_/(da_*sqrtPi);

         volScalarField trD = tr(D);
         volScalarField tr2D = sqr(trD);
         volScalarField trD2 = tr(D & D);

         volScalarField t1 = K1*alpha_ + rhoa_;
         volScalarField l1 = -t1*trD;
         volScalarField l2 = sqr(t1)*tr2D;
         volScalarField l3 = 4.0*K4*max(alpha_, 1e-6)*(2.0*K3*trD2 + K2*tr2D);

         Theta_ = sqr((l1 + sqrt(l2 + l3))/(2.0*(alpha_ + 1.0e-4)*K4));
     }

     Theta_.max(1.0e-15);
     Theta_.min(1.0e+3);

     volScalarField pf = frictionalStressModel_->frictionalPressure
     (
         alpha_,
         alphaMinFriction_,
         alphaMax_,
         Fr_,
         eta_,
         p_
     );

     PsCoeff += pf/(Theta_+Tsmall);

     PsCoeff.min(1.0e+10);
     PsCoeff.max(-1.0e+10);

     // update particle pressure
     pa_ = PsCoeff*Theta_;

     // frictional shear stress, Eq. 3.30, p. 52
     volScalarField muf = frictionalStressModel_->muf
     (
         alpha_,
         alphaMax_,
         pf,
         D,
         phi_
     );

    // add frictional stress
     mua_ += muf;
     mua_.min(1.0e+2);
     mua_.max(0.0);

     Info<< "kinTheory: max(Theta) = " << max(Theta_).value() << endl;

     volScalarField ktn = mua_/rhoa_;

     Info<< "kinTheory: min(nua) = " << min(ktn).value()
         << ", max(nua) = " << max(ktn).value() << endl;

     Info<< "kinTheory: min(pa) = " << min(pa_).value()
         << ", max(pa) = " << max(pa_).value() << endl;
 //}

 /*
 volScalarField& Foam::kineticTheoryModel::ppMagf(const volScalarField& alphaUpdate)
 {
     volScalarField alpha = alphaUpdate;

     gs0_ = radialModel_->g0(min(alpha, alphaMinFriction_), alphaMax_); 
     gs0Prime_ = radialModel_->g0prime(min(alpha, alphaMinFriction_), alphaMax_);

     // Computing ppMagf
     ppMagf_ = Theta_*granularPressureModel_->granularPressureCoeffPrime
     (
	 alpha, 
	 gs0_, 
	 gs0Prime_, 
	 rhoa_, 
	 e_
     );

     volScalarField ppMagfFriction = frictionalStressModel_->frictionalPressurePrime
     (
	 alpha, 
	 alphaMinFriction_, 
	 alphaMax_,
         Fr_,
         eta_,
         p_
     );

     // NOTE: this might not be appropriate if J&J model is used (verify)
     forAll(alpha, cellI)
     {
	 if(alpha[cellI] >= alphaMinFriction_.value())
	 {
	     ppMagf_[cellI] = ppMagfFriction[cellI];
	 }
     }

     ppMagf_.correctBoundaryConditions();

     return ppMagf_;
 }
 */
 
 }
 else if(mofidiedKineticTheoryPU_)
 {
     //if (!mofidiedKineticTheoryPU_)
     //{
     //    return;
     //}
     Info << " " << endl;
     Info << "Modified kinetic theory model - Chialvo-Sundaresan " << endl;
     const scalar Pi = constant::mathematical::pi;
     const scalar sqrtPi = sqrt(constant::mathematical::pi);
     const scalar constSMALL = 1.e-06;

     // Mass flux of solid@surface	
     surfaceScalarField phi(1.5*rhoa_*phia_*fvc::interpolate(alpha_));

     // Read from dictionary
     muFric_ = readScalar(kineticTheoryProperties_.lookup("muFriction"));
     eEff_ = e_ - 3.0 / 2.0 * muFric_ * exp(-3.0 * muFric_);
     alphaf_ = readScalar(kineticTheoryProperties_.lookup("alphaDiluteInertialUpperLimit"));
     alphac_ = readScalar(kineticTheoryProperties_.lookup("alphaCritical"));
     upsilons_ = readScalar(kineticTheoryProperties_.lookup("yieldStressRatio"));

     // NB, drag = K*alpha1*alpha2,
     // (the alpha1 and alpha2 has been extracted from the drag function for
     // numerical reasons)
     volScalarField Ur(mag(Ua_ - Ub_));
     volScalarField alpha2Prim(alpha_*(1.0 - alpha_)*draga_.K(Ur));

     // Calculating the radial distribution function (solid volume fraction is
     //  limited close to the packing limit, but this needs improvements)
     //  The solution is higly unstable close to the packing limit.

     gs0_ = radialModel_->g0jamming
     (
         Ua_.mesh(),
	 //min(max(alpha_, scalar(constSMALL)), alphaMax_ - 0.01),
         min(max(alpha_, scalar(constSMALL)), alphac_ - 0.01),
	 alphaMax_,
	 alphaf_,
	 alphac_
     );

     // particle pressure - coefficient in front of T (Eq. 1, p. 3)
     volScalarField PsCoeff	// -> rho_p * H 
     (
         granularPressureModel_->granularPressureCoeff
         (
             alpha_,
             gs0_,
             rhoa_,
             e_
         )
     );    

     // Model parameters
     dimensionedScalar I0(0.2); // Table 2, p.15
     dimensionedScalar const_alpha(0.36); // Table 2, p.15
     dimensionedScalar const_alpha1(0.06); // Table 2, p.15

     // Solid kinetic+collisional viscosity mua_ = nu_k^star + nu_c^star, Eq. 8,9, p.4
     // If Garzo-Dufty viscosity is used (viscosity is dimensionless), there is issue with dimension of mu1
     // Create dimensionedScalar
     dimensionedScalar mu1Dim("zero", dimensionSet(1, -1, -1, 0, 0), 1.0);     
     mua_ = mu1Dim * viscosityModel_->mua(alpha_, Theta_, gs0_, rhoa_, da_, e_);

     // Solid bulk viscosity mua_ = nu_k^star + nu_c^star, Eq. 10, p.4
     // If Garzo-Dufty viscosity is used (viscosity is dimensionless), there is issue with dimension of mu1
     // Create dimensionedScalar
     dimensionedScalar lambdaDim("zero", dimensionSet(1, -1, -1, 0, 0), 1.0);     
     lambda_ = lambdaDim * 384.0 / ( 25.0 * Pi ) * ( 1.0 + e_ ) * alpha_ * alpha_ * gs0_ ;  
     //lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0+e_)*sqrt(Theta_)/sqrtPi;

     // J Eq.5, p3
     volScalarField J_( 5.0 * sqrtPi / 96.0 *( mua_ + lambda_ ) /( mu1Dim ) ); // Dimension issue ( mu1Dim )

     // K Eq.6, p3
     volScalarField K_(12.0/sqrtPi*alpha_*alpha_*gs0_*(1.0-e_*e_));

     // K' Eq.26, p8 modified dissipation due to friction
     volScalarField Kmod_(K_*(1.0 - eEff_*eEff_)/(1.0 - e_*e_));
     dimensionedScalar KmodSmall("KmodSmall",dimensionSet(0 , 0 , 0 , 0 , 0, 0, 0), constSMALL);    

     // M Eq.30 p.9
     volScalarField M_( max( J_ / max( Kmod_, KmodSmall) , const_alpha1 / sqrt( max(alphac_ - alpha_, constSMALL) ) ) ); 

     // Shear stress rate tensor
     volTensorField dU(gradUat.T());    //fvc::grad(Ua_);
     volSymmTensorField D(symm(dU));    
     //volTensorField dU(fvc::grad(Ua_));
     //volSymmTensorField D(symm(fvc::grad(Ua_)));

     // Shear stress rate (gammaDot)
     volScalarField gammaDot(sqrt(magSqr(D)));
     dimensionedScalar gammaDotSmall("gammaDotSmall",dimensionSet(0 , 0 , -1 , 0 , 0, 0, 0), constSMALL);    

     // Dilute inertia temperature Eq.24, p8    
     volScalarField ThetaDil_ = ( J_ / max ( Kmod_ , KmodSmall ) ) * ( gammaDot * da_ ) * ( gammaDot * da_ );

     // Dense inertia temperature Eq.27, p8    
     volScalarField ThetaDense_ =   const_alpha1 * ( gammaDot * da_ ) * ( gammaDot * da_ )
                                  / sqrt( max(alphac_ - alpha_, constSMALL) ); 


     // Following equations for transport of equation of temperature
     // grad(Us)
     //volTensorField dU(gradUat.T());    //fvc::grad(Ua_);
     // 'thermal' conductivity (Table 3.3, p. 49) Van Wachem
     kappa_ = conductivityModel_->kappa(alpha_, Theta_, gs0_, rhoa_, da_, e_);

     dimensionedScalar Tsmall
     (
         "small",
         dimensionSet(0 , 2 ,-2 ,0 , 0, 0, 0),
         constSMALL
     );

     dimensionedScalar TsmallSqrt = sqrt(Tsmall);
     volScalarField ThetaSqrt(sqrt(Theta_));

     // Eq. 3.25, p. 50 Js = J1 - J2
     volScalarField J1(3.0*alpha2Prim);
     volScalarField J2
     (
         0.25*sqr(alpha2Prim)*da_*sqr(Ur)
	/(max(alpha_, scalar(1e-6))*rhoa_*sqrtPi*(ThetaSqrt + TsmallSqrt))
     );

     // Psi Eq.32, p.12
     dimensionedScalar psi(1.0 + 3.0/10.0*pow((1.0-e_*e_),-1.5)*(1.0-exp(-8.0*muFric_)));

     // Shear stress ratio in dilute regime, Eq.33, p.12
     dimensionedScalar paSmall("paSmall",dimensionSet(1, -1, -2, 0, 0), constSMALL);    
     //volScalarField inertiaNumber( gammaDot * da_ / sqrt( (pa_ + paSmall) / rhoa_ ) );
     // Use pressure at the previous time-step
     volScalarField inertiaNumber( gammaDot * da_ / sqrt( (pa_ + paSmall) / rhoa_ ) );

     // Modified inertia number Eq.35, p.13
     volScalarField modInertiaNumber( inertiaNumber /  max( alpha_, constSMALL ) ); 

     // Model parameters
     volScalarField chi( 1.0 / ( pow( I0 / max( modInertiaNumber,constSMALL ) , 1.5 ) + 1.0 ));

     // Beta + Sigma_tau Eq.49 p.14
     volScalarField beta(alpha_ * psi * J_ * sqrt( K_ /( max ( (Kmod_ * ( PsCoeff / rhoa_)), constSMALL ) ) ) ); 

     volScalarField sigmaTau( const_alpha / max( beta, constSMALL )  + ( 1 - const_alpha / max( beta, constSMALL ) ) * chi);

     // Sigma_gamma Eq.51 p.14
     volScalarField sigmaGamma( beta * sqrt(PsCoeff/rhoa_) / max( ( Kmod_ * M_), constSMALL ) * sigmaTau);

     // dissipation
     volScalarField gammaCoeff
     (
         // van Wachem  (Eq. 3.24, p.50) 12.0*(1.0 - sqr(e_))*sqr(alpha_)*rhoa_*gs0_*(1.0/da_)*ThetaSqrt/sqrtPi
         // Chialvo & Sundaresan Eq.50 p.14 
         //rhoa_ / da_ * Kmod_ * Theta_ * sqrt(Theta_) * sigmaGamma
         rhoa_ / da_ * Kmod_ * sqrt(Theta_) * sigmaGamma    
     );

     // stress tensor, Definitions, Table 3.1, p. 43
     volSymmTensorField tau(2.0*mua_*D + (lambda_ - (2.0/3.0)*mua_)*tr(D)*I);

     // Temperature
     if (!equilibrium_)
     {
         // construct the granular temperature equation (Eq. 3.20, p. 44)
         // NB. note that there are two typos in Eq. 3.20
         // no grad infront of Ps
         // wrong sign infront of laplacian
         fvScalarMatrix ThetaEqn
         (
             fvm::ddt(1.5*alpha_*rhoa_, Theta_)
           + fvm::div(phi, Theta_, "div(phi,Theta)")
          ==
             fvm::SuSp(-((PsCoeff*I) && dU), Theta_)
           + (tau && dU)
           + fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
           + fvm::Sp(-gammaCoeff, Theta_)
           + fvm::Sp(-J1, Theta_)
           + fvm::Sp(J2/(Theta_ + Tsmall), Theta_)
         );

         ThetaEqn.relax();
         ThetaEqn.solve();
     }else
     {
	// Relaxation of Theta
	Theta_ = max(ThetaDil_,ThetaDense_) ;
     }

     // Limit granular temperature
     Theta_.max(1.0e-15);
     Theta_.min(1.0e+3);

     // update particle pressure
     // Relaxation of pa_
     pa_ = PsCoeff * Theta_;

     // Blending function    
     volScalarField func_B( const_alpha + ( beta-const_alpha ) * chi );

     // Shear stress ratio
     upsilon_ = upsilons_ * (1 - chi) + func_B * modInertiaNumber;

     // Limit shear stress ratio
     //upsilon_.min(1.e-06);

     // Shear stress
     volSymmTensorField S( D - 1./3.*tr(D)*I);    
     volSymmTensorField hatS( S / max( gammaDot, gammaDotSmall ) );

     tau_ = pa_ * upsilon_ * hatS;

     // Limit solid viscosity 
     mua_.min(1.0e+2);
     mua_.max(0.0);

     Info<< "kinTheory: min(g0) = " << min(gs0_).value()
         << ", max(g0) = "          << max(gs0_).value() << endl;

     Info<< "kinTheory: min(g0Prime) = " << min(gs0Prime_).value()
         << ", max(g0Prime) = "          << max(gs0Prime_).value() << endl;

     Info<< "kinTheory: min(J) = " << min(J_).value()
         << ", max(J) = "          << max(J_).value() << endl;

     Info<< "kinTheory: min(Kmod) = " << min(Kmod_).value()
         << ", max(Kmod) = "          << max(Kmod_).value() << endl;

     Info<< "kinTheory: min(gammadot) = " << min(gammaDot).value()
         << ", max(gammadot) = "          << max(gammaDot).value() << endl;

     Info<< "kinTheory: min(PsCoeff) = " << min(PsCoeff).value()
         << ", max(PsCoeff) = "          << max(PsCoeff).value() << endl;

     Info<< "kinTheory: min(ThetaDil) = "   << min(ThetaDil_).value()
         << ", max(ThetaDil) = "            << max(ThetaDil_).value() << endl;

     Info<< "kinTheory: min(ThetaDense) = " << min(ThetaDense_).value()
         << ", max(ThetaDense) = "          << max(ThetaDense_).value() << endl;	

     Info<< "kinTheory: min(S) = " << min(S).value()
         << ", max(S) = "          << max(S).value() << endl;	

     Info<< "kinTheory: min(hatS) = " << min(hatS).value()
         << ", max(hatS) = "          << max(hatS).value() << endl;

     Info<< "kinTheory: min(modInertiaNumber) = " << min(modInertiaNumber).value()
         << ", max(modInertiaNumber) = "          << max(modInertiaNumber).value() << endl;					

     Info<< "kinTheory: min(beta) = " << min(beta).value()
         << ", max(beta) = "          << max(beta).value() << endl;

     Info<< "kinTheory: min(func_B) = " << min(func_B).value()
         << ", max(func_B) = "          << max(func_B).value() << endl;

     Info<< "kinTheory: min(sigmaTau) = " << min(sigmaTau).value()
         << ", max(sigmaTau) = "          << max(sigmaTau).value() << endl;	

     Info<< "kinTheory: min(chi) = " << min(chi).value()
         << ", max(chi) = "          << max(chi).value() << endl;

     volScalarField ktn(mua_/rhoa_);

     Info<< "kinTheory: min(nu1) = " << min(ktn).value()
         << ", max(nu1) = "          << max(ktn).value() << endl;

     Info<< "kinTheory: min(lamdba) = " << min(lambda_).value()
         << ", max(lamdba) = "          << max(lambda_).value() << endl;

     Info<< "kinTheory: min(upsilon) = " << min(upsilon_).value()
         << ", max(upsilon) = "          << max(upsilon_).value() << endl;	

     Info<< "kinTheory: min(Theta) = "   << min(Theta_).value()
         << ", max(Theta) = "           << max(Theta_).value() << endl;

     Info<< "kinTheory: min(pa) = " << min(pa_).value()
         << ", max(pa) = "          << max(pa_).value() << endl; 
      
     Info<< "kinTheory: min(tau) = " << min(tau_).value()
         << ", max(tau) = "          << max(tau_).value() << endl;

     Info << " " << endl;

 }
 else
 {
   return;
 }

}


// ************************************************************************* //
