Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<basicPsiThermo> thermo
    (
        basicPsiThermo::New(mesh)
    );
  
    	
    dimensionedScalar Cp(GEA::PerfectGasConstant::CpDim);
    dimensionedScalar Cv(GEA::PerfectGasConstant::CvDim);
    dimensionedScalar R(GEA::PerfectGasConstant::RDim);
    scalar gamma = GEA::PerfectGasConstant::gamma; 

    volScalarField& h = thermo->h(); // Used to update the cell-center temperature.

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
	IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
   
    Info<< "Reading field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
 
   volScalarField thetaP
   (
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T*scalar(0.0),
    Foam::zeroGradientFvPatchField<Foam::scalar>::typeName
   );
   

    volScalarField rho0("rho0",rho);
    volScalarField T0("T0",T);
    volScalarField p0("p0",(rho0*T0*R));
    
 
   #include "createBackground.H"
// #include "createStratifiedAtmosphere.H"  
    p0 = rho0*R*T0;

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Conservative variables

    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    volScalarField Phi
    (
     IOobject
     (
       "Phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
        ),
        GEA::PerfectGasConstant::gDim*mesh.C().component(1)
    );
    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(Cv*T + 0.5*magSqr(U)) 
    );


    // Create numeric flux
    autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New
    (
        rho,
        U,
        T,
        thermo()
    );
    basicNumericFlux& dbnsFlux = dbnsFluxPtr();

    // Create mass flux alias for easier coupling with other code components
    const surfaceScalarField& phi = dbnsFlux.rhoFlux();
    
    //- Reference to the pressure : 
    volScalarField& p = dbnsFlux.p(); // 
    
    

    
    
