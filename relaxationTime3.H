//////////////////////////////////////////////////////////////////////////////////////

Info<<"Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
	"transportProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ_IF_MODIFIED,
	IOobject::NO_WRITE
    )
);

//////////////////////////////////////////////////////////////////////////////////////

Info<< "Reading droplet diameter\n" << endl;
dimensionedScalar dW
(
    "dW",
    transportProperties.lookup("dW")
);
Info<< "Reading water density\n" << endl;
dimensionedScalar rhoW
(
    "rhoW",
    transportProperties.lookup("rhoW")
);
Info<< "Reading air density\n" << endl;
dimensionedScalar rho
(
    "rho",
    transportProperties.lookup("rho")
);
Info<< "Reading air viscosity\n" << endl;
dimensionedScalar mu
(
    "mu",
    transportProperties.lookup("mu")
);
Info<< "Reading air kinematic viscosity\n" << endl;
dimensionedScalar nu
(
    "nu",
    transportProperties.lookup("nu")
);

//////////////////////////////////////////////////////////////////////////////////////

const dimensionedScalar dummyV("dummyV",dimensionSet(0, 1, -1, 0, 0),1e-10);

Info<< "Computing Reloc\n" << endl;
volScalarField Reloc ("Reloc", ((max(mag(dummyV),mag(U)))*dW*rhoW/mu));		//local reynolds number//[adim]

Info<< "Computing Cw\n" << endl;
volScalarField Cw ("Cw", (24/Reloc+4*pow(Reloc,-0.333)));	//drag coefficient//[adim]

//////////////////////////////////////////////////////////////////////////////////////

Info<< "\nReading gravitationalProperties" << endl;
IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g")); //create gravity vector
const dimensionedScalar magg = mag(g);				//create scalar gravity magnitude
const dimensionedVector gHat = g/magg;				//create vector (0,0,-1)

Info<< "Reading field grav\n" << endl;
volVectorField grav
(
    IOobject
    (
        "grav",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    g
);

//////////////////////////////////////////////////////////////////////////////////////

/*
Settling velocity is the product of the scalar FF (the magnitude of the settling velocity)
and a dummy vector gHat (with the same direction of gravity)
*/

volScalarField FF ("FF",pow(((4*mag(g)*dW)*pow(3*Cw,-1))*((rhoW-rho)/rho),0.5)); 

Info<< "Computing settling velocity\n" << endl;
volVectorField Uw
(
    IOobject
    (
        "Uw",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    FF*gHat
);


Info<< "computing face flux field phiW\n" << endl;
surfaceScalarField phiW
(
    IOobject
    (
        "phiW",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(Uw)
);
//////////////////////////////////////////////////////////////////////////////////////
