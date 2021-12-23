    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
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

    Info<< "Creating previous physical velocity field U_Prev" << endl;
    volVectorField U_Prev
    (
        IOobject
        (
            "U_Prev",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );

    //mod by alice
    Info<< "Reading physical velocity field U" << endl;
    Info<< "Note: only if voidfraction at boundary is 1, U is superficial velocity!!!\n" << endl;
    volVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    /*Info<< "\nCreating outside field\n" << endl;
    volScalarField outside
    (
        IOobject
        (
            "outside",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        scalar(0)
    );*/

    Info<< "Reading force field Fs" << endl;
    volVectorField Fs
    (
        IOobject
        (
            "Fs",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading force field FInertial" << endl;
    volVectorField FInertial
    (
        IOobject
        (
            "FInertial",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Fs
    );

//========================
// drag law modelling
//========================

    Info<< "\nCreating dummy density field rho = 1\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field phiIB\n" << endl;
    volScalarField phiIB
    (
        IOobject
        (
            "phiIB",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    //mod by alice
    Info<< "Reading field phiIB\n" << endl;
    volScalarField voidfraction
    (
        IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
//========================

#   include "createPhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

//=== dyM ===================

    Info<< "Reading field interFace\n" << endl;
    volScalarField interFace
    (
        IOobject
        (
            "interFace",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        //dimensionedScalar("0", dimensionSet(0, -1, 0, 0, 0), 0.0)
        dimensionedScalar("0", dimensionSet(0, 0, 0, 0, 0), 0.0)
    );

//===========================
