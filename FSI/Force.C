#include "Force.H"

using namespace Foam;

preciceAdapter::FSI::Force::Force
(
    const Foam::fvMesh& mesh,
    const fileName& timeName,
	const std::string nameRho,
	const std::string nameU
)
:
mesh_(mesh),
nameRho_(nameRho),
nameU_(nameU)
{
    dataType_ = vector;

    Force_ = new volVectorField
    (
        IOobject
        (
            "Force",
            timeName,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "fdim",
            dimensionSet(1,1,-2,0,0,0,0),
            Foam::vector::zero
        )
    );
}

Foam::tmp<Foam::volSymmTensorField> preciceAdapter::FSI::Force::devRhoReff() const
{
	typedef compressible::turbulenceModel cmpTurbModel;
	typedef incompressible::turbulenceModel icoTurbModel;

	//- Turbulent, compressible
	if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
	{
		const cmpTurbModel& turb =
				mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

		return turb.devRhoReff();
	}

	//- Turbulent, incompressible
	else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
	{
		const incompressible::turbulenceModel& turb =
				mesh_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

		return rho()*turb.devReff();
	}

	//- Laminar, compressible
	else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
	{
		const fluidThermo& thermo =
				mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

		const volVectorField& U = mesh_.lookupObject<volVectorField>(nameU_);

		return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
	}

	//- Laminar transport model
	else if(mesh_.foundObject<transportModel>("transportProperties"))
	{
		const transportModel& laminarT =
				mesh_.lookupObject<transportModel>("transportProperties");

		const volVectorField& U = mesh_.lookupObject<volVectorField>(nameU_);

		return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
	}

	//- Direct calculation of laminar viscous stresses
	else if (mesh_.foundObject<dictionary>("transportProperties"))
	{
		const dictionary& transportProperties =
				mesh_.lookupObject<dictionary>("transportProperties");

		dimensionedScalar nu(transportProperties.lookup("nu"));

		const volVectorField& U = mesh_.lookupObject<volVectorField>(nameU_);

		return -rho() * nu * dev(twoSymm(fvc::grad(U)));
	}

	else
	{
	    FatalErrorInFunction
	        << "No valid model for viscous stress calculation available."
	        << exit(FatalError);
	}

}

tmp<volScalarField> preciceAdapter::FSI::Force::rho() const
{
	if (mesh_.foundObject<volScalarField>(nameRho_))
	{
		return mesh_.lookupObject<volScalarField>(nameRho_);
	}

	else
	{
		const dictionary& transportProperties =
				mesh_.lookupObject<IOdictionary>("transportProperties");

		//doubleScalar rhoConst = transportProperties.lookup(nameRho_);

		//Info << "Scalar read: " << rhoConst << endl;

		return tmp<volScalarField>
		(
				new volScalarField
				(
						IOobject
						(
								nameRho_,
								mesh_.time().timeName(),
								mesh_
						),
						mesh_,
						dimensionedScalar(nameRho_, dimDensity, transportProperties)
				)
		);
	}
}

doubleScalar preciceAdapter::FSI::Force::rho(const volScalarField& p) const
{
	if (p.dimensions() == dimPressure)
	{
		return 1.0;
	}
	else
	{
		if (!mesh_.foundObject<volScalarField>(nameRho_))
		{

			Info << "rho(p) transportProperties.lookup:"<< endl;
			const dictionary& transportProperties =
					mesh_.lookupObject<IOdictionary>("transportProperties");

			return dimensionedScalar(nameRho_, dimDensity, transportProperties).value();
		}
		else
		{
			FatalErrorInFunction
				<< "Rho cannot be a field if kinematic pressure is provided."
				<< exit(FatalError);
		}
	}
}

void preciceAdapter::FSI::Force::write(double * buffer, bool meshConnectivity)
{
    // Compute forces. See the Forces function object.
    // Normal vectors on the boundary, multiplied with the face areas
    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    // Stress tensor
    tmp<volSymmTensorField> tdevRhoReff = devRhoReff();

    // Stress tensor boundary field
    const volSymmTensorField::Boundary& devRhoReffb
        = tdevRhoReff().boundaryField();
   
    // Pressure
    const volScalarField& p =
        mesh_.lookupObject<volScalarField>("p");

    // TODO: Add reference pressure

    int bufferIndex = 0;
    
    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Pressure forces
        Force_->boundaryFieldRef()[patchID] =
            rho(p) * Sfb[patchID] * p.boundaryField()[patchID];

        // Viscous forces
        Force_->boundaryFieldRef()[patchID] +=
            Sfb[patchID] & devRhoReffb[patchID];

        // Write the forces to the preCICE buffer
        // For every cell of the patch
        forAll(Force_->boundaryFieldRef()[patchID], i)
        {
            // Copy the force into the buffer
            // x-dimension
            buffer[bufferIndex++]
            = 
            Force_->boundaryFieldRef()[patchID][i].x();

            // y-dimension
            buffer[bufferIndex++]
            =
            Force_->boundaryFieldRef()[patchID][i].y();

            // z-dimension
            buffer[bufferIndex++]
            =
            Force_->boundaryFieldRef()[patchID][i].z();
        }
    }
}

void preciceAdapter::FSI::Force::read(double * buffer)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE readBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Reading forces is not supported."
        << exit(FatalError);
}

preciceAdapter::FSI::Force::~Force()
{
    delete Force_;
}
