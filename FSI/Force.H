#ifndef FSI_FORCE_H
#define FSI_FORCE_H

#include "CouplingDataUser.H"

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads force
class Force : public CouplingDataUser
{

private:

    //- OpenFOAM fvMesh object (we need to access the objects' registry multiple times)
    const Foam::fvMesh& mesh_;
    
    //- Force field
    Foam::volVectorField * Force_;

    //- Rho name
    const std::string nameRho_;

    //- U name
    const std::string nameU_;

    //- Stress tensor (see the OpenFOAM "Forces" function object)
    Foam::tmp<Foam::volSymmTensorField> devRhoReff() const;

    //-  return Density field
    Foam::tmp<Foam::volScalarField> rho() const;

    //- return scalar density
    Foam::scalar rho(const Foam::volScalarField& p) const;

public:

    //- Constructor
    Force
    (
        const Foam::fvMesh& mesh,
        const fileName& timeName, // We create an IOobject and we need the time directory
		const std::string nameRho,
		const std::string nameU
    );

    //- Write the displacement values into the buffer
    void write(double * buffer, bool meshConnectivity);

    //- Read the displacement values from the buffer
    void read(double * buffer);

    //- Destructor
    ~Force();

};

}
}

#endif
