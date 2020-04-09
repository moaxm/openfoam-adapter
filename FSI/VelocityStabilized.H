#ifndef FSI_VELOCITYSTABILIZED_H
#define FSI_VELOCITYSTABILIZED_H

#include "CouplingDataUser.H"

#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads Velocity
class VelocityStabilized : public CouplingDataUser
{

private:
	//- OpenFOAM fvMesh object (we need to access the objects' registry multiple times)
	const Foam::fvMesh& mesh_;

    // Velocity pointVectorField
    Foam::pointVectorField * pointVelocity_;

    // stabilization factor
    Foam::scalar velocStabFac_;

public:

    //- Constructor
    VelocityStabilized
    (
        const Foam::fvMesh& mesh,
        const std::string namePointVelocity,
		const Foam::scalar velocStabFac
    );

    //- Write the Velocity values into the buffer
    void write(double * buffer, bool meshConnectivity);

    //- Read the Velocity values from the buffer
    void read(double * buffer);
};

}
}

#endif
