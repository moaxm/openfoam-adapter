#ifndef FSI_VELOCITY_H
#define FSI_VELOCITY_H

#include "CouplingDataUser.H"

#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads Velocity
class Velocity : public CouplingDataUser
{

private:

    // Velocity pointVectorField
    Foam::pointVectorField * pointVelocity_;

public:

    //- Constructor
    Velocity
    (
        const Foam::fvMesh& mesh,
        const std::string namePointVelocity
    );

    //- Write the Velocity values into the buffer
    void write(double * buffer, bool meshConnectivity, const unsigned int dim);

    //- Read the Velocity values from the buffer
    void read(double * buffer, const unsigned int dim);
};

}
}

#endif
