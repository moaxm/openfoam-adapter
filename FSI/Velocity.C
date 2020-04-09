#include "Velocity.H"

using namespace Foam;

preciceAdapter::FSI::Velocity::Velocity
(
    const Foam::fvMesh& mesh,
    const std::string namePointVelocity
)
:
pointVelocity_(
    const_cast<pointVectorField*>
    (
        &mesh.lookupObject<pointVectorField>(namePointVelocity)
    )
)
{
    dataType_ = vector;
}

void preciceAdapter::FSI::Velocity::write(double * buffer, bool meshConnectivity)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE writeBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Writing Velocities is not supported."
        << exit(FatalError);
}

// return the Velocity to use later in the velocity?
void preciceAdapter::FSI::Velocity::read(double * buffer)
{
    // For every element in the buffer
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the Velocity on the patch
        fixedValuePointPatchVectorField& pointVelocityFluidPatch =
            refCast<fixedValuePointPatchVectorField>
            (
                pointVelocity_->boundaryFieldRef()[patchID]
            );


		// For every cell of the patch
		forAll(pointVelocity_->boundaryFieldRef()[patchID], i)
		{
			// Set the Velocity to the received one
			pointVelocityFluidPatch[i][0] = buffer[bufferIndex++];
			pointVelocityFluidPatch[i][1] = buffer[bufferIndex++];
			pointVelocityFluidPatch[i][2] = buffer[bufferIndex++];
		}
    }
}
