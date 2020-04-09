#include "VelocityStabilized.H"

using namespace Foam;

preciceAdapter::FSI::VelocityStabilized::VelocityStabilized
(
    const Foam::fvMesh& mesh,
    const std::string namePointVelocity,
	const Foam::scalar velocStabFac
)
:
mesh_(mesh),
pointVelocity_(
    const_cast<pointVectorField*>
    (
        &mesh.lookupObject<pointVectorField>(namePointVelocity)
    )
),
velocStabFac_(velocStabFac)
{
    dataType_ = vector;
}

void preciceAdapter::FSI::VelocityStabilized::write(double * buffer, bool meshConnectivity)
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

// return the VelocityStabilized to use later in the VelocityStabilized?
void preciceAdapter::FSI::VelocityStabilized::read(double * buffer)
{
    // For every element in the buffer
    int bufferIndex = 0;

    // necessary fields
    const volScalarField& p =
            mesh_.lookupObject<volScalarField>("p");

    //parameters for penalization
    volScalarField dp = (p.oldTime()-p.oldTime().oldTime())/mesh_.time().deltaT();
    Foam::scalar gamma0 = 1;
    Foam::scalar mu = 1e-3;
    // gamma0 = 1
    // gamma = smallest_cell_size / time_step
    // pressure_time_derivative = (p_n-p_nm1)/time_step

    //Info << dp << endl;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // Get the VelocityStabilized on the patch
        fixedValuePointPatchVectorField& pointVelocityFluidPatch =
            refCast<fixedValuePointPatchVectorField>
            (
            		pointVelocity_->boundaryFieldRef()[patchID]
            );

		//- get mesh normals at patch
		const vectorField& Nf =
		        mesh_.boundary()[patchID].nf();

		// calculate penalization term (simplified!!)
		vectorField penal = (gamma0 / mu)* velocStabFac_ * Nf * dp.boundaryFieldRef()[patchID];

		//- set-up interpolator
		primitivePatchInterpolation patchInterpolator
		(
			dp.mesh().boundaryMesh()[patchID]
		);

		// perform interpolation
		vectorField penalPoint = patchInterpolator.faceToPointInterpolate(penal);


		// no penalty in first steps
		if (mesh_.time().value() < 20 * mesh_.time().deltaT().value())
		{
			forAll (penalPoint, i)
			{
				penalPoint[i][0] = 0.;
				penalPoint[i][1] = 0.;
				penalPoint[i][2] = 0.;
			}
		}

//		Info << "penalPoint: " << endl;
//		Info << penalPoint << endl;

		// For every cell of the patch
		forAll(pointVelocity_->boundaryFieldRef()[patchID], i)
		{
//			Info << "INDICES OUTPUT" << endl;
//			Info << i << endl;
			// Calculate the stabilized velocity
			pointVelocityFluidPatch[i][0] = buffer[bufferIndex++] + penalPoint[i][0];
			pointVelocityFluidPatch[i][1] = buffer[bufferIndex++] + penalPoint[i][1];
			pointVelocityFluidPatch[i][2] = buffer[bufferIndex++] + penalPoint[i][2];
		}


    }
}
