#ifndef FSI_H
#define FSI_H

#include "Interface.H"

#include "FSI/Displacement.H"
#include "FSI/DisplacementDelta.H"
#include "FSI/Force.H"

#include "fvCFD.H"

// YAML reader - Used to read the adapter's configuration file.
#include "yaml-cpp/yaml.h"

namespace preciceAdapter
{
namespace FSI
{

class FluidStructureInteraction
{

protected:

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    //- OpenFOAM Time object (we need to access the timestep to compute the pointMotionU)
    const Foam::Time& runTime_;

    //- Name of the pointDisplacement field
    std::string namePointDisplacement_ = "pointDisplacement";

    //- Name of density field, applies for incompressible multiphase solvers and compressible solvers
    //  rho value for single phase incompressible solvers is directly determined in Force.C
    std::string nameRho_ = "rho";

    //- Name of velocity field
    std::string nameU_ = "U";

    //- reference pressure
    const Foam::scalar pRef_;

protected:

    //- Read the FSI-related options from the adapter's configuration file
    bool readConfig(const YAML::Node adapterConfig);

public:

    //- Constructor
    //  We need also the runTime to access the timeName and the deltaT
    FluidStructureInteraction(const Foam::fvMesh& mesh, const Foam::Time& runTime, const Foam::scalar pRef);

    //- Configure
    bool configure(const YAML::Node adapterConfig);

    //- Add coupling data writers
    void addWriters(std::string dataName, Interface * interface);

    //- Add coupling data readers
    void addReaders(std::string dataName, Interface * interface);
};

}
}

#endif
