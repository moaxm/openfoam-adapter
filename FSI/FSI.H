#ifndef FSI_H
#define FSI_H

#include "Interface.H"

#include "FSI/Displacement.H"
#include "FSI/DisplacementDelta.H"
#include "FSI/Velocity.H"
#include "FSI/VelocityStabilized.H"
#include "FSI/Force.H"

#include "fvCFD.H"

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
    
    //- Solver type
    std::string solverType_ = "none";
    
    //- Name of the pointDisplacement field
    std::string namePointDisplacement_ = "pointDisplacement";

    //- Name of the pointDisplacement field
    std::string namePointVelocity_ = "pointMotionU";

    //~ //- stabilization factors
    //~ Foam::scalar velocStabFac_ = 0.;
    //~ Foam::scalar forceStabFac_ = 0.;

protected:
    
    //- Determine the solver type
    std::string determineSolverType();

    //- Read the FSI-related options from the adapter's configuration file
    bool readConfig(const IOdictionary& adapterConfig);

public:

    //- Constructor
    //  We need also the runTime to access the timeName and the deltaT
    FluidStructureInteraction(const Foam::fvMesh& mesh, const Foam::Time& runTime);

    //- Configure
    bool configure(const IOdictionary& adapterConfig);

    //- Add coupling data writers
    void addWriters(std::string dataName, Interface * interface);

    //- Add coupling data readers
    void addReaders(std::string dataName, Interface * interface);
};

}
}

#endif
