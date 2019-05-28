#include "FF.H"

#include "Utilities.H"

using namespace Foam;

preciceAdapter::FF::FluidFluid::FluidFluid
(
    const Foam::fvMesh& mesh
)
:
mesh_(mesh)
{}


std::string preciceAdapter::FF::FluidFluid::determineSolverType()
{
    // NOTE: When coupling a different variable, you may want to
    // add more cases here. Or you may provide the solverType in the config.

    std::string solverType;
    
    // TODO: Implement properly
    solverType = "incompressible";
    
    return solverType;
}

void preciceAdapter::FF::FluidFluid::addWriters(std::string dataName, Interface * interface)
{
    if (dataName.find("VelocityGradient") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new VelocityGradient(mesh_, nameU_)
        );
        DEBUG(adapterInfo("Added writer: Velocity Gradient."));
    }
    else if (dataName.find("Velocity") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new Velocity(mesh_, nameU_)
        );
        DEBUG(adapterInfo("Added writer: Velocity."));
    }
    else if (dataName.find("PressureGradient") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new PressureGradient(mesh_, nameP_)
        );
        DEBUG(adapterInfo("Added writer: Pressure Gradient."));
    }
    else if (dataName.find("Pressure") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new Pressure(mesh_, nameP_)
        );
        DEBUG(adapterInfo("Added writer: Pressure."));
    }
    else
    {
        adapterInfo("Unknown data type - cannot add " + dataName +".", "error");
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a reader below).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.
}

void preciceAdapter::FF::FluidFluid::addReaders(std::string dataName, Interface * interface)
{
    if (dataName.find("VelocityGradient") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new VelocityGradient(mesh_, nameU_)
        );
        DEBUG(adapterInfo("Added reader: VelocityGradient."));
    }
    else if (dataName.find("Velocity") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new Velocity(mesh_, nameU_)
        );
        DEBUG(adapterInfo("Added reader: Velocity."));
    }
    else if (dataName.find("PressureGradient") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new PressureGradient(mesh_, nameP_)
        );
        DEBUG(adapterInfo("Added reader: Pressure Gradient."));
    }
    else if (dataName.find("Pressure") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new Pressure(mesh_, nameP_)
        );
        DEBUG(adapterInfo("Added reader: Pressure."));
    }
    else
    {
        adapterInfo("Unknown data type - cannot add " + dataName +".", "error");
    }

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // reader here (and as a writer above).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.
}


