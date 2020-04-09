#include "FSI.H"

#include "Utilities.H"

using namespace Foam;

preciceAdapter::FSI::FluidStructureInteraction::FluidStructureInteraction
(
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime,
	const Foam::scalar pRef
)
:
mesh_(mesh),
runTime_(runTime),
pRef_(pRef)
{}

bool preciceAdapter::FSI::FluidStructureInteraction::configure(const YAML::Node adapterConfig)
{
    DEBUG(adapterInfo("Configuring the FSI module..."));

    // Read the FSI-specific options from the adapter's configuration file
    if (!readConfig(adapterConfig)) return false;

    /* TODO: If we need different solver types,
    /  here is the place to determine it.
    */

    return true;
}

bool preciceAdapter::FSI::FluidStructureInteraction::readConfig(const YAML::Node adapterConfig)
{
    /* TODO: Read the solver type, if needed.
    /  If you want to determine it automatically, implement a method
    /  as in CHT/CHT.C
    */

    /* TODO: Read the names of any needed fields and parameters.
    * Include the force here?
    */

    // Read the name of the pointDisplacement field (if different)
    if (adapterConfig["namePointDisplacement"])
    {
        namePointDisplacement_ = adapterConfig["namePointDisplacement"].as<std::string>();
    }
    DEBUG(adapterInfo("    pointDisplacement field name : " + namePointDisplacement_));

    if (adapterConfig["namePointVelocity"])
	{
		namePointVelocity_ = adapterConfig["namePointVelocity"].as<std::string>();
	}
	DEBUG(adapterInfo("    pointVelocity field name : " + namePointVelocity_));

	//~ if (adapterConfig["velocityStabilizationFactor"])
	//~ {
		//~ velocStabFac_ = adapterConfig["velocityStabilizationFactor"].as<double>();
	//~ }
	//~ DEBUG(adapterInfo("     velocStabFac_ value is : " + std::to_string(velocStabFac_)));

	//~ if (adapterConfig["forceStabilizationFactor"])
	//~ {
		//~ forceStabFac_ = adapterConfig["forceStabilizationFactor"].as<double>();
	//~ }
	//~ DEBUG(adapterInfo("     forceStabFac_ value is : " + std::to_string(forceStabFac_)));

    return true;
}

void preciceAdapter::FSI::FluidStructureInteraction::addWriters(std::string dataName, Interface * interface)
{
    if (dataName.find("Force") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new Force(mesh_, runTime_.timeName(), nameRho_, nameU_, pRef_) //, forceStabFac_)
        );
        DEBUG(adapterInfo("Added writer: Force."));
    }
    else if (dataName.find("DisplacementDelta") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new DisplacementDelta(mesh_, namePointDisplacement_)
        );
        DEBUG(adapterInfo("Added writer: DisplacementDelta."));
    }
    else if (dataName.find("Displacement") == 0)
    {
        interface->addCouplingDataWriter
        (
            dataName,
            new Displacement(mesh_, namePointDisplacement_)
        );
        DEBUG(adapterInfo("Added writer: Displacement."));
    }
    //~ else if (dataName.find("VelocityStabilized") == 0)
	//~ {
		//~ interface->addCouplingDataWriter
		//~ (
			//~ dataName,
			//~ new VelocityStabilized(mesh_, namePointVelocity_, velocStabFac_)
		//~ );
		//~ DEBUG(adapterInfo("Added writer: VelocityStabilized."));
	//~ }
    else if (dataName.find("Velocit") == 0)
	{
		interface->addCouplingDataWriter
		(
			dataName,
			new Velocity(mesh_, namePointVelocity_)
		);
		DEBUG(adapterInfo("Added writer: Velocity."));
	}

    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a reader below).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.
}

void preciceAdapter::FSI::FluidStructureInteraction::addReaders(std::string dataName, Interface * interface)
{
	Info << "DATA NAMEE: " << dataName << endl;
    if (dataName.find("Force") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new Force(mesh_, runTime_.timeName(), nameRho_, nameU_, pRef_) //, forceStabFac_)
        );
        DEBUG(adapterInfo("Added reader: Force."));
    }
    else if (dataName.find("DisplacementDelta") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new DisplacementDelta(mesh_, namePointDisplacement_)
        );
        DEBUG(adapterInfo("Added reader: DisplacementDelta."));
    }
    else if (dataName.find("Displacement") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            new Displacement(mesh_, namePointDisplacement_)
        );
        DEBUG(adapterInfo("Added reader: Displacement."));
    }
    //~ else if (dataName.find("VelocityStabilized") == 0)
	//~ {
		//~ interface->addCouplingDataReader
		//~ (
			//~ dataName,
			//~ new VelocityStabilized(mesh_, namePointVelocity_, velocStabFac_)
		//~ );
		//~ DEBUG(adapterInfo("Added reader: VelocityStabilized."));
	//~ }
    else if (dataName.find("Velocit") == 0)
	{
		interface->addCouplingDataReader
		(
			dataName,
			new Velocity(mesh_, namePointVelocity_)
		);
		DEBUG(adapterInfo("Added reader: Velocity."));
	}
    // NOTE: If you want to couple another variable, you need
    // to add your new coupling data user as a coupling data
    // writer here (and as a writer above).
    // The argument of the dataName.compare() needs to match
    // the one provided in the adapter's configuration file.
}
