#include "FSI.H"

#include "Utilities.H"

using namespace Foam;

preciceAdapter::FSI::FluidStructureInteraction::FluidStructureInteraction
(
    const Foam::fvMesh& mesh,
    const Foam::Time& runTime
)
:
mesh_(mesh),
runTime_(runTime)
{}

bool preciceAdapter::FSI::FluidStructureInteraction::configure(const IOdictionary& adapterConfig)
{
    DEBUG(adapterInfo("Configuring the FSI module..."));

    // Read the FSI-specific options from the adapter's configuration file
    if (!readConfig(adapterConfig)) return false;

    // NOTE: If you want to add a new solver type, which you can manually
    // specify in the configuration, add it here. See also the methods
    // addWriters() and addReaders().
    // Check the solver type and determine it if needed
    if (
        solverType_.compare("compressible") == 0 ||
        solverType_.compare("incompressible") == 0 ||
        solverType_.compare("basic") == 0
    )
    {
        DEBUG(adapterInfo("Known solver type: " + solverType_));
    }
    else if (solverType_.compare("none") == 0)
    {
        DEBUG(adapterInfo("Determining the solver type..."));
        solverType_ = determineSolverType();
    }
    else
    {
        DEBUG(adapterInfo("Unknown solver type. Determining the solver type..."));
        solverType_ = determineSolverType();
    }

    return true;
}

bool preciceAdapter::FSI::FluidStructureInteraction::readConfig(const IOdictionary& adapterConfig)
{
    const dictionary FSIdict = adapterConfig.subOrEmptyDict("FSI");
  
    // Read the solver type (if not specified, it is determined automatically)
    solverType_ = FSIdict.lookupOrDefault<word>("solverType", "");
    DEBUG(adapterInfo("    user-defined solver type : " + solverType_));

    /* TODO: Read the names of any needed fields and parameters.
    * Include the force here?
    */

    // Read the name of the pointDisplacement field (if different)
    namePointDisplacement_ = FSIdict.lookupOrDefault<word>("namePointDisplacement", "pointDisplacement");
    namePointVelocity_ = FSIdict.lookupOrDefault<word>("namePointVelocity", "pointMotionU");
    DEBUG(adapterInfo("    pointDisplacement field name : " + namePointDisplacement_));

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

// NOTE: This is exactly the same as in the CHT module.
std::string preciceAdapter::FSI::FluidStructureInteraction::determineSolverType()
{
    // NOTE: When coupling a different variable, you may want to
    // add more cases here. Or you may provide the solverType in the config.

    std::string solverType;

    // Determine the solver type: Compressible, Incompressible or Basic.
    // Look for the files transportProperties, turbulenceProperties,
    // and thermophysicalProperties
    bool transportPropertiesExists = false;
    bool turbulencePropertiesExists = false;
    bool thermophysicalPropertiesExists = false;

    if (mesh_.foundObject<IOdictionary>("transportProperties"))
    {
        transportPropertiesExists = true;
        DEBUG(adapterInfo("Found the transportProperties dictionary."));
    }
    else
    {
        DEBUG(adapterInfo("Did not find the transportProperties dictionary."));
    }

    
    if (mesh_.foundObject<IOdictionary>(turbulenceModel::propertiesName))
    {
        turbulencePropertiesExists = true;
        DEBUG(adapterInfo("Found the " + turbulenceModel::propertiesName
            + " dictionary."));
    }
    else
    {
        DEBUG(adapterInfo("Did not find the " + turbulenceModel::propertiesName
            + " dictionary."));
    }
    

    if (mesh_.foundObject<IOdictionary>("thermophysicalProperties"))
    {
        thermophysicalPropertiesExists = true;
        DEBUG(adapterInfo("Found the thermophysicalProperties dictionary."));
    }
    else
    {
        DEBUG(adapterInfo("Did not find the thermophysicalProperties dictionary."));
    }

    if (turbulencePropertiesExists)
    {
        if (thermophysicalPropertiesExists)
        {
            solverType = "compressible";
            DEBUG(adapterInfo("This is a compressible flow solver, "
                "as turbulence and thermophysical properties are provided."));
        }
        else if (transportPropertiesExists)
        {
            solverType = "incompressible";
            DEBUG(adapterInfo("This is an incompressible flow solver, "
            "as turbulence and transport properties are provided."));
        }
        else
        {
            adapterInfo("Could not determine the solver type, or this is not "
            "a compatible solver: although turbulence properties are provided, "
            "neither transport or thermophysical properties are provided.",
            "error");
        }
    }
    else
    {
        if (transportPropertiesExists)
        {
            solverType = "basic";
            DEBUG(adapterInfo("This is a basic solver, as transport properties "
            "are provided, while turbulence or transport properties are not "
            "provided."));
        }
        else
        {
            adapterInfo("Could not determine the solver type, or this is not a "
            "compatible solver: neither transport, nor turbulence properties "
            "are provided.",
            "error");
        }
    }

    return solverType;
}


void preciceAdapter::FSI::FluidStructureInteraction::addWriters(std::string dataName, Interface * interface)
{
    if (dataName.find("Force") == 0)
    {        
            interface->addCouplingDataWriter
            (
                dataName,
                new Force(mesh_, runTime_.timeName(), solverType_) /* TODO: Add any other arguments here */
		//new Force(mesh_, runTime_.timeName(), nameRho_, nameU_, pRef_) //, forceStabFac_)   
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
    if (dataName.find("Force") == 0)
    {
        interface->addCouplingDataReader
        (
            dataName,
            // new Force(mesh_, runTime_.timeName(), nameRho_, nameU_, pRef_) //, forceStabFac_)
            new Force(mesh_, runTime_.timeName(), solverType_) /* TODO: Add any other arguments here */
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
