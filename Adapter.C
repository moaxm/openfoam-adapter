#include "Adapter.H"
#include "Interface.H"
#include "Utilities.H"

#include "IOstreams.H"
#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"

using namespace Foam;

preciceAdapter::Adapter::Adapter(const Time& runTime, const fvMesh& mesh)
:
runTime_(runTime),
mesh_(mesh)
{
    adapterInfo("The preciceAdapter was loaded.", "info");

    #ifdef ADAPTER_DEBUG_MODE
        Info<< "Registered objects: " << mesh_.names() << endl;
    #endif

    return;
}


bool preciceAdapter::Adapter::configFileRead()
{

    // preCICE participant name
    participantName_ = "Fluid1";
    DEBUG(adapterInfo("  participant : " + participantName_));

    // preCICE configuration file name
    preciceConfigFilename_ = "precice-config.xml";
    DEBUG(adapterInfo("  precice-config-file : " + preciceConfigFilename_));

    DEBUG(adapterInfo("  interfaces : "));
        struct InterfaceConfig interfaceConfig;
        interfaceConfig.meshName = "Fluid1-Mesh";
        DEBUG(adapterInfo("  - mesh      : " + interfaceConfig.meshName));

        interfaceConfig.locationsType = "faceCenters";
        DEBUG(adapterInfo("    locations : " + interfaceConfig.locationsType));

        DEBUG(adapterInfo("    patches   : "));
            interfaceConfig.patchNames.push_back("outlet");
            DEBUG(adapterInfo("       outlet"));

        // write-data
            DEBUG(adapterInfo("    write-data : "));
                interfaceConfig.writeData.push_back("Velocity");
                interfaceConfig.writeData.push_back("PressureGradient");
                DEBUG(adapterInfo("      Velocity, PressureGradient"));

        // read-data
                interfaceConfig.readData.push_back("Pressure");
                DEBUG(adapterInfo("      Pressure"));

        interfacesConfig_.push_back(interfaceConfig);


        FFenabled_ = true;
    DEBUG(adapterInfo("    FF module enabled : " + std::to_string(FFenabled_)));

    // NOTE: set the switch for your new module here

    if (FFenabled_)
    {
        FF_ = new FF::FluidFluid(mesh_);
    }

    // TODO: Loading modules should be implemented in more general way,
    // in order to avoid code duplication. See issue #16 on GitHub.

    return true;
}

void preciceAdapter::Adapter::configure()
{
    // Read the adapter's configuration file
    if (!configFileRead())
    {
        // This method is called from the functionObject's read() method,
        // which is called by the Foam::functionObjectList::read() method.
        // All the exceptions triggered in this method are caught as
        // warnings and the simulation continues simply without the
        // functionObject. However, we want the simulation to exit with an
        // error in case something is wrong. We store the information that
        // there was an error and it will be handled by the first call to
        // the functionObject's execute(), which can throw errors normally.
        errorsInConfigure = true;

        return;
    }

try{
    // Check the timestep type (fixed vs adjustable)
    DEBUG(adapterInfo("Checking the timestep type (fixed vs adjustable)..."));
    adjustableTimestep_ = runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);

    if (adjustableTimestep_) {
        DEBUG(adapterInfo("  Timestep type: adjustable."));
    } else {
        DEBUG(adapterInfo("  Timestep type: fixed."));
    }

    // Initialize preCICE
    DEBUG(adapterInfo("Creating the preCICE solver interface..."));
    DEBUG(adapterInfo("  Number of processes: " + std::to_string(Pstream::nProcs())));
    DEBUG(adapterInfo("  MPI rank: " + std::to_string(Pstream::myProcNo())));
    precice_ = new precice::SolverInterface(participantName_, Pstream::myProcNo(), Pstream::nProcs());
    DEBUG(adapterInfo("  preCICE solver interface was created."));

    DEBUG(adapterInfo("Configuring preCICE..."));
    precice_->configure(preciceConfigFilename_);
    DEBUG(adapterInfo("  preCICE was configured."));

    // Create interfaces
    DEBUG(adapterInfo("Creating interfaces..."));
    for (uint i = 0; i < interfacesConfig_.size(); i++)
    {
        Interface * interface = new Interface(*precice_, mesh_, interfacesConfig_.at(i).meshName, interfacesConfig_.at(i).locationsType, interfacesConfig_.at(i).patchNames);
        interfaces_.push_back(interface);
        DEBUG(adapterInfo("Interface created on mesh " + interfacesConfig_.at(i).meshName));

        DEBUG(adapterInfo("Adding coupling data writers..."));
        for (uint j = 0; j < interfacesConfig_.at(i).writeData.size(); j++)
        {
            std::string dataName = interfacesConfig_.at(i).writeData.at(j);

            // Add FF-related coupling data writers
            if (FFenabled_)
            {
                FF_->addWriters(dataName, interface);
            }

            // NOTE: Add any coupling data writers for your module here.
        } // end add coupling data writers

        DEBUG(adapterInfo("Adding coupling data readers..."));
        for (uint j = 0; j < interfacesConfig_.at(i).readData.size(); j++)
        {
            std::string dataName = interfacesConfig_.at(i).readData.at(j);

            // Add FF-related coupling data readers
            if (FFenabled_)
            {
                FF_->addReaders(dataName, interface);
            }

            // NOTE: Add any coupling data readers for your module here.
        } // end add coupling data readers

        // Create the interface's data buffer
        interface->createBuffer();
    }

    // Initialize preCICE and exchange the first coupling data
    initialize();

    // Read the received coupling data
    readCouplingData();

    // If checkpointing is required, specify the checkpointed fields
    // and write the first checkpoint
    if (isWriteCheckpointRequired())
    {
        checkpointing_ = true;

        // Setup the checkpointing (find and add fields to checkpoint)
        if (!disableCheckpointing_)
        {
            setupCheckpointing();
        }

        // Write checkpoint (for the first iteration)
        writeCheckpoint();
        fulfilledWriteCheckpoint();
    }

    // Adjust the timestep for the first iteration, if it is fixed
    if (!adjustableTimestep_)
    {
        adjustSolverTimeStep();
    }

    // If the solver tries to end before the coupling is complete,
    // e.g. because the solver's endTime was smaller or (in implicit
    // coupling) equal with the max-time specified in preCICE,
    // problems may occur near the end of the simulation,
    // as the function object may be called only once near the end.
    // See the implementation of Foam::Time::run() for more details.
    // To prevent this, we set the solver's endTime to "infinity"
    // and let only preCICE control the end of the simulation.
    // This has the side-effect of not triggering the end() method
    // in any function object normally. Therefore, we trigger it
    // when preCICE dictates to stop the coupling.
    // However, the user can disable this behavior in the configuration.
    if (preventEarlyExit_)
    {
        adapterInfo
        (
            "Setting the solver's endTime to infinity to prevent early exits. "
            "Only preCICE will control the simulation's endTime. "
            "Any functionObject's end() method will be triggered by the adapter. "
            "You may disable this behavior in the adapter's configuration.",
            "info"
       );
        const_cast<Time&>(runTime_).setEndTime(GREAT);
    }

} catch (const Foam::error &e) {
    adapterInfo(e.message(), "info");
    errorsInConfigure = true;
}

    return;
}

void preciceAdapter::Adapter::execute()
{
    if (errorsInConfigure)
    {
        // Handle any errors during configure().
        // See the comments in configure() for details.
        adapterInfo
        (
            "There was a problem while configuring the adapter. "
            "See the log for details.",
            "error"
       );
    }

    // The solver has already solved the equations for this timestep.
    // Now call the adapter's methods to perform the coupling.

    // TODO add a function which checks if all fields are checkpointed. 
    // if (ncheckpointed is nregisterdobjects. )

    // Write the coupling data in the buffer
    writeCouplingData();

    // Advance preCICE
    advance();

    // Read checkpoint if required
    if (isReadCheckpointRequired())
    {
        readCheckpoint();
        fulfilledReadCheckpoint();
    }

    // Read the received coupling data from the buffer
    readCouplingData();

    // Adjust the timestep, if it is fixed
    if (!adjustableTimestep_)
    {
        adjustSolverTimeStep();
    }

    // Write checkpoint if required
    if (isWriteCheckpointRequired())
    {
        writeCheckpoint();
        fulfilledWriteCheckpoint();
    }

    // As soon as OpenFOAM writes the results, it will not try to write again
    // if the time takes the same value again. Therefore, during an implicit
    // coupling, we write again when the coupling timestep is complete.
    // Check the behavior e.g. by using watch on a result file:
    //     watch -n 0.1 -d ls --full-time Fluid/0.01/T.gz
    if (checkpointing_ && isCouplingTimestepComplete())
    {
        // Check if the time directory already exists
        // (i.e. the solver wrote results that need to be updated)
        if (runTime_.timePath().type() == fileName::DIRECTORY)
        {
            adapterInfo
            (
                "The coupling timestep completed. "
                "Writing the updated results.",
                "info"
           );
            const_cast<Time&>(runTime_).writeNow();
        }
    }

    // If the coupling is not going to continue, tear down everything
    // and stop the simulation.
    if (!isCouplingOngoing())
    {
        adapterInfo("The coupling completed.", "info");

        // Finalize the preCICE solver interface and delete data
        finalize();

        // Tell OpenFOAM to stop the simulation.
        // Set the solver's endTime to now. The next evaluation of
        // runTime.run() will be false and the solver will exit.
        const_cast<Time&>(runTime_).setEndTime(runTime_.value());

        if (preventEarlyExit_)
        {
            adapterInfo
            (
                "The simulation was ended by preCICE. "
                "Calling the end() methods of any functionObject explicitly.",
                "info"
           );
            const_cast<Time&>(runTime_).functionObjects().end();
        }
    }

    return;
}

void preciceAdapter::Adapter::adjustTimeStep()
{
    adjustSolverTimeStep();

    return;
}

void preciceAdapter::Adapter::readCouplingData()
{
    DEBUG(adapterInfo("Reading coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->readCouplingData();
    }

    return;
}

void preciceAdapter::Adapter::writeCouplingData()
{
    DEBUG(adapterInfo("Writing coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->writeCouplingData();
    }

    return;
}

void preciceAdapter::Adapter::initialize()
{
    DEBUG(adapterInfo("Initalizing the preCICE solver interface..."));
    timestepPrecice_ = precice_->initialize();

    preciceInitialized_ = true;

    if (precice_->isActionRequired(precice::constants::actionWriteInitialData()))
    {
        writeCouplingData();
        precice_->fulfilledAction(precice::constants::actionWriteInitialData());
    }

    DEBUG(adapterInfo("Initializing preCICE data..."));
    precice_->initializeData();

    adapterInfo("preCICE was configured and initialized", "info");

    return;
}

void preciceAdapter::Adapter::finalize()
{
    if (NULL != precice_ && preciceInitialized_ && !isCouplingOngoing())
    {
        DEBUG(adapterInfo("Finalizing the preCICE solver interface..."));

        // Finalize the preCICE solver interface
        precice_->finalize();

        preciceInitialized_ = false;

        // Delete the solver interface and all the related data
        teardown();
    }
    else
    {
        adapterInfo("Could not finalize preCICE.", "error");
    }

    return;
}

void preciceAdapter::Adapter::advance()
{
    DEBUG(adapterInfo("Advancing preCICE..."));

    timestepPrecice_ = precice_->advance(timestepSolver_);

    return;
}

void preciceAdapter::Adapter::adjustSolverTimeStep()
{
    DEBUG(adapterInfo("Adjusting the solver's timestep..."));

    // The timestep size that the solver has determined that it wants to use
    double timestepSolverDetermined;

    /* In this method, the adapter overwrites the timestep used by OpenFOAM.
       If the timestep is not adjustable, OpenFOAM will not try to re-estimate
       the timestep or read it again from the controlDict. Therefore, store
       the value that the timestep has is the beginning and try again to use this
       in every iteration.
       // TODO Treat also the case where the user modifies the timestep
       // in the controlDict during the simulation.
    */

    // Is the timestep adjustable or fixed?
    if (!adjustableTimestep_)
    {
        // Have we already stored the timestep?
        if (!useStoredTimestep_)
        {
            // Show a warning if runTimeModifiable is set
            if (runTime_.runTimeModifiable())
            {
                adapterInfo
                (
                    "You have enabled 'runTimeModifiable' in the "
                    "controlDict. The preciceAdapter does not yet "
                    "fully support this functionality when "
                    "'adjustableTimestep' is not enabled. "
                    "If you modify the 'deltaT' in the controlDict "
                    "during the simulation, it will not be updated.",
                    "warning"
               );
            }

            // Store the value
            timestepStored_ = runTime_.deltaT().value();

            // Ok, we stored it once, we will use this from now on
            useStoredTimestep_ = true;
        }

        // Use the stored timestep as the determined solver's timestep
        timestepSolverDetermined = timestepStored_;
    }
    else
    {
        // The timestep is adjustable, so OpenFOAM will modify it
        // and therefore we can use the updated value
        timestepSolverDetermined = runTime_.deltaT().value();
    }

    /* If the solver tries to use a timestep smaller than the one determined
       by preCICE, that means that the solver is trying to subcycle.
       This may not be allowed by the user.
       If the solver tries to use a bigger timestep, then it needs to use
       the same timestep as the one determined by preCICE.
    */

    if (timestepSolverDetermined < timestepPrecice_)
    {
        if (!subcyclingAllowed_)
        {
            adapterInfo
            (
                "The solver's timestep cannot be smaller than the "
                "coupling timestep, because subcycling is disabled. ",
                "error"
           );
        }
        else
        {
            // Add a bool 'subCycling = true' which is checked in the storeMeshPoints() function. 
            adapterInfo
            (
                "The solver's timestep is smaller than the "
                "coupling timestep. Subcycling...",
                "info"
           );
            timestepSolver_ = timestepSolverDetermined;
        }
    }
    else if (timestepSolverDetermined > timestepPrecice_)
    {
        adapterInfo
        (
            "The solver's timestep cannot be larger than the coupling timestep."
            " Adjusting from " +
            std::to_string(timestepSolverDetermined) +
            " to " +
            std::to_string(timestepPrecice_),
            "warning"
       );
        timestepSolver_ = timestepPrecice_;
    }
    else
    {
        DEBUG(adapterInfo("The solver's timestep is the same as the "
                            "coupling timestep."));
        timestepSolver_ = timestepPrecice_;
    }

    // Update the solver's timestep (but don't trigger the adjustDeltaT(),
    // which also triggers the functionObject's adjustTimeStep())
    // TODO: Keep this in mind if any relevant problem appears.
    const_cast<Time&>(runTime_).setDeltaT(timestepSolver_, false);

    return;
}

bool preciceAdapter::Adapter::isCouplingOngoing()
{
    bool isCouplingOngoing = false;

    // If the coupling ends before the solver ends,
    // the solver would try to access this method again,
    // giving a segmentation fault if precice_
    // was not available.
    if (NULL != precice_)
    {
        isCouplingOngoing = precice_->isCouplingOngoing();
    }

    return isCouplingOngoing;
}

bool preciceAdapter::Adapter::isCouplingTimestepComplete()
{
    return precice_->isTimestepComplete();
}

bool preciceAdapter::Adapter::isReadCheckpointRequired()
{
    return precice_->isActionRequired(precice::constants::actionReadIterationCheckpoint());
}

bool preciceAdapter::Adapter::isWriteCheckpointRequired()
{
    return precice_->isActionRequired(precice::constants::actionWriteIterationCheckpoint());
}

void preciceAdapter::Adapter::fulfilledReadCheckpoint()
{
    precice_->fulfilledAction(precice::constants::actionReadIterationCheckpoint());

    return;
}

void preciceAdapter::Adapter::fulfilledWriteCheckpoint()
{
    precice_->fulfilledAction(precice::constants::actionWriteIterationCheckpoint());

    return;
}

void preciceAdapter::Adapter::storeCheckpointTime()
{
    couplingIterationTimeIndex_ = runTime_.timeIndex();
    couplingIterationTimeValue_ = runTime_.value();
    DEBUG(adapterInfo("Stored time value t = " + std::to_string(runTime_.value())));

    return;
}

void preciceAdapter::Adapter::reloadCheckpointTime()
{
    const_cast<Time&>(runTime_).setTime(couplingIterationTimeValue_, couplingIterationTimeIndex_);
    // TODO also reset the current iteration?!
    DEBUG(adapterInfo("Reloaded time value t = " + std::to_string(runTime_.value())));

    return;
}

void preciceAdapter::Adapter::storeMeshPoints()
{
    DEBUG(adapterInfo("Storing mesh points..."));
    // TODO: In foam-extend, we would need "allPoints()". Check if this gives the same data.
    meshPoints_ = mesh_.points();
    oldMeshPoints_ = mesh_.oldPoints();

    /*  
    // TODO  This is only required for subcycling. It should not be called when not subcycling!!
    // Add a bool 'subcycling' which can be evaluated every timestep. 
    if ( !oldVolsStored && mesh_.foundObject<volScalarField::Internal>("V00") ) // For Ddt schemes which use one previous timestep
    {  
        setupMeshVolCheckpointing();
        oldVolsStored = true;
    }
    // Update any volume fields from the buffer to the checkpointed values (if already exists.)
    */

    DEBUG(adapterInfo("Stored mesh points."));    
    if (mesh_.moving())
    {
        if (!meshCheckPointed)
        {
            // Set up the checkpoint for the mesh flux: meshPhi
            setupMeshCheckpointing();
            meshCheckPointed = true;
        }
        writeMeshCheckpoint();
        writeVolCheckpoint(); // Does not write anything unless subcycling.
    }
}

void preciceAdapter::Adapter::reloadMeshPoints()
{
    // In Foam::polyMesh::movePoints.
    // TODO: The function movePoints overwrites the pointer to the old mesh. 
    // Therefore, if you revert the mesh, the oldpointer will be set to the points, which are the new values. 
    DEBUG(adapterInfo("Moving mesh points to their previous locations..."));

    // TODO
    // Switch oldpoints on for pure physics. (is this required?). Switch off for better mesh deformation capabilities?    
    // const_cast<pointField&>(mesh_.points()) = oldMeshPoints_;
    const_cast<fvMesh&>(mesh_).movePoints(meshPoints_);
   
    DEBUG(adapterInfo("Moved mesh points to their previous locations."));
    
    // TODO The if statement can be removed in this case, but it is still included for clarity 
    if ( meshCheckPointed )
    {   
        readMeshCheckpoint();
    }

    /*  // TODO This part should only be used when sybcycling. See the description in 'storeMeshPoints()'
        // The if statement can be removed in this case, but it is still included for clarity
    if ( oldVolsStored )
    {       
        readVolCheckpoint();
    }
    */ 
}

void preciceAdapter::Adapter::setupMeshCheckpointing()
{
    // The other mesh <type>Fields: 
    //      C
    //      Cf
    //      Sf
    //      magSf
    //      delta
    // are updated by the function fvMesh::movePoints. Only the meshPhi needs checkpointing. 
    DEBUG(adapterInfo("Creating a list of the mesh checkpointed fields..."));

    // Add meshPhi to the checkpointed fields
    addMeshCheckpointField
    (
        const_cast<surfaceScalarField&>
        (
            mesh_.phi()
        )
    );
    #ifdef ADAPTER_DEBUG_MODE
    adapterInfo
    (
        "Added " + mesh_.phi().name() +
        " in the list of checkpointed fields."
    );
    #endif
    
}

void preciceAdapter::Adapter::setupMeshVolCheckpointing()
{
    DEBUG(adapterInfo("Creating a list of the mesh volume checkpointed fields..."));
    // Add the V0 and the V00 to the list of checkpointed fields. 
    // For V0
    addVolCheckpointField
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V0()
        )
    );
    #ifdef ADAPTER_DEBUG_MODE
    adapterInfo
    (
        "Added " + mesh_.V0().name() +
        " in the list of checkpointed fields."
    );
    #endif
    // For V00
    addVolCheckpointField
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V00()
        )
    );
    #ifdef ADAPTER_DEBUG_MODE
    adapterInfo
    (
        "Added " + mesh_.V00().name() +
        " in the list of checkpointed fields."
    );
    #endif

    // Also add the buffer fields.
    // TODO For V0
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V0()
        )
    ); */
    #ifdef ADAPTER_DEBUG_MODE
    adapterInfo
    (
        "Added " + mesh_.V0().name() +
        " in the list of buffer checkpointed fields."
    );
    #endif
    // TODO For V00
    /* addVolCheckpointFieldBuffer
    (
        const_cast<volScalarField::Internal&>
        (
            mesh_.V00()
        )
    );*/
    #ifdef ADAPTER_DEBUG_MODE
    adapterInfo
    (
        "Added " + mesh_.V00().name() +
        " in the list of buffer checkpointed fields."
    );
    #endif
}


void preciceAdapter::Adapter::setupCheckpointing()
{
    // Add fields in the checkpointing list
    DEBUG(adapterInfo("Creating a list of checkpointed fields..."));

    /* Find and add all the registered objects in the mesh_
       of type volScalarField
    */

    // Print the available objects of type volScalarField
    DEBUG(adapterInfo("Available objects of type volScalarField : "));
    #ifdef ADAPTER_DEBUG_MODE
        Info << mesh_.lookupClass<volScalarField>() << nl << nl;
    #endif

    wordList objectNames_ = mesh_.lookupClass<volScalarField>().toc();

    forAll(objectNames_, i)
    {
        if (mesh_.foundObject<volScalarField>(objectNames_[i]))
        {
            addCheckpointField
            (
                const_cast<volScalarField&>
                (
                    mesh_.lookupObject<volScalarField>(objectNames_[i])
                )
            );

            #ifdef ADAPTER_DEBUG_MODE
            adapterInfo
            (
                "Added " + objectNames_[i] +
                " in the list of checkpointed fields."
           );
            #endif

            // TODO: Known bug, see readCheckpoint()
            if ("epsilon" == objectNames_[i])
            {
                DEBUG(adapterInfo("Known bug: after reading a checkpoint, "
                        "the boundaries for epsilon will not be corrected.",
                        "warning"));
            }
        }
        else
        {
            adapterInfo("Could not checkpoint " + objectNames_[i], "warning");
        }
    }

    /* Find and add all the registered objects in the mesh_
       of type volVectorField
    */

    // Print the available objects of type volVectorField
    DEBUG(adapterInfo("Available objects of type volVectorField : "));
    #ifdef ADAPTER_DEBUG_MODE
        Info << mesh_.lookupClass<volVectorField>() << nl << nl;
    #endif

    objectNames_ = mesh_.lookupClass<volVectorField>().toc();

    forAll(objectNames_, i)
    {
        if (mesh_.foundObject<volVectorField>(objectNames_[i]))
        {
            addCheckpointField
            (
                const_cast<volVectorField&>
                (
                    mesh_.lookupObject<volVectorField>(objectNames_[i])
                )
            );

            #ifdef ADAPTER_DEBUG_MODE
            adapterInfo
            (
                "Added " + objectNames_[i] +
                " in the list of checkpointed fields."
           );
            #endif
        }
        else
        {
            adapterInfo("Could not checkpoint " + objectNames_[i], "warning");
        }
    }

    // Print the available objects of type surfaceScalarField
    DEBUG(adapterInfo("Available objects of type surfaceScalarField : "));
    #ifdef ADAPTER_DEBUG_MODE
        Info << mesh_.lookupClass<surfaceScalarField>() << nl << nl;
    #endif

    /* Find and add all the registered objects in the mesh_
       of type surfaceScalarField
    */
    objectNames_ = mesh_.lookupClass<surfaceScalarField>().toc();

    forAll(objectNames_, i)
    {
        if (mesh_.foundObject<surfaceScalarField>(objectNames_[i]))
        {
            addCheckpointField
            (
                const_cast<surfaceScalarField&>
                (
                    mesh_.lookupObject<surfaceScalarField>(objectNames_[i])
               )
           );

            #ifdef ADAPTER_DEBUG_MODE
            adapterInfo
            (
                "Added " + objectNames_[i] +
                " in the list of checkpointed fields."
           );
            #endif
        }
        else
        {
            adapterInfo("Could not checkpoint " + objectNames_[i], "warning");
        }
    }

    /* Find and add all the registered objects in the mesh_
       of type surfaceVectorField
    */

    // Print the available objects of type surfaceVectorField
    DEBUG(adapterInfo("Available objects of type surfaceVectorField : "));
    #ifdef ADAPTER_DEBUG_MODE
        Info << mesh_.lookupClass<surfaceVectorField>() << nl << nl;
    #endif

    objectNames_ = mesh_.lookupClass<surfaceVectorField>().toc();

    forAll(objectNames_, i)
    {
        if (mesh_.foundObject<surfaceVectorField>(objectNames_[i]))
        {
            addCheckpointField
            (
                const_cast<surfaceVectorField&>
                (
                    mesh_.lookupObject<surfaceVectorField>(objectNames_[i])
               )
           );

            #ifdef ADAPTER_DEBUG_MODE
            adapterInfo
            (
                "Added " + objectNames_[i] +
                " in the list of checkpointed fields."
           );
            #endif
        }
        else
        {
            adapterInfo("Could not checkpoint " + objectNames_[i], "warning");
        }
    }

    return;
}


// All mesh checkpointed fields
void preciceAdapter::Adapter::addMeshCheckpointField(surfaceScalarField & field)
{
    surfaceScalarField * copy = new surfaceScalarField(field);
    meshSurfaceScalarFields_.push_back(&field);
    meshSurfaceScalarFieldCopies_.push_back(copy);
    return;
}

void preciceAdapter::Adapter::addMeshCheckpointField(surfaceVectorField & field)
{
    surfaceVectorField * copy = new surfaceVectorField(field);
    meshSurfaceVectorFields_.push_back(&field);
    meshSurfaceVectorFieldCopies_.push_back(copy);
    return;
}

void preciceAdapter::Adapter::addMeshCheckpointField(volVectorField & field)
{
    volVectorField * copy = new volVectorField(field);
    meshVolVectorFields_.push_back(&field);
    meshVolVectorFieldCopies_.push_back(copy);
    return;
}

// TODO Internal field for the V0 (volume old) and V00 (volume old-old) fields
void preciceAdapter::Adapter::addVolCheckpointField(volScalarField::Internal & field)
{
    volScalarField::Internal * copy = new volScalarField::Internal(field);
    volScalarInternalFields_.push_back(&field);
    volScalarInternalFieldCopies_.push_back(copy);
    return;
}

// All checkpointed fields
void preciceAdapter::Adapter::addCheckpointField(volScalarField & field)
{
    volScalarField * copy = new volScalarField(field);
    volScalarFields_.push_back(&field);
    volScalarFieldCopies_.push_back(copy);
    return;
}

void preciceAdapter::Adapter::addCheckpointField(volVectorField & field)
{
    volVectorField * copy = new volVectorField(field);
    volVectorFields_.push_back(&field);
    volVectorFieldCopies_.push_back(copy);
    return;
}

void preciceAdapter::Adapter::addCheckpointField(surfaceScalarField & field)
{
    surfaceScalarField * copy = new surfaceScalarField(field);
    surfaceScalarFields_.push_back(&field);
    surfaceScalarFieldCopies_.push_back(copy);
    return;
}

void preciceAdapter::Adapter::addCheckpointField(surfaceVectorField & field)
{
    surfaceVectorField * copy = new surfaceVectorField(field);
    surfaceVectorFields_.push_back(&field);
    surfaceVectorFieldCopies_.push_back(copy);
    return;
}

// NOTE: Add here methods to add other object types to checkpoint, if needed.

void preciceAdapter::Adapter::readCheckpoint()
{
    // TODO: To increase efficiency: only the oldTime() fields of the quantities which are used in the time 
    //  derivative are necessary. (In general this is only the velocity). Also old information of the mesh
    //  is required. 
    //  Therefore, loading the oldTime() and oldTime().oldTime() fields for the other fields can be excluded 
    //  for efficiency.  
    DEBUG(adapterInfo("Reading a checkpoint..."));

    // Reload the runTime
    reloadCheckpointTime();

    // Reload all the fields of type volScalarField
    for (uint i = 0; i < volScalarFields_.size(); i++)
    {
        // Load the volume field 
        *(volScalarFields_.at(i)) == *(volScalarFieldCopies_.at(i));
        // TODO: Do we need this?
        // *(volScalarFields_.at(i))->boundaryField() = *(volScalarFieldCopies_.at(i))->boundaryField();

        int nOldTimes(volScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            volScalarFields_.at(i)->oldTime() == volScalarFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            volScalarFields_.at(i)->oldTime().oldTime() == volScalarFieldCopies_.at(i)->oldTime().oldTime();
        }

        // Evaluate the boundaries, if supported
        if (evaluateBoundaries_)
        {
            try{
                // TODO Check if these fields require adding besides only epsilon.
                // (from Max Mueller's fork)
                /* 
                if ( ("epsilon"  != volScalarFields_.at(i)->name()) &&
                ("epsilon_0"!= volScalarFields_.at(i)->name()) &&
                ("omega"    != volScalarFields_.at(i)->name()) &&
                ("omega_0"  != volScalarFields_.at(i)->name()) &&
                ("cellDisplacementx"!=volScalarFields_.at(i)->name()) &&
                ("cellDisplacementy"!=volScalarFields_.at(i)->name()) &&
                ("cellDisplacementz"!=volScalarFields_.at(i)->name()))
                */
                if ("epsilon" != volScalarFields_.at(i)->name())
                {
                    volScalarFields_.at(i)->correctBoundaryConditions();
                }
                // TODO: Known bug: cannot find "volScalarField::Internal kEpsilon:G"
                // Currently it is skipped. Before it was not corrected at all.
                // A warning for this is thrown when adding epsilon to the checkpoint.
            } catch (Foam::error) {
                DEBUG(adapterInfo("Could not evaluate the boundary for" + volScalarFields_.at(i)->name(), "warning"));
            }
        }
    }

    // Reload all the fields of type volVectorField
    for (uint i = 0; i < volVectorFields_.size(); i++)
    {
        // Load the volume field
        *(volVectorFields_.at(i)) == *(volVectorFieldCopies_.at(i));

        int nOldTimes(volVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            volVectorFields_.at(i)->oldTime() == volVectorFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            volVectorFields_.at(i)->oldTime().oldTime() == volVectorFieldCopies_.at(i)->oldTime().oldTime();
        }

        // TODO. Derek: Should the switch evaluateBoundaries not be implemented here?
        // Find also similar parts below.
        // Evaluate the boundaries
        try{
            DEBUG(adapterInfo("Evaluating the volVector boundary conditions for " + volVectorFields_.at(i)->name()));
            volVectorFields_.at(i)->correctBoundaryConditions();
        } catch (...) {
            DEBUG(adapterInfo("Could not evaluate the boundary for" + volVectorFields_.at(i)->name(), "warning"));
        }
    }

    // Reload all the fields of type surfaceScalarField
    for (uint i = 0; i < surfaceScalarFields_.size(); i++)
    {
        *(surfaceScalarFields_.at(i)) == *(surfaceScalarFieldCopies_.at(i));

        int nOldTimes(surfaceScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            surfaceScalarFields_.at(i)->oldTime() == surfaceScalarFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            surfaceScalarFields_.at(i)->oldTime().oldTime() == surfaceScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
        // no boundary to evaluate
    }

    // Reload all the fields of type surfaceVectorField
    for (uint i = 0; i < surfaceVectorFields_.size(); i++)
    {
        *(surfaceVectorFields_.at(i)) == *(surfaceVectorFieldCopies_.at(i));

        int nOldTimes(surfaceVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            surfaceVectorFields_.at(i)->oldTime() == surfaceVectorFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            surfaceVectorFields_.at(i)->oldTime().oldTime() == surfaceVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
        // no boundary to evaluate
    }

    // NOTE: Add here other field types to read, if needed.

    #ifdef ADAPTER_DEBUG_MODE
        adapterInfo
        (
            "Checkpoint was read. Time = " + std::to_string(runTime_.value())
       );
    #endif

    return;
}


void preciceAdapter::Adapter::writeCheckpoint()
{
    DEBUG(adapterInfo("Writing a checkpoint..."));

    // Store the runTime
    storeCheckpointTime();

    // Store all the fields of type volScalarField
    for (uint i = 0; i < volScalarFields_.size(); i++)
    {
        *(volScalarFieldCopies_.at(i)) == *(volScalarFields_.at(i));
    }

    // Store all the fields of type volVectorField
    for (uint i = 0; i < volVectorFields_.size(); i++)
    {
        *(volVectorFieldCopies_.at(i)) == *(volVectorFields_.at(i));
    }

    // Store all the fields of type volTensorField
    for (uint i = 0; i < volTensorFields_.size(); i++)
    {
        *(volTensorFieldCopies_.at(i)) == *(volTensorFields_.at(i));
    }

    // Store all the fields of type volSymmTensorField
    for (uint i = 0; i < volSymmTensorFields_.size(); i++)
    {
        *(volSymmTensorFieldCopies_.at(i)) == *(volSymmTensorFields_.at(i));
    }

    // Store all the fields of type surfaceScalarField
    for (uint i = 0; i < surfaceScalarFields_.size(); i++)
    {
        *(surfaceScalarFieldCopies_.at(i)) == *(surfaceScalarFields_.at(i));
    }

    // Store all the fields of type surfaceVectorField
    for (uint i = 0; i < surfaceVectorFields_.size(); i++)
    {
        *(surfaceVectorFieldCopies_.at(i)) == *(surfaceVectorFields_.at(i));
    }

    // Store all the fields of type surfaceTensorField
    for (uint i = 0; i < surfaceTensorFields_.size(); i++)
    {
        *(surfaceTensorFieldCopies_.at(i)) == *(surfaceTensorFields_.at(i));
    }

    // NOTE: Add here other types to write, if needed.

    #ifdef ADAPTER_DEBUG_MODE
        adapterInfo
        (
            "Checkpoint for time t = " + std::to_string(runTime_.value()) +
            " was stored."
       );
    #endif

    return;
}

void preciceAdapter::Adapter::readMeshCheckpoint()
{
    DEBUG(adapterInfo("Reading a mesh checkpoint..."));

    //TODO only the meshPhi field is here, which is a surfaceScalarField. The other fields can be removed. 
    // Reload all the fields of type mesh surfaceScalarField
    for (uint i = 0; i < meshSurfaceScalarFields_.size(); i++)
    {    
        // Load the volume field
        *(meshSurfaceScalarFields_.at(i)) == *(meshSurfaceScalarFieldCopies_.at(i));

        int nOldTimes(meshSurfaceScalarFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            meshSurfaceScalarFields_.at(i)->oldTime() == meshSurfaceScalarFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            meshSurfaceScalarFields_.at(i)->oldTime().oldTime() == meshSurfaceScalarFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type mesh surfaceVectorField
    for (uint i = 0; i < meshSurfaceVectorFields_.size(); i++)
    {
        // Load the volume field
        *(meshSurfaceVectorFields_.at(i)) == *(meshSurfaceVectorFieldCopies_.at(i));

        int nOldTimes(meshSurfaceVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            meshSurfaceVectorFields_.at(i)->oldTime() == meshSurfaceVectorFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            meshSurfaceVectorFields_.at(i)->oldTime().oldTime() == meshSurfaceVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }

    // Reload all the fields of type mesh volVectorField
    for (uint i = 0; i < meshVolVectorFields_.size(); i++)
    {
        // Load the volume field
        *(meshVolVectorFields_.at(i)) == *(meshVolVectorFieldCopies_.at(i));

        int nOldTimes(meshVolVectorFields_.at(i)->nOldTimes());
        if (nOldTimes >= 1)
        {
            meshVolVectorFields_.at(i)->oldTime() == meshVolVectorFieldCopies_.at(i)->oldTime();        
        }
        if (nOldTimes == 2)
        {
            meshVolVectorFields_.at(i)->oldTime().oldTime() == meshVolVectorFieldCopies_.at(i)->oldTime().oldTime();
        }
    }
    
    #ifdef ADAPTER_DEBUG_MODE
        adapterInfo
        (
            "Mesh checkpoint was read. Time = " + std::to_string(runTime_.value())
       );
    #endif

    return;
}

void preciceAdapter::Adapter::writeMeshCheckpoint()
{
    DEBUG(adapterInfo("Writing a mesh checkpoint..."));

    // Store all the fields of type mesh surfaceScalar
    for (uint i = 0; i < meshSurfaceScalarFields_.size(); i++)
    {
        *(meshSurfaceScalarFieldCopies_.at(i)) == *(meshSurfaceScalarFields_.at(i));
    }

    // Store all the fields of type mesh surfaceVector
    for (uint i = 0; i < meshSurfaceVectorFields_.size(); i++)
    {
        *(meshSurfaceVectorFieldCopies_.at(i)) == *(meshSurfaceVectorFields_.at(i));
    }

    // Store all the fields of type mesh volVector
    for (uint i = 0; i < meshVolVectorFields_.size(); i++)
    {
        *(meshVolVectorFieldCopies_.at(i)) == *(meshVolVectorFields_.at(i));
    }

    #ifdef ADAPTER_DEBUG_MODE
        adapterInfo
        (
            "Mesh checkpoint for time t = " + std::to_string(runTime_.value()) +
            " was stored."
       );
    #endif

    return;
}

// TODO for the volumes of the mesh, check this part for subcycling. 
void preciceAdapter::Adapter::readVolCheckpoint()
{
    DEBUG(adapterInfo("Reading the mesh volumes checkpoint..."));

    // Reload all the fields of type mesh volVectorField::Internal
    for (uint i = 0; i < volScalarInternalFields_.size(); i++)
    {
        // Load the volume field
        *(volScalarInternalFields_.at(i)) = *(volScalarInternalFieldCopies_.at(i));
        // There are no old times for the internal fields. 
    }

    #ifdef ADAPTER_DEBUG_MODE
        adapterInfo
        (
            "Mesh volumes were read. Time = " + std::to_string(runTime_.value())
       );
    #endif

    return;
}

void preciceAdapter::Adapter::writeVolCheckpoint()
{
    DEBUG(adapterInfo("Writing a mesh volumes checkpoint..."));

    // Store all the fields of type mesh volScalarField::Internal
    for (uint i = 0; i < volScalarInternalFields_.size(); i++)
    {
        *(volScalarInternalFieldCopies_.at(i)) = *(volScalarInternalFields_.at(i));
    }

    #ifdef ADAPTER_DEBUG_MODE
        adapterInfo
        (
            "Mesh volumes checkpoint for time t = " + std::to_string(runTime_.value()) +
            " was stored."
       );
    #endif

    return;
}


void preciceAdapter::Adapter::end()
{
    // Throw a warning if the simulation exited before the coupling was complete
    if (NULL != precice_ && isCouplingOngoing())
    {
        adapterInfo("The solver exited before the coupling was complete.", "warning");
    }

    return;
}

void preciceAdapter::Adapter::teardown()
{
    // If the solver interface was not deleted before, delete it now.
    // Normally it should be deleted when isCouplingOngoing() becomes false.
    if (NULL != precice_)
    {
        DEBUG(adapterInfo("Destroying the preCICE solver interface..."));
        delete precice_;
        precice_ = NULL;
    }

    // Delete the preCICE solver interfaces
    if (interfaces_.size() > 0)
    {
        DEBUG(adapterInfo("Deleting the interfaces..."));
        for (uint i = 0; i < interfaces_.size(); i++)
        {
            delete interfaces_.at(i);
        }
        interfaces_.clear();
    }

    // Delete the copied fields for checkpointing
    if (checkpointing_)
    {
        DEBUG(adapterInfo("Deleting the checkpoints... "));

        // Fields
            // volScalarFields
            for (uint i = 0; i < volScalarFieldCopies_.size(); i++)
            {
                delete volScalarFieldCopies_.at(i);
            }
            volScalarFieldCopies_.clear();
            // volVector
            for (uint i = 0; i < volVectorFieldCopies_.size(); i++)
            {
                delete volVectorFieldCopies_.at(i);
            }
            volVectorFieldCopies_.clear();
            // surfaceScalar
            for (uint i = 0; i < surfaceScalarFieldCopies_.size(); i++)
            {
                delete surfaceScalarFieldCopies_.at(i);
            }
            surfaceScalarFieldCopies_.clear();
            // surfaceVector
            for (uint i = 0; i < surfaceVectorFieldCopies_.size(); i++)
            {
                delete surfaceVectorFieldCopies_.at(i);
            }
            surfaceVectorFieldCopies_.clear();

        //TODO for the internal volume 
            // volScalarInternal
            for (uint i = 0; i < volScalarInternalFieldCopies_.size(); i++)
            {
                delete volScalarInternalFieldCopies_.at(i);
            }
            volScalarInternalFieldCopies_.clear();

        // NOTE: Add here delete for other types, if needed

        checkpointing_ = false;
    }

    // Delete the FF module
    if(NULL != FF_)
    {
        DEBUG(adapterInfo("Destroying the FF module..."));
        delete FF_;
        FF_ = NULL;
    }

    // NOTE: Delete your new module here

    return;
}

preciceAdapter::Adapter::~Adapter()
{
    teardown();

    return;
}
