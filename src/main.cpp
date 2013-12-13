#include <mpi.h>
#include "RootEngine.h"

enum {rootTASK = 0};

static void master(std::string cfgname);
static void slave(std::string cfgname, int id);

int main(int argc, char** argv)
{
    int taskid;

    /*-------------mpi section-------------*/
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    if (taskid == rootTASK)
    {
        master("default.cfg");
    }
    else
    {
        slave("default.cfg", taskid);
    }

    MPI_Finalize();
	/*-------------end of mpi section-------------*/
    return 0;
}

static void master(std::string filename)
{
    int ntasks;
    RootEngine *  engine;
    //MCCalculator::Data * curData, * cumData;
    int buffsize;
    unsigned long nbSteps;
    double * recvbuff, *sendbuff;
    ProgramSettings * settings;

    settings = new ProgramSettings();
    try
    {
    	settings->read(filename);
    	settings->print();
    }catch(const ProgramSettings::Exception & ex)
    {
    	std::cout << "Exception caught:\t" << ex.what() << std::endl;
    	delete settings;
    	MPI_Abort(MPI_COMM_WORLD, -1);
    }
    catch(...)
    {
    	std::cout << "Exception caught:\t unknown" << std::endl;
    	delete settings;
    	MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /* Find out how many processes there are in the default
     communicator */
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    /*allocations and initializations*/
    engine = new RootEngine(settings, 0, ntasks);
    buffsize = engine->getDataSize();
    recvbuff = new double[buffsize];
    sendbuff = new double[buffsize];

    /*wait while all preparations are finished*/
    MPI_Barrier(MPI_COMM_WORLD);

    /* Loop over getting new work requests until there is no more work
     to be done */
    while (!engine->done())
    {
    	/*how many MC steps to perform this time?*/
    	nbSteps = engine->getNbSteps();
        /* Send the slave a new work unit */
        MPI_Bcast(&nbSteps, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        /*master also works*/
        engine->run();
        engine->transferDataTo(sendbuff);

        /* Receive results from slaves */
        MPI_Reduce(sendbuff, recvbuff, buffsize, MPI_DOUBLE, MPI_SUM,
					rootTASK, MPI_COMM_WORLD );

        /**append cummulative data and save result**/
        engine->appendData(recvbuff);
        engine->saveData();
    }
    /* Tell all the slaves to exit by sending 0 nb of steps to perform*/
    nbSteps = engine->getNbSteps();
    /* Send the slave a new work unit */
    MPI_Bcast(&nbSteps, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    delete[] recvbuff;
    delete[] sendbuff;
    delete engine;
    delete settings;

    std::cout << "Done." << std::endl;
}

static void slave(std::string filename, int taskid)
{
    Engine * engine;
    int buffsize;
    double * recvbuff, *sendbuff;
    unsigned long nbSteps;
    ProgramSettings * settings;

    settings = new ProgramSettings();
    try
    {
    	settings->read(filename);
    }catch(const ProgramSettings::Exception & ex)
    {
    	std::cout << "Exception caught:\t" << ex.what() << std::endl;
    	delete settings;
    	MPI_Abort(MPI_COMM_WORLD, -1);
    }
    catch(...)
    {
    	std::cout << "Exception caught:\t unknown" << std::endl;
    	delete settings;
    	MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /*allocations and initializations*/
    engine = new Engine(settings, taskid);
    buffsize = engine->getDataSize();
    recvbuff = new double[buffsize];
    sendbuff = new double[buffsize];

    MPI_Barrier(MPI_COMM_WORLD);
    while (true)
    {
        /* Receive a message from the master */
        MPI_Bcast(&nbSteps, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

        /* if nb steps received is equal to 0, stop working */
        if (nbSteps == 0)
            break;
        /*else reset nb of steps*/
        engine->setNbSteps(nbSteps);

        /* Do the work */
        engine->run();
        engine->transferDataTo(sendbuff);

        /* Send the result back */
        MPI_Reduce(sendbuff, recvbuff, buffsize, MPI_DOUBLE, MPI_SUM,
					rootTASK, MPI_COMM_WORLD );
    }
    delete[] recvbuff;
    delete[] sendbuff;
    delete engine;
    delete settings;
}
