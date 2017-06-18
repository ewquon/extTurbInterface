/*
 * Interface/"glue code" for the Vestas Turbine Simulator
 * - Created by Eliot Quon 2017-06-16
 *   based on glue code created by Ganesh Vijayakumar for SOWFA/OpenFAST coupling
 */

#ifndef VTS_h
#define VTS_h

#include <iostream>
#include <vector>
#include <fstream>
#include <set>
#include <map>

#include "sys/stat.h"
#include "H5Cpp.h"

//#include "math.h"
//#include <string>
//#include <cstring>
//#include <stdexcept>
//#include "dlfcn.h"
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"

// Definitions of VTS functions (from shared library, written in Pascal)
#include "vts_interface.H"

namespace externalTurbineModel
{

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

struct globalTurbineDataType
{
    std::string VTSMasterFile;
    std::string VTSLoadCase;

    int TurbID;
    std::vector<double> TurbineBasePos;
    std::vector<double> TurbineHubPos;
    int numForcePtsBlade;
    int numForcePtsTower;
};

enum simStartType {
  init = 0,
  trueRestart = 1,
  restartDriverInit = 2,
  simStartType_END
};


class vtsInputs
{

public: 

    // MPI
    MPI_Comm comm;
    int nTurbinesGlob;

    // Development
    bool dryRun; // If this is true, class will simply go through allocation and deallocation of turbine data
    bool debug;  // Write out extra information if this flags is turned on

    // Simulation paramters
    double tStart;
    simStartType simStart;
//  int nEveryCheckPoint; // check point files will be written every 'nEveryCheckPoint' time steps  
    double tMax;
    double dtVTS;  

    std::vector<globalTurbineDataType> globTurbineData;

    // Constructor 
    vtsInputs();
 
    // Destructor
    ~vtsInputs()
    {}

};


class VTS
{

private:

    bool dryRun;
    bool debug;

    bool timeZero;
    simStartType simStart;
    double tStart;
//  int nEveryCheckPoint;
    double tMax;
    double dtVTS;

    std::vector<globalTurbineDataType> globTurbineData;

    // array-level quantities
    int nTurbinesProc;
    int nTurbinesGlob;

    std::vector<int> TurbID;
    std::vector<std::vector<float> > TurbineBasePos;
//  std::vector<std::vector<float> > TurbineHubPos;

    std::vector<std::string > VTSMasterFileName;
    std::vector<std::string > VTSLoadCaseName;

    int nt_global; // Number of time steps performed
    int ntStart; // The time step to start the VTS simulation

    // turbine discretization
    std::vector<int> numBlades;           // Number of blades
    std::vector<int> numForcePtsBlade;
    std::vector<int> numForcePtsTower;
    std::vector<int> numVelPtsBlade;
    std::vector<int> numVelPtsTower;

    // Velocity at force nodes
    // - Store temporarily to interpolate to the velocity nodes
    std::vector<std::vector<std::vector<double> > > forceNodeVel;

    // Position and velocity data at the velocity (e.g., Aerodyn) nodes
    // - (nTurbines, nTimesteps * nPoints * 6)
    std::vector<std::vector<double> > velNodeData;

    hid_t velNodeDataFile; // HDF-5 tag of file containing velocity node data

    // Turbine Number is DIFFERENT from TurbID; turbine Number simply runs from
    // 0:n-1 locally and globally.
    // - Mapping global turbine number to processor number
    std::map<int, int> turbineMapGlobToProc;
    // - Mapping local to global turbine number
    std::map<int, int> turbineMapProcToGlob;
    // - Reverse Mapping global turbine number to local turbine number
    std::map<int, int> reverseTurbineMapProcToGlob;
    // - Set of processors containing at least one turbine 
    std::set<int> turbineSetProcs;
    // - Same as the turbineSetProcs, but as an integer array
    std::vector<int> turbineProcs;

    // MPI variables
    MPI_Comm mpiComm;

    int vtsMPIGroupSize;
    MPI_Group vtsMPIGroup;
    MPI_Comm  vtsMPIComm;
    int vtsMPIRank;

    MPI_Group worldMPIGroup;
    int worldMPIRank;

    // VTS return status
    int ErrStat;

public: 

    // Constructor 
    VTS();
  
    // Destructor
    ~VTS()
    {}

    void setInputs(const vtsInputs &);  

    // solution control
    void init();
    void solution0();
    void step(bool writeRestart = true);
    void end();

    // restart I/O
    hid_t openVelocityDataFile(bool createFile);
    void readVelocityData(int nTimesteps);
//  void writeVelocityData(hid_t h5file, int iTurb, int iTimestep, OpFM_InputType_t iData, OpFM_OutputType_t oData);
//  herr_t closeVelocityDataFile(int nt_global, hid_t velDataFile);

    // get/setters
//  void setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
    void allocateTurbinesToProcsSimple();
    void getApproxHubPos(std::vector<double> & currentCoords, int iTurbGlob);
    void getHubPos(std::vector<double> & currentCoords, int iTurbGlob);
    void getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob);

    void getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob);
    void setVelocity(std::vector<double> & velocity, int iNode, int iTurbGlob);
    void setVelocityForceNode(std::vector<double> & velocity, int iNode, int iTurbGlob);
    void interpolateVel_ForceToVelNodes();

    void getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob);
    void getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob);
    void getForce(std::vector<double> & force, int iNode, int iTurbGlob);
    double getChord(int iNode, int iTurbGlob);

    // access functions
    int get_ntStart() { return ntStart; }
    bool isDryRun() { return dryRun; }
    bool isDebug() { return debug; }
    simStartType get_simStartType() { return simStart; }
    bool isTimeZero() { return timeZero; }

//  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; } // Get processor number of a turbine with global id 'iTurbGlob'
//  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
//  int get_nTurbinesGlob() { return nTurbinesGlob; } 
//  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
//  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
//  int get_numVelPtsTower(int iTurbGlob) { return get_numVelPtsTowerLoc(get_localTurbNo(iTurbGlob)); }
//  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
//  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
//  int get_numForcePtsTower(int iTurbGlob) { return get_numForcePtsTowerLoc(get_localTurbNo(iTurbGlob)); }
//  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

    // processing functions
    void computeTorqueThrust
    (
        int iTurGlob,
        std::vector<double>& torque,
        std::vector<double>& thrust
    );

private:

//  void checkError(const int ErrStat, const char * ErrMsg);
    inline bool checkFileExists(const std::string& name);

    void allocateMemory();
  
    int get_numBladesLoc(int iTurbLoc) { return numBlades[iTurbLoc]; }
    int get_numVelPtsBladeLoc(int iTurbLoc) { return numVelPtsBlade[iTurbLoc]; }
    int get_numVelPtsTowerLoc(int iTurbLoc) { return numVelPtsTower[iTurbLoc]; }
    int get_numVelPtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numVelPtsBlade[iTurbLoc] + numVelPtsTower[iTurbLoc]; }
    int get_numForcePtsBladeLoc(int iTurbLoc) { return numForcePtsBlade[iTurbLoc]; }
    int get_numForcePtsTowerLoc(int iTurbLoc) { return numForcePtsTower[iTurbLoc]; }
    int get_numForcePtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numForcePtsBlade[iTurbLoc] + numForcePtsTower[iTurbLoc]; }
    
//  void applyVelocityData(int iPrestart, int iTurb, OpFM_OutputType_t cDriver_Output_to_FAST, std::vector<double> & velData) ;
};


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

} // namespace vts

#endif
