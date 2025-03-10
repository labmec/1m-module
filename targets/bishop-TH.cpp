#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzskylstrmatrix.h>
#include <TPZMultiphysicsCompMesh.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include "pzfstrmatrix.h"
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include "Elasticity/TPZHybridMixedElasticityUP.h"
#include "TPZMatInterfaceHybridElasticityStokes.h"
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "meshpath_config.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "ProblemData.h"
#include "pzintel.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include <pzelementgroup.h>
#include <pzcondensedcompel.h>
#include "pzskylmat.h"
#include "TPZMatrixSolver.h"
#include "TPZPardisoSolver.h"
#include "TPZSYSMPPardiso.h"
#include "TPZSparseMatRed.h"
#include "Elasticity/TPZElasticityTH.h"

const int global_nthread = 64;
const int global_pord_bc = 10;

using namespace std;

// functions declaration
TPZGeoMesh *ReadMeshFromGmsh(std::string file_name, ProblemData *problem_data);
TPZCompMesh *CreateCMeshU(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshP(ProblemData *simData, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh, TPZAnalyticSolution* sol);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data);

#ifdef PZ_LOG
static TPZLogger logger("pz.1mmodule");
#endif

int main(int argc, char *argv[])
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG(std::string(INPUTDIR) + "/log4cxx.cfg");
#endif

    std::cout << "--------- Starting simulation ---------" << std::endl;

    // Reading problem data from json
    std::string jsonfilename = "conv-bishop-";
    int meshref = 1;
    if(argc > 1) meshref = atoi(argv[1]);
    jsonfilename += to_string(meshref) + "-tet.json";
    // jsonfilename = "bishop-beam-UP.json";
    
    ProblemData problemdata;
    std::cout << "json input filename: " << jsonfilename << std::endl;
    problemdata.ReadJson(std::string(MESHES_DIR) + "/" + jsonfilename);
    
    // Create gmesh
    const int pord = problemdata.DisppOrder();
    TPZGeoMesh *gmesh = nullptr;
    std::string filename = problemdata.MeshName();
    gmesh = ReadMeshFromGmsh(std::string(MESHES_DIR) + "/" + filename, &problemdata);
    if (0)
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
    }

    const REAL young = problemdata.DomainVec()[0].E;
    const REAL poisson = argc > 2? atof(argv[2]) : problemdata.DomainVec()[0].nu;
    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    elas->fE = young;
    elas->fPoisson = poisson;
    elas->fProblemType = TElasticity3DAnalytic::ETestShearMoment;

    TPZCompMesh *cmesh_u = CreateCMeshU(&problemdata, gmesh);
    if (0)
    {
        std::ofstream out("cmesh_u.txt");
        cmesh_u->Print(out);
    }
    TPZCompMesh *cmesh_p = CreateCMeshP(&problemdata, gmesh);
    if (0)
    {
        std::ofstream out("cmesh_p.txt");
        cmesh_u->Print(out);
    }
    problemdata.MeshVector().resize(2);

    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiphysicsMesh(&problemdata, gmesh, elas);
    if (0)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
    }
    int64_t ndofs_total = cmesh_m->NEquations();
    std::cout << "ndofs_total = " << ndofs_total << std::endl;
    int64_t ndofs_condensed = cmesh_m->NEquations();
    std::cout << "ndofs_condensed = " << ndofs_condensed << std::endl;

    std::cout << "Number of equations: " << cmesh_m->NEquations() << std::endl;

    // Analysis
    // Solve Multiphysics
    RenumType renum = RenumType::EMetis;
    
    TPZSimpleTimer time_band;        
    TPZLinearAnalysis an(cmesh_m, renum);
    std::cout << "Total time optimize band = " << time_band.ReturnTimeDouble() / 1000. << " s" << std::endl;
    SolveProblemDirect(an, cmesh_m);
    if (0){
        std::ofstream out("cmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out);
        std::ofstream out2("cmesh.txt");
        cmesh_m->Print(out2);
    }

    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh_m, &problemdata);

    // Calculating error
    if (elas->fProblemType != 0)
    {
        an.SetExact(elas->ExactSolution());
        an.SetThreadsForError(global_nthread);
        std::ofstream out("TH-bishop-convergence.txt",std::ios::app);
        out << "\n----------------- Starting new simulation -----------------" << std::endl;
        std::cout << "\n----------------- Starting error computation -----------------" << std::endl;
        out << "meshref: " << meshref << ", poisson: " << poisson << std::endl;
        out << "ndofs: " << ndofs_total << ", after condensation: " << ndofs_condensed << std::endl;

        TPZMaterial *mat = cmesh_m->FindMaterial(problemdata.DomainVec()[0].matID);
        TPZMatErrorCombinedSpaces<STATE> *materr = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
        TPZManVector<REAL, 10> Errors(materr->NEvalErrors());
        
        bool store_errors = false;
        std::ofstream ErroOut("myerrors.txt", std::ios::app);
        std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
        an.PostProcessError(Errors, store_errors, ErroOut);
        std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
        std::cout << "Time PostProc Error = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count()/1000. << " s" << std::endl;

        std::cout << "Computed errors." << std::endl;
        std::cout << "L2 p error | L2 p_ex | L2 u error | L2 u_ex | L2 divu error | L2 divu_ex | L2 sigma error | L2 sigma_Ex" << std::endl;
        std::cout << "Errors = ";
        std::cout << Errors << std::endl;
        out << "Errors = ";
        out << Errors << std::endl;
        out << "ErrorsFixedPrecision = ";
        out.precision(15);
        out.setf(std::ios::fixed);
        out << Errors << std::endl;
    }

    // deleting stuff
    if (cmesh_m)
        delete cmesh_m;
    if (cmesh_u)
        delete cmesh_u;
    if (cmesh_p)
        delete cmesh_p;
    if (gmesh)
        delete gmesh;

    std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh *ReadMeshFromGmsh(std::string file_name, ProblemData *problem_data)
{
    // read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    TPZGmshReader reader;
    reader.GeometricGmshMesh(file_name, gmesh);
    gmesh->BuildConnectivity();
    return gmesh;
}

TPZCompMesh *CreateCMeshU(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_u = new TPZCompMesh(gmesh);
    cmesh_u->SetName("CMesh_U");

    std::set<int> materialIDs;

    // domain's material (2D or 3D)
    if (simData->DomainVec().size() != 0)
    {
        cmesh_u->SetDefaultOrder(simData->DisppOrder());

        cmesh_u->SetDimModel(simData->Dim());

        cmesh_u->SetAllCreateFunctionsContinuous();        

        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        mat->SetNStateVariables(simData->Dim());
        cmesh_u->InsertMaterialObject(mat);

        materialIDs.insert(simData->DomainVec()[0].matID);

        // boundary conditions' material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;
            int type = 0;
            if (bc.matID == 2) type = 0; //imposed displacement at all directions
            else if (bc.matID == 3 || bc.matID == 4) type = 1; //surface traction
            auto BCmat = mat->CreateBC(mat, bc.matID, type, val1, val2);
            cmesh_u->InsertMaterialObject(BCmat);
            materialIDs.insert(bc.matID);
        }

        cmesh_u->AutoBuild(materialIDs);
        gmesh->ResetReference();
    }

    // expanding the solution vector
    cmesh_u->ExpandSolution();

    simData->MeshVector()[0] = cmesh_u;

    return cmesh_u;
}

TPZCompMesh *CreateCMeshP(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("CMesh_P");

    std::set<int> materialIDs;

    if (simData->DomainVec().size() != 0) // if true, there is a 2D/3D domain
    {
        cmesh_p->SetDimModel(simData->Dim());

        cmesh_p->SetDefaultOrder(simData->DisppOrder() - 1);
        cmesh_p->SetAllCreateFunctionsContinuous();
        
        // domain's material
        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(mat);

        materialIDs.insert(simData->DomainVec()[0].matID);

        cmesh_p->AutoBuild(materialIDs);

        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        gmesh->ResetReference();

        materialIDs.clear();
    }

    cmesh_p->ExpandSolution();

    simData->MeshVector()[1] = cmesh_p;

    return cmesh_p;
}

TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh, TPZAnalyticSolution* sol)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    cmesh_m->SetName("CMesh_M");

    cmesh_m->SetDefaultOrder(simData->DisppOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // Creating Materials

    // 1. For domain
    if (simData->DomainVec().size() != 0)
    {
        REAL young = simData->DomainVec()[0].E;
        REAL poisson = simData->DomainVec()[0].nu;
        TPZElasticityTH::AnalysisType analysistype = TPZElasticityTH::AnalysisType::EGeneral;

        if (dynamic_cast<TElasticity3DAnalytic*>(sol))
        {
            TElasticity3DAnalytic* elas = dynamic_cast<TElasticity3DAnalytic*>(sol);
            young = elas->fE;
            poisson = elas->fPoisson;
            analysistype = TPZElasticityTH::AnalysisType::EGeneral;
            if (elas->fProblemType == TElasticity3DAnalytic::ENone) sol = nullptr;
        }

        TPZElasticityTH *mat = new TPZElasticityTH(simData->DomainVec()[0].matID, simData->Dim(), young, poisson, analysistype);
        if (sol) mat->SetExactSol(sol->ExactSolution(), 4);
        cmesh_m->InsertMaterialObject(mat);

        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;
            int type = 0;
            if (bc.matID == 2) type = 0; //imposed displacement at all directions
            else if (bc.matID == 3 || bc.matID == 4) type = 1; //surface traction

            val2 = bc.value;

            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            if (sol) matBC2->SetForcingFunctionBC(sol->ExactSolution(), global_pord_bc);
            cmesh_m->InsertMaterialObject(matBC);
        }
    }

    TPZManVector<int, 2> active_approx_spaces(simData->MeshVector().size(), 1);

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, simData->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();

    return cmesh_m;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    // TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);

    /// Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); // ELU //ECholesky // ELDLt
    an.SetSolver(step);

    // assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble() / 1000. << " s" << std::endl;

    /// solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble() / 1000. << " s" << std::endl;

    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data)
{

    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";

    TPZVec<std::string> fields = {
        "Displacement",
        "Pressure",
        "Stress",
        "Strain",
        "VonMises"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, problem_data->Resolution(), problem_data->Dim());
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;

    return;
}