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

const int global_nthread = 64;
const int global_pord_bc = 4;

using namespace std;

enum HdivType
{
    EStandard,
    EConstant
};

enum EMatid
{
    ENone,
    EDomain,
    EZeroNormalDisp,
    EZeroNormalStress,
    EZeroTangentialDisp,
    EZeroTangentialStress,
    EPressure,
    EStiffner
};

// functions declaration
TPZGeoMesh *ReadMeshFromGmsh(std::string file_name, ProblemData *problem_data);
void InsertLagrandeMultipliers(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshU(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshP(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshG(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshPm(ProblemData *simData, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh, TElasticity3DAnalytic* elas);
void CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh *cmesh_m);
void InsertBCInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh);
void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data);
void SolveProblemSparseMatRed(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData* problem_data);

#ifdef PZ_LOG
static TPZLogger logger("pz.1mmodule");
#endif

int main(int argc, char *argv[])
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG(std::string(MESHES_DIR) + "/" + "log4cxx.cfg");
#endif

    std::cout << "--------- Starting simulation ---------" << std::endl;

    // Reading problem data from json
    std::string jsonfilename = "test-fully-hybrid.json";
    int meshref = 1;
    
    ProblemData problemdata;
    std::cout << "json input filename: " << jsonfilename << std::endl;
    problemdata.ReadJson(std::string(MESHES_DIR) + "/" + jsonfilename);
    
    // Create gmesh
    const int pord = problemdata.DisppOrder();
    TPZGeoMesh *gmesh = nullptr;
    std::string filename = problemdata.MeshName();
    gmesh = ReadMeshFromGmsh(std::string(MESHES_DIR) + "/" + filename, &problemdata);
    
    InsertLagrandeMultipliers(&problemdata, gmesh);
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

    // Create compmeshes
    if (problemdata.DomainVec().size() > 1)
        DebugStop(); // Please implement the next lines correctly if many domains

    TPZCompMesh *cmesh_u = CreateCMeshU(&problemdata, gmesh);
    {
        std::ofstream out("cmesh_u.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_u, out);
        std::ofstream out2("cmesh_u.txt");
        cmesh_u->Print(out2);
    }
    TPZCompMesh *cmesh_p = CreateCMeshP(&problemdata, gmesh);
    {
        std::ofstream out("cmesh_p.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_p, out);
        std::ofstream out2("cmesh_p.txt");
        cmesh_p->Print(out2);
    }
    TPZCompMesh *cmesh_pm = nullptr;
    TPZCompMesh *cmesh_g = nullptr;
    if (problemdata.CondensedElements() and fabs(poisson - 0.5) < 1.e-3) {
        cmesh_g = CreateCMeshG(&problemdata, gmesh);
        cmesh_pm = CreateCMeshPm(&problemdata, gmesh);
    }
    else {
        problemdata.MeshVector().resize(2);
    }
    
    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiphysicsMesh(&problemdata, gmesh, elas);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
    }
    {
        std::ofstream out("cmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh_m, out);
        std::ofstream out2("cmesh.txt");
        cmesh_m->Print(out2);
    }

    if (problemdata.CondensedElements())
        CondenseElements(&problemdata, cmesh_m);

    // Analysis
    // Solve Multiphysics
    TPZLinearAnalysis an(cmesh_m, RenumType::EMetis);
    SolveProblemDirect(an, cmesh_m, &problemdata);
    // SolveProblemSparseMatRed(an, cmesh_m, &problemdata);

    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh_m, &problemdata);
    
    // Calculating error
    an.SetExact(elas->ExactSolution());
    an.SetThreadsForError(global_nthread);
    std::ofstream out("bishop-convergence.txt",std::ios::app);
    out << "\n----------------- Starting new simulation -----------------" << std::endl;
    std::cout << "\n----------------- Starting error computation -----------------" << std::endl;
    out << "meshref: " << meshref << ", poisson: " << poisson << std::endl;
    
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
    // error_sigma - error_energy - error_div_sigma - error_u - error_r - error_as - energy_norm_exact_sol
    std::cout << "L2 p error | L2 p_ex | L2 u error | L2 u_ex | L2 divu error | L2 divu_ex | L2 sigma error | L2 sigma_Ex" << std::endl;
//    std::cout.precision(15);
//    std::cout.setf(std::ios::fixed);
    std::cout << "Errors = ";
    std::cout << Errors << std::endl;
    out << "Errors = ";
    out << Errors << std::endl;
    out << "ErrorsFixedPrecision = ";
    out.precision(15);
    out.setf(std::ios::fixed);
    out << Errors << std::endl;

    // deleting stuff
    if (cmesh_m)
        delete cmesh_m;
    if (cmesh_u)
        delete cmesh_u;
    if (cmesh_p)
        delete cmesh_p;
    if (cmesh_g)
        delete cmesh_g;
    if (cmesh_pm)
        delete cmesh_pm;
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

void InsertLagrandeMultipliers(ProblemData *simData, TPZGeoMesh *gmesh)
{
    int64_t nEl = gmesh->NElements();

    // We look for two domain neighbour elements 
    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);

        if (!geoEl)
            continue;
        if (geoEl->HasSubElement())
            continue;
        if (geoEl->Dimension() != gmesh->Dimension())
            continue;

        int nside = geoEl->NSides();

        if (simData->DomainVec().size() != 0)
        {
            for (int side = 0; side < nside; side++)
            {
                if (geoEl->SideDimension(side) != gmesh->Dimension() - 1)
                    continue;

                TPZGeoElSide geoElSide(geoEl, side);
                TPZGeoElSide neighbour = geoElSide.Neighbour();

                if (neighbour == geoElSide)
                    continue;
                if (neighbour.Element()->HasSubElement())
                    continue;

                TPZGeoElSide neighbourOfneighbour = neighbour; //This was introduced to keep track of the neighbour element because in the fully-hybrid formulation we need to add a pair of tangential stress
                while (neighbourOfneighbour != geoElSide)
                {

                    if (neighbourOfneighbour.Element()->Dimension() == gmesh->Dimension() - 1) //in this case, we already inserted the lagrange multipliers 
                    {
                        break;
                    }

                    neighbourOfneighbour = neighbourOfneighbour.Neighbour();
                }

                if (neighbourOfneighbour == geoElSide)
                {
                    TPZGeoElBC(geoElSide, 10); //Adding a tangential stress Lagrange multiplier as a neighbour to the element side
                    TPZGeoElBC(neighbour, 10); //Adding a tangential stress Lagrange multiplier as a neighbour to the neighbour element side
                    neighbourOfneighbour = neighbour.Neighbour();

                    neighbour = geoElSide.Neighbour(); //this will return the tangential stress Lagrange multiplier
                    TPZGeoElBC(neighbour, 15); //Adding a tangential displacement Lagrange multiplier as a neighbour to the tangential stress element
                    neighbour = neighbour.Neighbour();
                }
            }
        }
    }
    gmesh->BuildConnectivity();
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

        if (simData->HdivType() == HdivType::EConstant)
        {
            cmesh_u->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
        }
        else if (simData->HdivType() == HdivType::EStandard)
        {
            cmesh_u->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
        }

        cmesh_u->SetAllCreateFunctionsHDiv();

        auto *mat_normal = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_u->InsertMaterialObject(mat_normal);

        materialIDs.insert(simData->DomainVec()[0].matID);

        // boundary conditions' material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;

            auto BCmat = mat_normal->CreateBC(mat_normal, bc.matID, bc.type, val1, val2);
            cmesh_u->InsertMaterialObject(BCmat);
            materialIDs.insert(bc.matID);
        }

        cmesh_u->AutoBuild(materialIDs);

        // Increasing internal function order
        int64_t ncEl = cmesh_u->NElements();
        for (int64_t cEl = 0; cEl < ncEl; cEl++)
        {
            TPZCompEl *compEl = cmesh_u->Element(cEl);

            // only in those elements whose dimension equals to the simulation dim
            if (compEl->Dimension() == simData->Dim())
            {
                // dynamic casting the compEl object to use the ForceSideOrder function
                TPZInterpolatedElement *intercEl = dynamic_cast<TPZInterpolatedElement *>(compEl);

                // checking if the dynamic cast exists
                if (!intercEl)
                    continue;

                intercEl->ForceSideOrder(compEl->Reference()->NSides() - 1, simData->DisppOrder() + 1);
            }
        }

        int64_t ncon = cmesh_u->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_u->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(2);
        }
        for (auto& cel : cmesh_u->ElementVec())
        {
            int mat = cel->Reference()->MaterialId();
            if (mat == 2) //big number
            {
                ncon = cel->NConnects();
                for (int i = 0; i < ncon; i++)
                {
                    TPZConnect& con = cel->Connect(i);
                    con.SetLagrangeMultiplier(0);
                }
            }
        }
        gmesh->ResetReference();

        // tangent displacement material
        auto mat_tan = new TPZNullMaterial<>(15);
        mat_tan->SetNStateVariables(simData->Dim() - 1); // In 3D, there are 2 state variables (one at each tangential direction)
        cmesh_u->InsertMaterialObject(mat_tan);

        materialIDs.insert(15);

        cmesh_u->ApproxSpace().CreateDisconnectedElements(true);
        cmesh_u->SetDefaultOrder(simData->LambdapOrder());
        cmesh_u->SetDimModel(simData->Dim() - 1);
        cmesh_u->AutoBuild(materialIDs);

        ncon = cmesh_u->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_u->ConnectVec()[i];
            if (newnod.LagrangeMultiplier() == 0)
                newnod.SetLagrangeMultiplier(1);
        }

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

        if (simData->HdivType() == HdivType::EConstant)
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->SetDefaultOrder(0);
        }
        else if (simData->HdivType() == HdivType::EStandard)
        {
            cmesh_p->SetDefaultOrder(simData->DisppOrder() + 1);
            cmesh_p->SetAllCreateFunctionsContinuous();
        }

        cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

        // domain's material
        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(mat);

        materialIDs.insert(simData->DomainVec()[0].matID);

        cmesh_p->AutoBuild(materialIDs);
        gmesh->ResetReference();

        materialIDs.clear();

        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(3);
        }

        // matlambda traction material
        auto matLambda = new TPZNullMaterial<>(simData->LambdaID());
        matLambda->SetNStateVariables(simData->Dim() - 1); // In 3D, lambda has 2 state variables (one at each tangential direction)
        cmesh_p->InsertMaterialObject(matLambda);

        materialIDs.insert(simData->LambdaID());

        // traction on boundary material
        for (const auto &bc : simData->TangentialBCs())
        {
            auto matLambdaBC = new TPZNullMaterial<>(bc.matID);
            matLambdaBC->SetNStateVariables(simData->Dim() - 1);
            cmesh_p->InsertMaterialObject(matLambdaBC);

            materialIDs.insert(bc.matID);
        }

        if (simData->LambdapOrder() > 0)
        {
            cmesh_p->SetAllCreateFunctionsContinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }

        cmesh_p->SetDefaultOrder(simData->LambdapOrder());
        cmesh_p->SetDimModel(simData->Dim() - 1);
        cmesh_p->AutoBuild(materialIDs);

        gmesh->ResetReference();

        ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            if (newnod.LagrangeMultiplier() == 0)
                newnod.SetLagrangeMultiplier(1);
        }
    }

    cmesh_p->ExpandSolution();

    simData->MeshVector()[1] = cmesh_p;

    return cmesh_p;
}

TPZCompMesh *CreateCMeshG(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_g = new TPZCompMesh(gmesh);

    cmesh_g->SetName("CMesh_g");
    cmesh_g->SetDefaultOrder(0);
    cmesh_g->SetDimModel(simData->Dim());

    cmesh_g->SetAllCreateFunctionsDiscontinuous();

    auto mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_g->InsertMaterialObject(mat);

    int64_t ncel = cmesh_g->NElements();

    cmesh_g->AutoBuild();
    cmesh_g->AdjustBoundaryElements();
    cmesh_g->CleanUpUnconnectedNodes();

    simData->MeshVector()[2] = cmesh_g;

    return cmesh_g;
}

TPZCompMesh *CreateCMeshPm(ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_pm = new TPZCompMesh(gmesh);

    cmesh_pm->SetName("CMesh_pm");
    cmesh_pm->SetDefaultOrder(0);
    cmesh_pm->SetDimModel(simData->Dim());

    cmesh_pm->SetAllCreateFunctionsDiscontinuous();

    auto mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
    cmesh_pm->InsertMaterialObject(mat);

    cmesh_pm->AutoBuild();
    cmesh_pm->AdjustBoundaryElements();
    cmesh_pm->CleanUpUnconnectedNodes();

    int64_t ncon = cmesh_pm->NConnects();
    for (int64_t i = 0; i < ncon; i++)
    {
        TPZConnect &newnod = cmesh_pm->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(4);
    }

    simData->MeshVector()[3] = cmesh_pm;

    return cmesh_pm;
}

TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh, TElasticity3DAnalytic* elas)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    cmesh_m->SetName("CMesh_M");

    cmesh_m->SetDefaultOrder(simData->DisppOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // Creating Materials

    // 1. For domain
    if (simData->DomainVec().size() != 0)
    {
        const REAL young = elas->fE;
        const REAL poisson = elas->fPoisson;

        TPZHybridMixedElasticityUP *mat = new TPZHybridMixedElasticityUP(simData->DomainVec()[0].matID, simData->Dim(), young, poisson, AnalysisType::EGeneral);
        mat->SetExactSol(elas->ExactSolution(), 4);
        cmesh_m->InsertMaterialObject(mat);

        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            matBC2->SetForcingFunctionBC(elas->ExactSolution(), global_pord_bc);
            cmesh_m->InsertMaterialObject(matBC);
        }
        
        for (const auto &bc : simData->TangentialBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            matBC2->SetForcingFunctionBC(elas->ExactSolution(), global_pord_bc);
            cmesh_m->InsertMaterialObject(matBC);
        }

        // 3 - Material for tangential displacement
        TPZNullMaterialCS<> *mat_tan = new TPZNullMaterialCS<>(15);
        double dim = simData->Dim() - 1;
        mat_tan->SetDimension(dim);
        mat_tan->SetNStateVariables(dim);
        cmesh_m->InsertMaterialObject(mat_tan);

        // 3 - Material for tangential traction
        TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->LambdaID());
        matLambda->SetDimension(dim);
        matLambda->SetNStateVariables(dim);
        cmesh_m->InsertMaterialObject(matLambda);

        // 4 - Material for interfaces (Inner)
        TPZMatInterfaceHybridElasticityStokes *matInterfaceLeft = new TPZMatInterfaceHybridElasticityStokes(simData->InterfaceID(), simData->Dim() - 1);
        matInterfaceLeft->SetMultiplier(1.);
        cmesh_m->InsertMaterialObject(matInterfaceLeft);
    }

    TPZManVector<int, 2> active_approx_spaces(simData->MeshVector().size(), 1);

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, simData->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();

    InsertBCInterfaces(cmesh_m, simData, gmesh);
    InsertInterfaces(cmesh_m, simData, gmesh);

    if (simData->CondensedElements())
    {
        cmesh_m->SetName("CMesh_M_BeforeCond");
        cmesh_m->ComputeNodElCon();
    }

    return cmesh_m;
}

void CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh *cmesh_m)
{
    int64_t ncompEl = cmesh_m->ElementVec().NElements();
    int dim = cmesh_m->Reference()->Dimension();

    std::set<int64_t> externalNode;
    std::vector<int64_t> groupIndex;
    groupIndex.reserve(ncompEl);
    TPZStack<TPZElementGroup *> elGroups;
    int count = 0;

    // Creating the element groups for the domain
    for (int64_t el = 0; el < ncompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);

        if (compEl->Dimension() != dim)
            continue;

        int numConnectExt = (simData->MeshVector().size() == 2)? 0 : 1; //If no meshes were created for distributed flux and mean pressure, we condense all pressures. Otherwise, we do not condense the distributed flux
        int nConnect = compEl->NConnects();

        for (int ic = nConnect - numConnectExt; ic < nConnect; ic++)
        {
            int64_t conIndex = compEl->ConnectIndex(ic);
            externalNode.insert(conIndex);
        }

        count++;
        groupIndex.push_back(compEl->Index());

        TPZElementGroup *groupEl = new TPZElementGroup(*cmesh_m);
        elGroups.Push(groupEl);
        elGroups[count - 1]->AddElement(compEl);
    }

    // Inserting interfaces, tangential stress and boundary conditions in the element groups

    for (int64_t el = 0; el < ncompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);

        TPZMultiphysicsInterfaceElement *interEl = dynamic_cast<TPZMultiphysicsInterfaceElement *>(compEl);

        if (interEl)
        {
            TPZCompEl *leftEl = interEl->LeftElement();
            TPZGeoEl *leftGel = leftEl->Reference();
            int left_matid = leftGel->MaterialId();

            if (leftEl->Dimension() != dim)
                continue;

            int64_t leftIndex = leftEl->Index();
            for (int64_t iEl = 0; iEl < groupIndex.size(); iEl++)
            {
                if (leftIndex == groupIndex[iEl])
                {
                    elGroups[iEl]->AddElement(compEl);
                }
            }
        }

        if (!compEl)
            continue;
    }

    cmesh_m->ComputeNodElCon();

    for (auto it = externalNode.begin(); it != externalNode.end(); it++)
    {
        int64_t coIndex = *it;
        cmesh_m->ConnectVec()[coIndex].IncrementElConnected();
    }

    //     Creating  condensed elements
    int64_t nenvel = elGroups.NElements();
    for (int64_t iEnv = 0; iEnv < nenvel; iEnv++)
    {
        TPZElementGroup *elGroup = elGroups[iEnv];
        new TPZCondensedCompElT<STATE>(elGroup);
    }

    cmesh_m->SetName("CMesh_M_Condensed");

    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
}

void InsertBCInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZManVector<int, 2> Interfaces(2, 0);
    Interfaces[0] = simData->InterfaceID();
    Interfaces[1] = -simData->InterfaceID();

    if (!gmesh)
        DebugStop();

    int64_t nel = gmesh->NElements();

    TPZVec<int> IDVec(simData->TangentialBCs().size(), 0); 
    for (int i = 0; i < simData->TangentialBCs().size(); i++)
    {
        IDVec[i] = simData->TangentialBCs()[i].matID;
    }

    // For tangential boundary conditions
    for (auto const &BcMatID : IDVec)
    {
        for (int64_t el = 0; el < nel; el++)
        {
            TPZGeoEl *gel = gmesh->Element(el);
            int meshDim = gmesh->Dimension();
            int matID = gel->MaterialId();

            if (matID != BcMatID)
                continue;

            int nsides = gel->NSides();
            TPZGeoElSide gelSide(gel, nsides - 1);
            TPZCompElSide celSide = gelSide.Reference();

            TPZStack<TPZGeoElSide> neighbourSet;
            gelSide.AllNeighbours(neighbourSet);

            int64_t nneighs = neighbourSet.size();

            TPZManVector<int64_t, 1> LeftElIndex(1, 0), RightElIndex(1, 0);
            LeftElIndex[0] = 0;
            RightElIndex[0] = 1;

            for (int stack_i = 0; stack_i < nneighs; stack_i++)
            {
                TPZGeoElSide neigh = neighbourSet[stack_i];
                int neighMatID = neigh.Element()->MaterialId();
                TPZCompElSide celNeigh = neigh.Reference();

                int64_t neighIndex = neigh.Element()->Index();

                if (neigh.Element()->Dimension() != meshDim)
                    continue;

                if (neigh.Element()->HasSubElement())
                {
                    // Check if it is working in the case with refined meshes
                    DebugStop();
                }
                else
                {
                    TPZGeoElBC gbc(gelSide, Interfaces[0]);

                    TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celNeigh, celSide);
                    interElem->SetLeftRightElementIndices(LeftElIndex, RightElIndex);
                }
            }
        }
    }
}

void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh)
{
    TPZManVector<int64_t, 3> LeftElIndices(1, 0), RightElIndices(1, 1);

    int dim = cmesh_m->Dimension();

    if (!gmesh)
        DebugStop();

    int nInterfaceCreated = 0;

    int mat_lambda = simData->LambdaID();
    int mat_tan = 15;

    int64_t nel = gmesh->NElements();

    for (int64_t el = 0; el < nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        int meshdim = gmesh->Dimension();
        int matid = gel->MaterialId();

        if (gel->Dimension() != meshdim)
            continue;
        if (gel->HasSubElement())
            continue;
        
        int nsides = gel->NSides();

        for (int side = 0; side < nsides; side++)
        {
            if (gel->SideDimension(side) != gmesh->Dimension() - 1)
                continue;
            
            TPZGeoElSide gel_side(gel,side);
            TPZCompElSide cel_side = gel_side.Reference();
            if (!cel_side) DebugStop();

            TPZGeoElSide gel_neighbour = gel_side.Neighbour();
            int neighbour_matid = gel_neighbour.Element()->MaterialId();

            if (neighbour_matid == simData->LambdaID()) //If it is an internal element, the first neighbour should be a tangential stress element due to the way nehgbouring was constructed
            {
                TPZCompElSide cel_neighbour = gel_neighbour.Reference();
                if (!cel_neighbour) DebugStop();

                if (gel_neighbour.Element()->HasSubElement()) DebugStop();

                TPZGeoElBC geo_interface(gel_side, simData->InterfaceID()); //creating the geometric interface element

                TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, geo_interface.CreatedElement(), cel_side, cel_neighbour);
                interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                nInterfaceCreated++;
                
                TPZGeoElSide gel_neighbour2 = gel_neighbour;
                for (; gel_neighbour2 != gel_side; gel_neighbour2++) //Now we look for the neighbour of the neighbour element (tangential stress) until we find the tangential displacement element
                {
                    int neighbour_matid2 = gel_neighbour2.Element()->MaterialId();

                    if (neighbour_matid2 != 15) continue;

                    TPZCompElSide cel_side2 = gel_neighbour2.Reference(); //Tangential displacement element
                    if (!cel_side2) DebugStop();

                    if (gel_neighbour2.Element()->HasSubElement()) DebugStop();

                    TPZGeoElBC geo_interface2(gel_neighbour, simData->InterfaceID()); //creating the geometric interface element

                    TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, geo_interface2.CreatedElement(), cel_side2, cel_neighbour);
                    interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                    nInterfaceCreated++;

                    break;
                }
            }
        }
    }

    std::cout << __PRETTY_FUNCTION__ << "Number of Interfaces Created " << nInterfaceCreated << std::endl;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data)
{
    // TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    // TPZFStructMatrix<STATE> matskl(cmesh);
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

    // int64_t nrow_null_pivot = 1488;
    // TPZFMatrix<REAL> matSingular(nrow_null_pivot+1,nrow_null_pivot+1,0.);
    // auto mat = an.MatrixSolver<STATE>().Matrix();
    // mat->GetSub(0,0,nrow_null_pivot+1,nrow_null_pivot+1,matSingular);
    // matSingular.Print("mat", std::cout, EMathematicaInput);
    // matSingular.PutVal(nrow_null_pivot,nrow_null_pivot,1.);

    // TPZFMatrix<REAL> rhsSingular(nrow_null_pivot+1,1,0.);
    // rhsSingular.PutVal(nrow_null_pivot,0,1.);

    // matSingular.SolveDirect(rhsSingular, ELDLt);

    // TPZFMatrix<REAL> rigidbody(an.Rhs().Rows(),1,0.);
    // rigidbody.PutSub(0,0,rhsSingular);

    // an.LoadSolution(rigidbody);
    // cmesh->TransferMultiphysicsSolution();

    // PrintResults(an, cmesh, problem_data);
    // {
    //     std::ofstream out("cmesh_solve.txt");
    //     cmesh->Print(out);
    // }
    // auto *pardiso = an.MatrixSolver<STATE>().GetPardisoControl();

    // if (pardiso)
    // {
    //     pardiso->SetMessageLevel(1);
    //     constexpr auto sys_type = SymProp::Sym;
    //     constexpr auto prop = TPZPardisoSolver<STATE>::MProperty::EIndefinite;
    //     pardiso->SetMatrixType(sys_type,prop);
    //     TPZManVector<long long, 64> param = pardiso->GetParam();
    //     param[3] = 0;
    //     param[4] = 0;
    //     param[9] = 8;
    //     param[26] = 1;
    //     param[59] = 0;
    //     pardiso->SetParam(param);
    // }
    /// solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    // auto matK = an.MatrixSolver<STATE>().Matrix();
    // auto& rhs = an.Rhs();
    // {
    //     std::ofstream out("system.txt");
    //     (*matK).Print("mat", out, EMathematicaInput); out << std::endl;
    //     rhs.Print("rhs", out, EMathematicaInput); out << std::endl;
    // }
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble() / 1000. << " s" << std::endl;

    // auto res = rhs;
    // matK->MultAdd(an.Solution(), rhs, res, 1.0, -1.0);

    // REAL vecnorm = 0.;
    // TPZFMatrix<STATE>& a = res;
    // an.Solution() = res;
    // an.LoadSolution();
    // cmesh->TransferMultiphysicsSolution();

    // PrintResults(an, cmesh, problem_data);
    // {
    //     std::ofstream out("cmesh_solve.txt");
    //     cmesh->Print(out);
    // }

    // for (int64_t i = 0; i < res.Rows(); i++)
    //     vecnorm += a(i,0)*a(i,0);
    // vecnorm = sqrt(vecnorm);
    // {
    //     std::ofstream out("residual.txt");
    //     res.Print("res", out, EMathematicaInput);
    //     out << std::endl;
    //     out << "VecNorm: " << vecnorm << std::endl;
    // }
    return;
}

void SolveProblemSparseMatRed(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data)
{
    int dimension = cmesh->Dimension();

    // Primeiro cria a matriz auxiliar K00 - que será decomposta
    TPZSYsmpMatrixPardiso<REAL> K00;

    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); // ELU //ECholesky // ELDLt
    an.SetSolver(step);
    step.SetMatrix(&K00);

    // Cria a matriz esparsa
    std::set<int> firstMaterials = {2};
    TPZSparseMatRed<STATE> *matRed = new TPZSparseMatRed<STATE>(cmesh, firstMaterials, TPZSparseMatRed<STATE>::EMaterial);

    {
        std::ofstream out("cmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
        std::ofstream out2("cmesh.txt");
        cmesh->Print(out2);
    }

    std::cout << "Allocating Sub Matrices ...\n";
    
    // Transfere as submatrizes da matriz auxiliar para a matriz correta.
    matRed->SetSolver(&step);
    K00.Resize(matRed->Dim0(), matRed->Dim0());
    matRed->AllocateSubMatrices(cmesh);

    // Compute the number of equations in the system
    int64_t nEqFull = cmesh->NEquations();
    int64_t nEq00 = matRed->Dim0();
    int64_t nEq11 = matRed->Dim1();

    std::cout << "NUMBER OF EQUATIONS:\n "
              << "Full problem = " << nEqFull << ", K00 = " << nEq00 << ", K11 = " << nEq11 << std::endl;

    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull, 1, 0.);
    TPZFMatrix<STATE> rhs1(nEq11, 1, 0.);

    // Creates the problem matrix
    TPZSSpStructMatrix<STATE, TPZStructMatrixOR<STATE>> Stiffness(cmesh);
    Stiffness.SetNumThreads(0);

    TPZAutoPointer<TPZGuiInterface> guiInterface;

    an.SetStructuralMatrix(Stiffness);

    // Monta a matriz
    rhsFull.Zero();
    Stiffness.Assemble(*matRed, rhsFull, guiInterface);

    matRed->SetF(rhsFull);
    matRed->SetReduced();

    // Decomposes the reduced matrix
    matRed->F1Red(rhs1);

    TPZSimpleTimer time_sol;
    TPZSYsmpMatrixPardiso<REAL> K11 = matRed->K11();
    K11.SolveDirect(rhs1, ELDLt);
    std::cout << "Total time = " << time_sol.ReturnTimeDouble() / 1000. << " s" << std::endl;

    TPZFMatrix<STATE> result(nEqFull,1,0.);
    matRed->UGlobal(rhs1, result);

    an.Solution() = result;
    an.LoadSolution();
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, ProblemData *problem_data)
{

    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";

    TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "Stress",
        "Strain",
        "VonMises",
        "Pressure"
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, problem_data->Resolution());
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;

    return;
}
