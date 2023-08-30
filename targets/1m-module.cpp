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
#include "TPZMixedModelProblem.h"
#include <TPZNullMaterial.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzskylstrmatrix.h>
#include <TPZMultiPhysicsCompMesh.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include <Elasticity/TPZElasticity3D.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "meshpath_config.h"

enum EMatid {ENone,EDomain,EBC,EinternalBC,ESymmetryBC,EFixedXZ,EFixedZ};
const int global_nthread = 16;

TPZGeoMesh* CreateGMesh(int ndiv);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);


int main() {
#ifdef PZ_LOG
//    TPZLogger::InitializePZLOG("log4cxx.cfg");
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    // Create gmesh
    int ndiv = 2;
    const int pord = 3;
    const bool readGMeshFromGmsh = true;
    TPZGeoMesh* gmesh = nullptr;
    if(readGMeshFromGmsh){
//        std::string filename = "geometry_shell_test.msh";
        std::string filename = "geometry_shell_new2.msh";
        gmesh = ReadMeshFromGmsh(std::string(MESHES_DIR) + "/" + filename);
        CreateBCs(gmesh);
    }
    else{
        gmesh = CreateGMesh(ndiv);
    }
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);
    }
    
    // Create compmeshes
    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    elas->fE = 250.;//206.8150271873455;
    elas->fPoisson = 0.;
    elas->fProblemType = TElasticity3DAnalytic::EStretchx;
    TPZCompMesh* cmeshH1 = CreateH1CMesh(gmesh,pord,elas);
    
    // Analysis
    //Solve Multiphysics
    TPZLinearAnalysis an(cmeshH1);
    an.SetExact(elas->ExactSolution());
    SolveProblemDirect(an,cmeshH1);
    
    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an,cmeshH1);

    
    // deleting stuff
    delete cmeshH1;
    delete gmesh;
    
    
    std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh* CreateGMesh(int ndiv) {
    TPZGeoMesh* gmesh = new TPZGeoMesh;
//    const TPZManVector<REAL, 3> x0 = {0., 0., 0.};
//    const TPZManVector<REAL, 3> x1 = {1., 1., 1.};
//    const TPZManVector<int, 3> ndivvec = {ndiv, ndiv, ndiv};
//    TPZGenGrid2D gen(ndivvec, x0, x1);
//    gen.Read(gmesh);
//
//    gen.SetBC(gmesh, 4, EBC); // bot
//    gen.SetBC(gmesh, 5, EBC); // right
//    gen.SetBC(gmesh, 6, EBC); // top
//    gen.SetBC(gmesh, 7, EBC); // left
//
//    return gmesh;
    
    MMeshType meshType = MMeshType::EHexahedral;
    int dim = 3;
    
    TPZManVector<REAL,3> minX = {-1,-1,-1};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;
    
    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,EBC);
    matIds[0] = EDomain;
    // matIds[1] = bcId;
    // matIds[2] = EBoundary1;
    // matIds[3] = EBoundary1;
    // matIds[4] = EBoundary1;
    
    TPZManVector<int,3> ndivvec = {ndiv,ndiv,ndiv};
    gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,matIds, ndivvec, meshType,createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);
    
    return gmesh;
}

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[3]["dom"] = EDomain;
//        stringtoint[2]["Surfaces"] = EBC;
        
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);
    }
    
    return gmesh;
}

void CreateBCs(TPZGeoMesh* gmesh) {
//    TPZManVector<int64_t,1> cornerindexes = {24};
    TPZManVector<int64_t,1> cornerindexes = {351};
    int64_t index = -1;
    gmesh->CreateGeoElement(EPoint, cornerindexes, EFixedXZ, index);
//    cornerindexes = {27};
    cornerindexes = {352};
    gmesh->CreateGeoElement(EPoint, cornerindexes, EFixedZ, index);
    gmesh->BuildConnectivity();
    
    const int64_t nel = gmesh->NElements();
    int count = 0;
    for(int64_t iel = 0 ; iel < nel ; iel++){
        TPZGeoEl* gel = gmesh->Element(iel);
        const int firstside = gel->FirstSide(2);
        const int lastside = gel->FirstSide(3);
        for(int iside = firstside ; iside < lastside ; iside++) {
            TPZGeoElSide gelside(gel,iside);
            if(gelside.Neighbour() == gelside){
                TPZManVector<REAL,3> centerX(3,0.), center(2,0.), normal(3,0.);
                gelside.CenterX(centerX);
                gelside.CenterPoint(center);
                REAL radius = sqrt(centerX[0]*centerX[0] + centerX[1]*centerX[1]);
                gelside.Normal(center, normal);
                if(radius < 343.1 && fabs(normal[2]) < 0.2){
                    TPZGeoElBC(gelside, EinternalBC);
                }
                else if(fabs(centerX[1]) < 1.e-3){
                    TPZGeoElBC(gelside, ESymmetryBC);
//                    if(count == 0){
//                        const int nnod = gelside.NSideNodes();
//                        if(nnod != 3) DebugStop();
//                        TPZGeoElBC(gel,gelside.SideNodeLocIndex(0),EFixedXZ);
//                        TPZGeoElBC(gel,gelside.SideNodeLocIndex(1),EFixedX);
//                        count++;
//                    }
                }
                else{
                    TPZGeoElBC(gelside, EBC);
                }
            }
        }
    }
    
}


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Domain elas mat
    const STATE E = elas->fE, nu = elas->fPoisson;
    TPZManVector<STATE> force = {0,0,0};
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
    mat->SetExactSol(elas->ExactSolution(), 2);
    mat->SetForcingFunction(elas->ForceFunc(), 2);
    cmesh->InsertMaterialObject(mat);
//    mat->SetDimension(dim);
    
    // BC null mat
    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    
    const int diri = 0, neu = 1, mixed = 2, normaltrac = 4;
    auto* BCCond0 = mat->CreateBC(mat, EBC, neu, val1, val2);
//    BCCond0->SetForcingFunctionBC(elas->ExactSolution(), 2);
    cmesh->InsertMaterialObject(BCCond0);
    
    val1(1,1) = mat->BigNumber();
    auto* BCCondSymm = mat->CreateBC(mat, ESymmetryBC, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCCondSymm);
    
    val1.Identity();
    auto* BCCondInt = mat->CreateBC(mat, EinternalBC, normaltrac, val1, val2);
    cmesh->InsertMaterialObject(BCCondInt);

    val1.Zero();
    val1(0,0) = mat->BigNumber();
    val1(2,2) = mat->BigNumber();
    auto* BCCondXZ = mat->CreateBC(mat, EFixedXZ, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCCondXZ);

    val1.Zero();
    val1(2,2) = mat->BigNumber();
    auto* BCCondZ = mat->CreateBC(mat, EFixedZ, mixed, val1, val2);
    cmesh->InsertMaterialObject(BCCondZ);
    
    // Constructs mesh
    cmesh->AutoBuild();
    
    return cmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

//    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    TPZSimpleTimer time_ass;
    an.Assemble();
    std::cout << "Total time = " << time_ass.ReturnTimeDouble()/1000. << " s" << std::endl;
    
//    extern TPZManVector<STATE,3> integratedforce;
//    std::cout << "\nintegratedforce = " << integratedforce << std::endl;
    
    ///solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";
    constexpr int vtkRes{1};
    
    TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "Stress",
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

