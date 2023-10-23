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
#include <TPZSimpleTimer.h>
#include "pzvisualmatrix.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include <Elasticity/TPZElasticity3D.h>
#include "TPZAnalyticSolution.h"
#include "TPZGeoMeshTools.h"
#include <TPZGmshReader.h>
#include "tpzchangeel.h"
#include "meshpath_config.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "ProblemData.h"

enum EMatid {ENone,EDomain,EBC,EinternalBC,ESymmetryBC,EFixedXZ,EFixedZ,EStiffner};
const int global_nthread = 16;

TPZGeoMesh* CreateGMesh(int ndiv);
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);
void CreateBCs(TPZGeoMesh* gmesh);
void ChangeElsToCylMap(TPZGeoMesh* gmesh);
void RefinePyrTo2Tets(TPZGeoMesh* gmesh);
void printVTKWJacInfo(std::string filename, TPZGeoMesh* gmesh);
TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas, const REAL& internalPressure);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

#ifdef PZ_LOG
static TPZLogger logger("pz.1mmodule");
#endif


int main(int argc, char *argv[]) {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "--------- Starting simulation ---------" << std::endl;
    const bool isUseDirectionalRef = false;
    const bool isUseCylMap = true;
    const int nUniformRef = 0;
    if(isUseDirectionalRef){
        gRefDBase.InitializeRefPatterns();
    }
    
    // Reading problem data from json
    std::string jsonfilename = "1m-module-init.json";
//    std::string jsonfilename = "1m-module-lid.json";
    if(argc > 1) jsonfilename = std::string(argv[1]);
    ProblemData problemdata;
    problemdata.ReadJson(std::string(MESHES_DIR) + "/" + jsonfilename);
        
    // Create gmesh
    const int pord = problemdata.DisppOrder();
    const bool readGMeshFromGmsh = true;
    TPZGeoMesh* gmesh = nullptr;
    if(readGMeshFromGmsh){
        std::string filename = problemdata.MeshName();
        gmesh = ReadMeshFromGmsh(std::string(MESHES_DIR) + "/" + filename);
        CreateBCs(gmesh);
    }
    else{ // use TPZGenGrid
        int ndiv = 2;
        gmesh = CreateGMesh(ndiv);
    }
    
    
    if(isUseCylMap){
        printVTKWJacInfo("gmesh_jac_before_cyl.vtk",gmesh);
        ChangeElsToCylMap(gmesh);
        printVTKWJacInfo("gmesh_jac_after_cyl.vtk",gmesh);
    }
       
    // Refine mesh directionally towards the singularity at the stiffner
    if(isUseDirectionalRef){
        std::set<int> matidstoref = {EStiffner};
        TPZRefPatternTools::RefineDirectional(gmesh, matidstoref);
        TPZRefPatternTools::RefineDirectional(gmesh, matidstoref);
        TPZRefPatternTools::RefineDirectional(gmesh, matidstoref);
    }
    

#ifdef PZ_LOG
//    if (logger.isDebugEnabled()) {
//        for(int iel = 0 ; iel < gmesh->NElements() ; iel++){
//            TPZGeoEl* geoel = gmesh->Element(iel);
//            TPZFNMatrix<9,REAL> gradx(3,3,0.);
//            TPZManVector<REAL,3> qsicenter(geoel->Dimension(),0.);
//            geoel->CenterPoint(geoel->NSides()-1, qsicenter);
//            geoel->GradX(qsicenter, gradx);
//            std::stringstream sout;
//            sout << "el = " << iel << std::endl;
//            gradx.Print(sout);
//            LOGPZ_DEBUG(logger, sout.str())
//        }
//    }
#endif
    
    if(isUseDirectionalRef){
        RefinePyrTo2Tets(gmesh);
    }
        
    if(nUniformRef){
        TPZCheckGeom geom(gmesh);
        geom.UniformRefine(nUniformRef);
    }
    
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);
    }
    
    // Create compmeshes
    if(problemdata.DomainVec().size() > 1) DebugStop(); // Please implement the next lines correctly if many domains
    
    TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
    elas->fE = problemdata.DomainVec()[0].E;
    elas->fPoisson = problemdata.DomainVec()[0].nu;
    elas->fProblemType = TElasticity3DAnalytic::ENone;
    
    TPZCompMesh* cmesh = nullptr;
    if(problemdata.HdivType() == 0){
        cmesh = CreateH1CMesh(gmesh,pord,elas,problemdata.InternalPressure());
    }
    else{
        std::cout << "Implement me!" << std::endl;
    }
    
    // Analysis
    //Solve Multiphysics
    TPZLinearAnalysis an(cmesh);
    an.SetExact(elas->ExactSolution());
    SolveProblemDirect(an,cmesh);
    
    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an,cmesh);

    
    // deleting stuff
    delete cmesh;
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
//        stringtoint[1]["stif"] = EStiffner;
        
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh,false);
    }
    
    return gmesh;
}

void CreateBCs(TPZGeoMesh* gmesh) {
//    TPZManVector<int64_t,1> cornerindexes = {24};
//    TPZManVector<int64_t,1> cornerindexes = {351};
//    TPZManVector<int64_t,1> cornerindexes = {0};
//    TPZManVector<int64_t,1> cornerindexes = {140};
    TPZManVector<int64_t,1> cornerindexes = {142};
    int64_t index = -1;
    gmesh->CreateGeoElement(EPoint, cornerindexes, EFixedXZ, index);
//    cornerindexes = {27};
//    cornerindexes = {352};
//    cornerindexes = {2};
//    cornerindexes = {141};
    cornerindexes = {143};
    gmesh->CreateGeoElement(EPoint, cornerindexes, EFixedZ, index);
    gmesh->BuildConnectivity();
    
    const int64_t nel = gmesh->NElements();
    int count = 0;
    for(int64_t iel = 0 ; iel < nel ; iel++){
        TPZGeoEl* gel = gmesh->Element(iel);
        if(gel->Dimension() == 1) continue;
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


TPZCompMesh* CreateH1CMesh(TPZGeoMesh* gmesh, const int pord, TElasticity3DAnalytic *elas, const REAL& internalPressure) {
    TPZCompMesh* cmesh = new TPZCompMesh(gmesh);
    const int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pord);
    cmesh->SetAllCreateFunctionsContinuous();
    
    // Domain elas mat
    const STATE E = elas->fE, nu = elas->fPoisson;
    TPZManVector<STATE> force = {0,0,0};
    TPZElasticity3D *mat = new TPZElasticity3D(EDomain, E, nu, force, 0., 0., 0.);
//    mat->SetExactSol(elas->ExactSolution(), 2);
//    mat->SetForcingFunction(elas->ForceFunc(), 2);
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
    val1 *= internalPressure;
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

    // TPZSkylineStructMatrix<STATE> matskl(cmesh);
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
    constexpr int vtkRes{2};
    
    TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "Stress",
        "Strain",
        "PrincipalStrain",
        "VonMises",
        "I2"
    };
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    
    return;
}

void ChangeElsToCylMap(TPZGeoMesh* gmesh) {
    const int64_t nel = gmesh->NElements();
    TPZManVector<REAL,3> xcenter(3,0.);
    TPZFNMatrix<9,REAL> axis(3,3,0.);
    axis.Identity();
    for(int64_t iel = 0 ; iel < nel ; iel++){
        TPZGeoEl* geoel = gmesh->Element(iel);
        TPZChangeEl::ChangeToCylinder(gmesh, iel, xcenter, axis);
    }
}

void printVTKWJacInfo(std::string filename, TPZGeoMesh* gmesh) {
    TPZVec<REAL> elData(gmesh->NElements(), -100);
    for (int i = 0; i < gmesh->NElements(); i++) {
        TPZGeoEl* gel = gmesh->Element(i);
        if(!gel) DebugStop();
        const int geldim = gel->Dimension();
        TPZManVector<REAL,3> qsi(geldim,0.);
        gel->CenterPoint(gel->NSides()-1, qsi);
        REAL detjac = -1000;
        TPZFMatrix<REAL> jac(3,3,0.), axes(3,3,0.), jacinv(3,3,0.);
        gel->Jacobian(qsi, jac, axes, detjac, jacinv);
        elData[i] = detjac;
    }
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, elData);
}

void RefinePyrTo2Tets(TPZGeoMesh* gmesh) {
    auto refpat = gRefDBase.FindRefPattern("PyrTwoTets");
    if(!refpat) DebugStop();
    for(auto& gel : gmesh->ElementVec()){
        if(gel->HasSubElement()) continue;
        if (gel->Type() == EPiramide) {
            gel->SetRefPattern(refpat);
            TPZVec<TPZGeoEl*> subels;
            gel->Divide(subels);
        }
    }
}
