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
#include <TPZMultiphysicsCompMesh.h>
#include <pzstepsolver.h>
#include <TPZLinearAnalysis.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
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

const int global_nthread = 16;

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
void CreateBCs(TPZGeoMesh *gmesh, const ProblemData *problem_data);
void InsertLambda(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshU(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshP(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshG(ProblemData *simData, TPZGeoMesh *gmesh);
TPZCompMesh *CreateCMeshPm(ProblemData *simData, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh);
void CondenseElements(ProblemData *simData, TPZMultiphysicsCompMesh *cmesh_m);
void InsertBCInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh);
void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m, ProblemData *simData, TPZGeoMesh *gmesh);
void printVTKWJacInfo(std::string filename, TPZGeoMesh *gmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

#ifdef PZ_LOG
static TPZLogger logger("pz.1mmodule");
#endif

int main(int argc, char *argv[])
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    std::cout << "--------- Starting simulation ---------" << std::endl;
    gRefDBase.InitializeRefPatterns();

    // Reading problem data from json
    std::string jsonfilename = "1m-module-init.json";
    if (argc > 1)
        jsonfilename = std::string(argv[1]);
    ProblemData problemdata;
    problemdata.ReadJson(std::string(MESHES_DIR) + "/" + jsonfilename);

    // Create gmesh
    const int pord = problemdata.DisppOrder();
    TPZGeoMesh *gmesh = nullptr;
    std::string filename = problemdata.MeshName();
    gmesh = ReadMeshFromGmsh(std::string(MESHES_DIR) + "/" + filename, &problemdata);

    CreateBCs(gmesh, &problemdata);

    gmesh->ElementVec()[10440]->SetMaterialId(EZeroNormalDisp);
    gmesh->ElementVec()[10046]->SetMaterialId(EZeroNormalDisp);

    for (auto gel : gmesh->ElementVec())
    {
        if (gel->Dimension() != 2)
            continue;
        
        int matid = gel->MaterialId();
        int64_t nnode = gel->NNodes();
        for (int i = 0; i < nnode; i++)
        {
            int64_t node_index = gel->NodeIndex(i);
            if (node_index == 142)
            {
                TPZGeoElSide geoside(gel);
                TPZVec<REAL> point(2,1./3.);
                TPZVec<REAL> normal(3,0.);
                geoside.Normal(point,normal);
                if (abs(abs(normal[2]) - 1.0) <= 1.0e-3)
                {
                    std::cout << "gel with node 142 is: " << gel->Index() << ", Normal: " << normal << ", matid: " << matid << std::endl;
                }

            }
            else if (node_index == 143)
            {
                TPZGeoElSide geoside(gel);
                TPZVec<REAL> point(2,1./3.);
                TPZVec<REAL> normal(3,0.);
                geoside.Normal(point,normal);
                if (abs(abs(normal[2]) - 1.0) <= 1.0e-3)
                {
                    std::cout << "gel with node 143 is: " << gel->Index() << ", Normal: " << normal << ", matid: " << matid << std::endl;
                }
            }
        }
    }

    printVTKWJacInfo("gmesh_jac_before_cyl.vtk", gmesh);

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

    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);
    }

    // Create compmeshes
    if (problemdata.DomainVec().size() > 1)
        DebugStop(); // Please implement the next lines correctly if many domains

    TPZCompMesh *cmesh_u = CreateCMeshU(&problemdata, gmesh);
    {
        std::ofstream out("cmesh_u.txt");
        cmesh_u->Print(out);
    }
    TPZCompMesh *cmesh_p = CreateCMeshP(&problemdata, gmesh);
    TPZCompMesh *cmesh_pm = nullptr;
    TPZCompMesh *cmesh_g = nullptr;
    if (problemdata.CondensedElements())
    {
        cmesh_pm = CreateCMeshPm(&problemdata, gmesh);
        cmesh_g = CreateCMeshG(&problemdata, gmesh);
    }

    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiphysicsMesh(&problemdata, gmesh);

    if (problemdata.CondensedElements())
        CondenseElements(&problemdata, cmesh_m);

    // Analysis
    // Solve Multiphysics
    TPZLinearAnalysis an(cmesh_m);
    SolveProblemDirect(an, cmesh_m);

    // Post Process
    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh_m);

    // deleting stuff
    if (cmesh_u)
        delete cmesh_u;
    if (cmesh_p)
        delete cmesh_p;
    if (cmesh_g)
        delete cmesh_g;
    if (cmesh_pm)
        delete cmesh_pm;
    if (cmesh_m)
        delete cmesh_m;

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
    InsertLambda(problem_data, gmesh);
    gmesh->BuildConnectivity();
    return gmesh;
}

void InsertLambda(ProblemData *simData, TPZGeoMesh *gmesh)
{
    int64_t nEl = gmesh->NElements();

    // We insert a lambda element between two neighbour elements
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

                while (neighbour != geoElSide)
                {

                    if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1)
                    {
                        int neighbourMatId = neighbour.Element()->MaterialId();

                        break;
                    }

                    neighbour = neighbour.Neighbour();
                }

                if (neighbour == geoElSide)
                    TPZGeoElBC(geoElSide, 10);
            }
        }
    }
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

        auto *mat = new TPZNullMaterial<>(simData->DomainVec()[0].matID);
        cmesh_u->InsertMaterialObject(mat);

        materialIDs.insert(simData->DomainVec()[0].matID);

        // boundary conditions' material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;

            auto BCmat = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            cmesh_u->InsertMaterialObject(BCmat);
            materialIDs.insert(bc.matID);
        }

        val2[0] = simData->InternalPressure();
        auto BCmat = mat->CreateBC(mat, EPressure, 2, val1, val2);
        cmesh_u->InsertMaterialObject(BCmat);
        materialIDs.insert(EPressure);

        val2[0] = 0.;
        auto BCmat2 = mat->CreateBC(mat, EZeroNormalDisp, 0, val1, val2);
        cmesh_u->InsertMaterialObject(BCmat2);
        materialIDs.insert(EZeroNormalDisp);

        auto BCmat3 = mat->CreateBC(mat, EZeroNormalStress, 2, val1, val2);
        cmesh_u->InsertMaterialObject(BCmat3);
        materialIDs.insert(EZeroNormalStress);

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
        gmesh->ResetReference();
    }

    // expanding the solution vector
    cmesh_u->ExpandSolution();

    if ((simData->CondensedElements() && simData->MeshVector().size() != 4) || simData->MeshVector().size() < 2)
        DebugStop();

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

        TPZNullMaterial<>* matLambdaBC = new TPZNullMaterial<>(EZeroTangentialStress);
        matLambdaBC->SetNStateVariables(simData->Dim() - 1);
        cmesh_p->InsertMaterialObject(matLambdaBC);
        materialIDs.insert(EZeroTangentialStress);

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

        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }

        gmesh->ResetReference();
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

    simData->MeshVector()[3] = cmesh_g;

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

    simData->MeshVector()[2] = cmesh_pm;

    return cmesh_pm;
}

TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(ProblemData *simData, TPZGeoMesh *gmesh)
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

        TPZHybridMixedElasticityUP *mat = new TPZHybridMixedElasticityUP(simData->DomainVec()[0].matID, simData->Dim(), young, poisson, AnalysisType::EGeneral);

        cmesh_m->InsertMaterialObject(mat);

        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        for (const auto &bc : simData->NormalBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            cmesh_m->InsertMaterialObject(matBC);
        }

        //This is hard coded to impose normal displacement and stress
        {
            val2[0] = -1.*simData->InternalPressure();
            auto BCmat = mat->CreateBC(mat, EPressure, 2, val1, val2);
            auto matBC = dynamic_cast<TPZBndCondT<STATE> *>(BCmat);
            cmesh_m->InsertMaterialObject(matBC);

            val2[0] = 0.;
            auto BCmat2 = mat->CreateBC(mat, EZeroNormalDisp, 0, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(BCmat2);
            cmesh_m->InsertMaterialObject(matBC2);

            auto BCmat3 = mat->CreateBC(mat, EZeroNormalStress, 0, val1, val2);
            auto matBC3 = dynamic_cast<TPZBndCondT<STATE> *>(BCmat3);
            cmesh_m->InsertMaterialObject(matBC3);
        }
        
        for (const auto &bc : simData->TangentialBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            cmesh_m->InsertMaterialObject(matBC);
        }

        //This is hard coded to imposed tangential stress
        if (0)
        {
            val2[0] = 0.;
            auto BCmat = mat->CreateBC(mat, EZeroTangentialStress, 3, val1, val2);
            auto matBC = dynamic_cast<TPZBndCondT<STATE> *>(BCmat);
            cmesh_m->InsertMaterialObject(matBC);
        }

        // 3 - Material for tangential traction
        TPZNullMaterialCS<> *matLambda = new TPZNullMaterialCS<>(simData->LambdaID());
        double dim = simData->Dim() - 1;
        matLambda->SetDimension(dim);
        matLambda->SetNStateVariables(dim);
        cmesh_m->InsertMaterialObject(matLambda);

        // 4 - Material for interfaces (Inner)
        TPZMatInterfaceHybridElasticityStokes *matInterfaceLeft = new TPZMatInterfaceHybridElasticityStokes(simData->InterfaceID(), simData->Dim() - 1);
        matInterfaceLeft->SetMultiplier(1.);
        cmesh_m->InsertMaterialObject(matInterfaceLeft);

        TPZMatInterfaceHybridElasticityStokes *matInterfaceRight = new TPZMatInterfaceHybridElasticityStokes(-simData->InterfaceID(), simData->Dim() - 1);
        matInterfaceRight->SetMultiplier(-1.);
        cmesh_m->InsertMaterialObject(matInterfaceRight);
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
    TPZStack<TPZElementGroup *> elGroups;
    int count = 0;

    auto cmesh_u = cmesh_m->MeshVector()[0];

    // Creating the element groups for the domain
    for (int64_t el = 0; el < ncompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);

        if (compEl->Dimension() != dim)
            continue;

        TPZMultiphysicsElement *multEl = dynamic_cast<TPZMultiphysicsElement *>(compEl);
        int64_t numSpaces = multEl->NMeshes();

        if (numSpaces < 4)
            DebugStop();

        int64_t numConnectExt = numSpaces - 3;
        int nConnect = multEl->NConnects();

        for (int ic = nConnect - numConnectExt; ic < nConnect; ic++)
        {
            int64_t conIndex = compEl->ConnectIndex(ic);
            externalNode.insert(conIndex);
        }

        count++;
        groupIndex.resize(count);
        groupIndex[count - 1] = compEl->Index();

        TPZElementGroup *groupEl = new TPZElementGroup(*cmesh_m);
        elGroups.Push(groupEl);
        elGroups[count - 1]->AddElement(compEl);
    }

    // Inserting interfaces and boundary conditions

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
    TPZManVector<int, 2> Interfaces(2, 0);
    Interfaces[0] = simData->InterfaceID();
    Interfaces[1] = -simData->InterfaceID();

    int dim = cmesh_m->Dimension();

    if (!gmesh)
        DebugStop();

    int nInterfaceCreated = 0;

    int matfrom = simData->LambdaID();

    int64_t nel = gmesh->NElements();

    for (int64_t el = 0; el < nel; el++)
    {
        TPZGeoEl *gel = gmesh->Element(el);
        int meshdim = gmesh->Dimension();
        int matid = gel->MaterialId();

        if (matid != matfrom)
            continue;
        if (gel->HasSubElement())
            continue;

        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel, nsides - 1);
        TPZCompElSide celside = gelside.Reference();

        TPZStack<TPZGeoElSide> neighbourSet;
        gelside.AllNeighbours(neighbourSet);

        gelside.LowerLevelCompElementList2(1);

        int64_t nneighs = neighbourSet.size();

        TPZManVector<int64_t, 3> LeftElIndices(1, 0.), RightElIndices(1, 0);
        LeftElIndices[0] = 0;
        RightElIndices[0] = 1;

        for (int stack_i = 0; stack_i < nneighs; stack_i++)
        {
            TPZGeoElSide neigh = neighbourSet[stack_i];
            int neighMatID = neigh.Element()->MaterialId();
            TPZCompElSide celneigh = neigh.Reference();

            if (!celside)
                DebugStop();

            if (neigh.Element()->HasSubElement())
            {
                DebugStop();
            }
            else
            {

                int64_t neighIndex = neigh.Element()->Index();

                if (neigh.Element()->Dimension() != meshdim)
                    continue;

                TPZGeoElBC gbc(gelside, Interfaces[stack_i]);

                TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh_m, gbc.CreatedElement(), celneigh, celside);
                interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                nInterfaceCreated++;
            }
        }
    }

    std::cout << __PRETTY_FUNCTION__ << "Number of Interfaces Created " << nInterfaceCreated << std::endl;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

    //    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE> matskl(cmesh);
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

    //    extern TPZManVector<STATE,3> integratedforce;
    //    std::cout << "\nintegratedforce = " << integratedforce << std::endl;

    /// solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    TPZSimpleTimer time_sol;
    an.Solve();
    std::cout << "Total time = " << time_sol.ReturnTimeDouble() / 1000. << " s" << std::endl;

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
        "I2"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;

    return;
}

// void ChangeElsToCylMap(TPZGeoMesh* gmesh) {
//     const int64_t nel = gmesh->NElements();
//     TPZManVector<REAL,3> xcenter(3,0.);
//     TPZFNMatrix<9,REAL> axis(3,3,0.);
//     axis.Identity();
//     for(int64_t iel = 0 ; iel < nel ; iel++){
//         TPZGeoEl* geoel = gmesh->Element(iel);
//         TPZChangeEl::ChangeToCylinder(gmesh, iel, xcenter, axis);
//     }
// }

void printVTKWJacInfo(std::string filename, TPZGeoMesh *gmesh)
{
    TPZVec<REAL> elData(gmesh->NElements(), -100);
    for (int i = 0; i < gmesh->NElements(); i++)
    {
        TPZGeoEl *gel = gmesh->Element(i);
        if (!gel)
            DebugStop();
        const int geldim = gel->Dimension();
        TPZManVector<REAL, 3> qsi(geldim, 0.);
        gel->CenterPoint(gel->NSides() - 1, qsi);
        REAL detjac = -1000;
        TPZFMatrix<REAL> jac(3, 3, 0.), axes(3, 3, 0.), jacinv(3, 3, 0.);
        gel->Jacobian(qsi, jac, axes, detjac, jacinv);
        elData[i] = detjac;
    }
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, elData);
}

void RefinePyrTo2Tets(TPZGeoMesh *gmesh)
{
    auto refpat = gRefDBase.FindRefPattern("PyrTwoTets");
    if (!refpat)
        DebugStop();
    for (auto &gel : gmesh->ElementVec())
    {
        if (gel->HasSubElement())
            continue;
        if (gel->Type() == EPiramide)
        {
            gel->SetRefPattern(refpat);
            TPZVec<TPZGeoEl *> subels;
            gel->Divide(subels);
        }
    }
}

void CreateBCs(TPZGeoMesh *gmesh, const ProblemData *problem_data)
{
    //    TPZManVector<int64_t,1> cornerindexes = {24};
    //    TPZManVector<int64_t,1> cornerindexes = {351};
    //    TPZManVector<int64_t,1> cornerindexes = {0};
    //    TPZManVector<int64_t,1> cornerindexes = {140};
    TPZManVector<int64_t, 1> cornerindexes = {142};
    int64_t index = -1;
    // gmesh->CreateGeoElement(EPoint, cornerindexes, EFixedXZ, index);
    //    cornerindexes = {27};
    //    cornerindexes = {352};
    //    cornerindexes = {2};
    //    cornerindexes = {141};
    cornerindexes = {143};
    // gmesh->CreateGeoElement(EPoint, cornerindexes, EFixedZ, index);
    gmesh->BuildConnectivity();

    const int64_t nel = gmesh->NElements();
    int count = 0;
    for (int64_t iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = gmesh->Element(iel);
        int matid = gel->MaterialId();
        if (gel->Dimension() == 1)
            continue;
        if (matid == problem_data->LambdaID())
            continue;
        const int firstside = gel->FirstSide(2);
        const int lastside = gel->FirstSide(3);
        for (int iside = firstside; iside < lastside; iside++)
        {
            TPZGeoElSide gelside(gel, iside);
            if (gelside.Neighbour() == gelside)
            {
                TPZManVector<REAL, 3> centerX(3, 0.), center(2, 0.), normal(3, 0.);
                gelside.CenterX(centerX);
                gelside.CenterPoint(center);
                REAL radius = sqrt(centerX[0] * centerX[0] + centerX[1] * centerX[1]);
                gelside.Normal(center, normal);
                if (radius < 343.1 && fabs(normal[2]) < 0.2)
                {
                    TPZGeoElBC(gelside, EPressure);
                    TPZGeoElBC(gelside, EZeroTangentialStress);
                }
                else if (fabs(centerX[1]) < 1.e-3)
                {
                    TPZGeoElBC(gelside, EZeroNormalDisp);
                    TPZGeoElBC(gelside, EZeroTangentialStress);
                    //                    if(count == 0){
                    //                        const int nnod = gelside.NSideNodes();
                    //                        if(nnod != 3) DebugStop();
                    //                        TPZGeoElBC(gel,gelside.SideNodeLocIndex(0),EFixedXZ);
                    //                        TPZGeoElBC(gel,gelside.SideNodeLocIndex(1),EFixedX);
                    //                        count++;
                    //                    }
                }
                else
                {
                    TPZGeoElBC(gelside, EZeroNormalStress);
                    TPZGeoElBC(gelside, EZeroTangentialStress);
                }
            }
        }
    }
}