#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
{
    // struct responsible to summarize all the data from every domain
    struct DomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        REAL E = -1.; // domain young modulus
        REAL nu = -1.; // poisson
    };
    
    // struct responsible to store boundary condition data
    struct BcData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<REAL,3>  value = {0., 0., 0.}; // bc value
        int matID = 0; // bc material ID
    };

private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none";
    
    bool fMshFile = false;
    
    int fHdivtype = 0; // 0 -> Standard & 1 -> Constant
    
    int fDisppOrder = -1; // polynomial approximation order for velocity
    
    int fLambdapOrder = -1; // polynomial approximation order for traction
    
    int fDim = -1;
    
    int fResolution = -1;

    bool fCondensedElement = false;
    
    std::vector<DomainData> fDomain; // vector containing every domain created

    std::vector<BcData> fBcNormalVec; // vector containg all the velocity bcs info
    
    std::vector<BcData> fBcTangentialVec; // vector containg all the traction bcs info
    
    int fInterfaceID = 20;
    
    int fLambdaID = 10;
    
    REAL fInternalPressure = -1;
    
    TPZVec<TPZCompMesh*> fMeshVector;
    
public:
    ProblemData();
    
    ~ProblemData();
    
    void ReadJson(std::string jsonfile);
    
    void Print(std::ostream& out = std::cout);
    // you can pass a file, so that the simulation data will be printed inside it. Otherwise, it will be displayed on terminal
    
    const std::string& MeshName() const {return fMeshName;} //setter using reference variable;
    void SetMeshName(const std::string& meshname) {fMeshName = meshname;}  //getter using reference variable
    
    const bool& CreateMshFile() const {return fMshFile;}
    void SetCreateMshFile(bool create) {fMshFile = create;}
    
    const int& HdivType() const {return fHdivtype;}
    void SetHdivType(int hdivtype) {fHdivtype = hdivtype;}
    
    const int& DisppOrder() const {return fDisppOrder;}
    void SetDisppOrder(int velp) {fDisppOrder = velp;}
    
    const int& LambdapOrder() const {return fLambdapOrder;}
    void SetLambdapOrder(int tracp) {fLambdapOrder = tracp;}
    
    const int& Dim() const {return fDim;}
    void SetDim(int dim) {fDim = dim;}
    
    const int& Resolution() const {return fResolution;}
    void SetResolution(int res) {fResolution = res;}
    
    const bool& CondensedElements() const {return fCondensedElement;}
    void SetCondensedElements(bool condense) {fCondensedElement = condense;}

    const REAL& InternalPressure() const {return fInternalPressure;}
    void SetInternalPressure(bool intpressure) {fInternalPressure = intpressure;}
    
    const std::vector<DomainData>& DomainVec() const {return fDomain;}
    void SetDomainVec(const std::vector<DomainData>& vec) {fDomain = vec;}
    
    const std::vector<BcData>& NormalBCs() const {return fBcNormalVec;}
    void SetNormalBCs(const std::vector<BcData>& bcs) {fBcNormalVec = bcs;}
    
    const std::vector<BcData>& TangentialBCs() const {return fBcTangentialVec;}
    void SetTangentialBCs(const std::vector<BcData>& bcs) {fBcTangentialVec = bcs;}
    
    const int& InterfaceID() const{return fInterfaceID;}
    void SetInterfaceID(int id) {fInterfaceID = id;}
    
    const int& LambdaID() const{return fLambdaID;}
    void SetLambdaID(int id ){fLambdaID = id;}
    
    TPZVec<TPZCompMesh*>& MeshVector() {return fMeshVector;}
    void SetMeshVector(const TPZVec<TPZCompMesh*>& vec) {fMeshVector = vec;}
};

#endif
