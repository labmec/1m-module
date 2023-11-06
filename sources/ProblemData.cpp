#include <iostream>
#include <string>
#include<fstream>
#include "json.hpp"
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData(){
    fBcNormalVec.clear(); fBcNormalVec.reserve(10);
    fBcTangentialVec.clear(); fBcTangentialVec.reserve(10);
    fDomain.clear(); fDomain.reserve(10);
}

// deconstructor
ProblemData::~ProblemData(){
    
}

// readjson function. takes a json function as parameter and completes the required simulation data
void ProblemData::ReadJson(std::string file){
    std::ifstream filejson(file);
    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
    
    // checking infos in the json file
    if(input.find("MeshName") == input.end()) DebugStop();
    if(input.find("HdivType")==input.end()) DebugStop();
    if(input.find("DisppOrder") == input.end()) DebugStop();
    if(input.find("Dim") == input.end()) DebugStop();
    if(input.find("Resolution")==input.end()) DebugStop();
    if(input.find("StaticCondensation") == input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
//    if(input.find("NormalBoundary") == input.end()) DebugStop();
//    if(input.find("TangentialBoundary") == input.end()) DebugStop();
        
    // accessing and assigning values
    fMeshName = input["MeshName"];    
   
    fHdivtype = input["HdivType"]; // if hdivtype == -1, then it is H1
    
    fDisppOrder = input["DisppOrder"];
    
    if(input.find("LambdapOrder") != input.end()){
        fLambdapOrder = input["LambdapOrder"];
    }

    if(input.find("CreateMsh") != input.end())
        fMshFile = input["CreateMsh"];
    
    fDim = input["Dim"];
    
    fResolution = input["Resolution"];
    
    fCondensedElement = input["StaticCondensation"];

    if(input.find("InternalPressure") != input.end())
        fInternalPressure = input["InternalPressure"];
    
    DomainData domaindata;
    for(auto& domainjson : input["Domain"]){
        if(domainjson.find("name") == domainjson.end()) DebugStop();
        if(domainjson.find("matID") == domainjson.end()) DebugStop();
        if(domainjson.find("E") == domainjson.end()) DebugStop();
        if(domainjson.find("nu") == domainjson.end()) DebugStop();
        
        domaindata.name = domainjson["name"];
        domaindata.matID = domainjson["matID"];
        domaindata.E = domainjson["E"];
        domaindata.nu = domainjson["nu"];
        
        fDomain.push_back(domaindata);
    }
    
    BcData bcNormaldata;
    for(auto& bcjson : input["NormalBoundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcNormaldata.name = bcjson["name"];
        bcNormaldata.type = bcjson["type"];
        bcNormaldata.matID = bcjson["matID"];
        bcNormaldata.value[0] = bcjson["value"];
        
        fBcNormalVec.push_back(bcNormaldata);
    }
    
    BcData bcTangentialdata;
    for(auto& bcjson : input["TangentialBoundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();

        bcTangentialdata.name = bcjson["name"];
        bcTangentialdata.type = bcjson["type"];
        for (int i = 0; i < fDim; i++)
            bcTangentialdata.value[i] = bcjson["value"][i];
        bcTangentialdata.matID = bcjson["matID"];
        
        fBcTangentialVec.push_back(bcTangentialdata);
    }
    
    if(input.find("InterfaceID") != input.end()) {
        fInterfaceID = input["InterfaceID"];
    }
    if(input.find("LambdaID") != input.end()) {
        fLambdaID = input["LambdaID"];
    }

    fCondensedElement? fMeshVector.resize(4) : fMeshVector.resize(2);
}

void ProblemData::Print(std::ostream& out){
    out << "\nA new simulation has been started: \n\n";
    out << "Mesh Name: " << fMeshName << std::endl << std::endl;
    
    out << "Hdiv Type: " << fHdivtype << std::endl << std::endl;
    
    out << "Velocity pOrder: " << fDisppOrder << std::endl << std::endl;
    
    out << "Traction pOrder: " << fLambdapOrder << std::endl << std::endl;
    
    out << "Dimension: " << fDim << std::endl << std::endl;
    
    out << "Resolution: " << fResolution << std::endl << std::endl;
    
    out << "Static Condensation: " << fCondensedElement << std::endl << std::endl;
    
    out << "Internal Pressure: " << fInternalPressure << std::endl << std::endl;
    
    out << "Domain: " << std::endl;
    
    for(const auto& domaindata : fDomain){
        out << "  Domain Name: " << domaindata.name << std::endl;
        out << "  Domain MatID: " << domaindata.matID << std::endl;
        out << "  Domain E: " << domaindata.E << std::endl;
        out << "  Domain Poisson: " << domaindata.nu << std::endl << std::endl;
    }
    
    out << "Normal Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fBcNormalVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    out << "Tangential Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fBcTangentialVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
    
    for(const auto& bcdata : fBcTangentialVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }

    out << "Interface Elements ID: " << fInterfaceID << std::endl << std::endl;
    
    out << "Lambda Elements ID: " << fLambdaID << std::endl << std::endl;
}
