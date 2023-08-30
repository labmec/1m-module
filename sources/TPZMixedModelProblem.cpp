#include "TPZMixedModelProblem.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"


TPZMixedModelProblem::TPZMixedModelProblem() : TBase() {
}


TPZMixedModelProblem::TPZMixedModelProblem(int mat_id, int dimension) : TBase(mat_id, dimension) {
    
    
}


TPZMixedModelProblem::TPZMixedModelProblem(const TPZMixedModelProblem & other) : TBase(other){
}


TPZMixedModelProblem & TPZMixedModelProblem::operator=(const TPZMixedModelProblem & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    return *this;
}


TPZMixedModelProblem::~TPZMixedModelProblem(){
    
}


void TPZMixedModelProblem::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
//        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}


void TPZMixedModelProblem::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
//        datavec[idata].fNeedsSol = true;
        datavec[idata].fDeformedDirections = true;
    }
}

void TPZMixedModelProblem::Print(std::ostream &out) const{
    TPZMaterial::Print(out);
}

int TPZMixedModelProblem::VariableIndex(const std::string &name) const{
    if(!strcmp("Flux",name.c_str()))            return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    return TPZMaterial::VariableIndex(name);
}


int TPZMixedModelProblem::NSolutionVariables(int var) const{
    if(var == 1) return 3;
    if(var == 2) return 1;
    return TBase::NSolutionVariables(var);
}


void TPZMixedModelProblem::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    TBase::Solution(datavec,var,Solout);
}


void TPZMixedModelProblem::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    // First do the usual Poisson equation contribution
    TBase::Contribute(datavec, weight, ek, ef);
    
    // Now, add the terms relative to u.utest

    int qb = 0;
    int pb = 1;

    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;

    int64_t nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int64_t nphi_p       = phi_ps.Rows();
    int64_t first_q      = 0;
    int64_t first_p      = nphi_q + first_q;
    STATE alpha = 10000.;
    for (int64_t ip = 0; ip < nphi_p; ip++) {
        for (int64_t jp = 0; jp < nphi_p; jp++) {
            ek(ip + first_p, jp + first_p) += alpha*weight * ( - phi_ps(ip,0) ) * phi_ps(jp,0);
        }
    }


}

void TPZMixedModelProblem::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    DebugStop(); // am I being used?
    this->Contribute(datavec, weight, ef);
}


void TPZMixedModelProblem::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}


void TPZMixedModelProblem::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    TBase::ContributeBC(datavec,weight,ek,ef,bc);
}

