#ifndef TPZMixedDarcyFlowOrtotropic_H
#define TPZMixedDarcyFlowOrtotropic_H

#include "DarcyFlow/TPZMixedDarcyFlow.h"

/**
 * @ingroup material
 * @brief This class implements a mixed approximation for the Darcy flow equation with a zero order in the conservation eq. The constitutive law remains the same:
 * div(sigma) + u  = f
 * sigma = -K.grad(u)
 */

class TPZMixedModelProblem : public TPZMixedDarcyFlow {
    using TBase = TPZMixedDarcyFlow;
    
public:
    
    /// Default constructor
    TPZMixedModelProblem();
    
    /// Constructor based on a material id
    TPZMixedModelProblem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TPZMixedModelProblem(const TPZMixedModelProblem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TPZMixedModelProblem &operator=(const TPZMixedModelProblem & other);
    
    /// Default destructor
    ~TPZMixedModelProblem();
    
    /// Set the required data at each integration point
    void FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec)const override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    /// Returns the name of the material
    std::string Name() const override {
        return "TPZMixedModelProblem";
    }
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TPZMixedModelProblem(*this);
    }
    
    /// Set data transfer object
    //    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) const override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) const override;
    
    /// returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var) const override;
    
    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    
    // Contribute Methods being used
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
        
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)override;
    
};

#endif
