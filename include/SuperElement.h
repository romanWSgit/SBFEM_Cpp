/**
 * \file   SuperElement.h
 * \author Roman Wallner- Silberhuber
 * \date   21.05.23
 * \brief  This head contains the SuperElement class
 */

#ifndef SBFEM_SUPERELEMENT_H
#define SBFEM_SUPERELEMENT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tuple>
#include <nlohmann/json.hpp>
#include "helper_functions.h"
#include <iostream>

/**
 * @class SuperElement
 * @brief A class representing a SuperElement with various properties.
 */
class SuperElement
{
private:

    // naming convention: m_SE_SE_ for member Variables of the SuperElement class
    std::string m_SE_type;
    std::string m_SE_sbfemDomainType;
    std::string m_SE_structName;
    std::string m_SE_meshPartition;
    std::string m_SE_shapeFct;
    int m_SE_polyOrd;
    std::string m_SE_geom;
    Eigen::VectorXd m_SE_centre;
    int m_SE_ielem;
    int m_SE_dim;
    int m_SE_nodedim;
    std::string m_SE_evMethod;
    std::string m_SE_blockPartition;
    std::string m_SE_schurAlgorithm;
    Eigen::MatrixXd m_SE_nodes;
    Eigen::MatrixXd m_SE_nodesC;
    Eigen::MatrixXd m_SE_elementsM;
    Eigen::MatrixXd m_SE_elementsMC;
    Eigen::VectorXd m_SE_offset;
    bool m_SE_bcSideFaceFlag;
    int m_SE_c, m_SE_a, m_SE_b, m_SE_l, m_SE_d;
    Eigen::MatrixXi m_SE_EZT;
    Eigen::MatrixXi m_SE_ltg;
    Eigen::MatrixXi m_SE_ltgC;
    Eigen::MatrixXi m_SE_ltgU;
    Eigen::MatrixXi m_SE_ltgCNew;
    Eigen::MatrixXi m_SE_ltgUNew;
    int m_SE_dimC, m_SE_dimU;







public:
    // Constructor to initialize from JSON
    explicit SuperElement(const nlohmann::json& json);

    // Getters
    [[nodiscard]] std::string getMSEType() const;
    [[nodiscard]] std::string getMSESbfemDomainType() const;
    [[nodiscard]] std::string getMSEStructName() const;
    [[nodiscard]] std::string getMSEMeshPartition() const;
    [[nodiscard]] std::string getMSEShapeFct() const;
    [[nodiscard]] int getMSEPolyOrd() const;
    [[nodiscard]] std::string getMSEGeom() const;
    [[nodiscard]] Eigen::VectorXd getMSECentre() const;
    [[nodiscard]] int getMSEIelem() const;
    [[nodiscard]] int getMSEDim() const;
    [[nodiscard]] int getMSENodedim() const;
    [[nodiscard]] std::string getMSEEvMethod() const;
    [[nodiscard]] std::string getMSEBlockPartition() const;
    [[nodiscard]] std::string getMSESchurAlgorithm() const;
    [[nodiscard]] Eigen::MatrixXd getMSENodes() const;
    [[nodiscard]] Eigen::MatrixXd getMSENodesC() const;
    [[nodiscard]] Eigen::MatrixXd getMSEElementsM() const;
    [[nodiscard]] Eigen::MatrixXd getMSEElementsMC() const;
    [[nodiscard]] Eigen::VectorXd getMSEOffset() const;
    [[nodiscard]] bool getMSEBcSideFaceFlag() const;
    [[nodiscard]] int getMSEC() const;
    [[nodiscard]] int getMSEA() const;
    [[nodiscard]] int getMSEB() const;
    [[nodiscard]] int getMSEL() const;
    [[nodiscard]] int getMSED() const;
    [[nodiscard]] Eigen::MatrixXi getMSEEZT() const;
    [[nodiscard]] Eigen::MatrixXi getMSELtg() const;
    [[nodiscard]] Eigen::MatrixXi getMSELtgC() const;
    [[nodiscard]] Eigen::MatrixXi getMSELtgU() const;
    [[nodiscard]] Eigen::MatrixXi getMSELtgCNew() const;
    [[nodiscard]] Eigen::MatrixXi getMSELtgUNew() const;
    [[nodiscard]] int getMSEDimC() const;
    [[nodiscard]] int getMSEDimU() const;

    // Setters
    void setMSEType (const std::string& type);
    void setMSESbfemDomainType (const std::string& sbfemDomainType);
    void setMSEStructName (const std::string& structName);
    void setMSEMeshPartition (const std::string& meshPartition);
    void setMSEShapeFct (const std::string& shapeFct);
    void setMSEPolyOrd (int polyOrd);
    void setMSEGeom (const std::string& geom);
    void setMSECentre (const Eigen::VectorXd& centre);
    void setMSEIelem (int ielem);
    void setMSEDim (int dim);
    void setMSENodedim (int nodedim);
    void setMSEEvMethod (const std::string& evMethod);
    void setMSEBlockPartition (const std::string& blockPartition);
    void setMSESchurAlgorithm (const std::string& schurAlgorithm);
    void setMSENodes (const Eigen::MatrixXd& nodes);
    void setMSENodesC (const Eigen::MatrixXd& nodesC);
    void setMSEElementsM (const Eigen::MatrixXd& elementsM);
    void setMSEElementsMC (const Eigen::MatrixXd& elementsMC);
    void setMSEOffset (const Eigen::VectorXd& offset);
    void setMSEBcSideFaceFlag (bool bcSideFaceFlag);
    void setMSEC (int c);
    void setMSEA (int a);
    void setMSEB (int b);
    void setMSEL (int l);
    void setMSED (int d);
    void setMSEEZT (const Eigen::MatrixXi& EZT);
    void setMSELtg (const Eigen::MatrixXi& ltg);
    void setMSELtgC (const Eigen::MatrixXi& ltgC);
    void setMSELtgU (const Eigen::MatrixXi& ltgU);
    void setMSELtgCNew (const Eigen::MatrixXi& ltgCNew);
    void setMSELtgUNew (const Eigen::MatrixXi& ltgUNew);
    void setMSEDimC (int dimC);
    void setMSEDimU (int dimU);


    /**
     * @brief Overloaded output stream operator for SuperElement class.
     *
     * This operator allows the SuperElement objects to be easily printed using the
     * standard output stream.
     *
     * @param os The output stream.
     * @param sEt The object of SuperElement class to be printed.
     * @return std::ostream& The modified output stream with the SuperElement printed.
     */
    friend std::ostream& operator << (std::ostream& os, const SuperElement& sEt);
};





#endif //SBFEM_SUPERELEMENT_H
