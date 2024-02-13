/**
 * \file   SuperElement.cpp
 * \author Roman Wallner- Silberhuber
 * \date   21.05.23
 * \brief  This implements the SuperElement header file
 */

#include "SuperElement.h"

//
//SuperElement::SuperElement()
//{
//    polygon_order = 1;
//    number_of_elements = 10;
//    nodes = Eigen::Matrix<double, 2, 20>::Zero();
//    nodes_C = Eigen::Matrix<double, 2, 20>::Zero();
//    elements_M = Eigen::Matrix<double, 4, 10>::Zero();
//    elements_M_C = Eigen::Matrix<double, 4, 10>::Zero();
//    ltg = Eigen::Matrix<int, 4, 10>::Zero();
//
//}
//
//

SuperElement::SuperElement(const nlohmann::json& json)
{
    m_SE_type = json["Type"];
    m_SE_sbfemDomainType = json["sbfemDomainType"];
    m_SE_structName = json["struct"];
    m_SE_meshPartition = json["meshPartition"];
    m_SE_shapeFct = json["shapeFct"];
    m_SE_polyOrd = json["polyOrd"];
    m_SE_geom = json["geom"];
    m_SE_centre = Eigen::VectorXd::Map(json["centre"].get<std::vector<double>>().data(), json["centre"].size());
    m_SE_ielem = json["ielem"];
    m_SE_dim = json["dim"];
    m_SE_nodedim = json["nodedim"];
    m_SE_evMethod = json["evMethod"];
    m_SE_blockPartition = json["blockPartition"];
    m_SE_schurAlgorithm = json["schurAlgorithm"];
    m_SE_nodes = jsonToEigenMatrix<Eigen::MatrixXd>(json["nodes"]);
    m_SE_nodesC = jsonToEigenMatrix<Eigen::MatrixXd>(json["nodesC"]);
    m_SE_elementsM = jsonToEigenMatrix<Eigen::MatrixXd>(json["elementsM"]);
    m_SE_elementsMC = jsonToEigenMatrix<Eigen::MatrixXd>(json["elementsMC"]);
    m_SE_offset = Eigen::VectorXd::Map(json["offset"].get<std::vector<double>>().data(), json["offset"].size());
    m_SE_bcSideFaceFlag = json["bcSideFaceFlag"];
    m_SE_c = json["c"];
    m_SE_a = json["a"];
    m_SE_b = json["b"];
    m_SE_l = json["l"];
    m_SE_d = json["d"];
    m_SE_EZT =jsonToEigenMatrix<Eigen::MatrixXi>(json["EZT"]);
    m_SE_ltg = jsonToEigenMatrix<Eigen::MatrixXi>(json["ltg"]);
    m_SE_ltgC = jsonToEigenMatrix<Eigen::MatrixXi>(json["ltgC"]);
    m_SE_ltgU = jsonToEigenMatrix<Eigen::MatrixXi>(json["ltgU"]);
    m_SE_ltgCNew = jsonToEigenMatrix<Eigen::MatrixXi>(json["ltgCNew"]);
    m_SE_ltgUNew = jsonToEigenMatrix<Eigen::MatrixXi>(json["ltgUNew"]);
    m_SE_dimC = json["dimC"];
    m_SE_dimU = json["dimU"];
}

    // Getters
    std::string SuperElement::getMSEType() const { return m_SE_type; }
    std::string SuperElement::getMSESbfemDomainType() const { return m_SE_sbfemDomainType; }
    std::string SuperElement::getMSEStructName() const { return m_SE_structName; }
    std::string SuperElement::getMSEMeshPartition() const { return m_SE_meshPartition; }
    std::string SuperElement::getMSEShapeFct() const { return m_SE_shapeFct; }
    int SuperElement::getMSEPolyOrd() const { return m_SE_polyOrd; }
    std::string SuperElement::getMSEGeom() const { return m_SE_geom; }
    Eigen::VectorXd SuperElement::getMSECentre() const { return m_SE_centre; }
    int SuperElement::getMSEIelem() const { return m_SE_ielem; }
    int SuperElement::getMSEDim() const { return m_SE_dim; }
    int SuperElement::getMSENodedim() const { return m_SE_nodedim; }
    std::string SuperElement::getMSEEvMethod() const { return m_SE_evMethod; }
    std::string SuperElement::getMSEBlockPartition() const { return m_SE_blockPartition; }
    std::string SuperElement::getMSESchurAlgorithm() const { return m_SE_schurAlgorithm; }
    Eigen::MatrixXd SuperElement::getMSENodes() const { return m_SE_nodes; }
    Eigen::MatrixXd SuperElement::getMSENodesC() const { return m_SE_nodesC; }
    Eigen::MatrixXd SuperElement::getMSEElementsM() const { return m_SE_elementsM; }
    Eigen::MatrixXd SuperElement::getMSEElementsMC() const { return m_SE_elementsMC; }
    Eigen::VectorXd SuperElement::getMSEOffset() const { return m_SE_offset; }
    bool SuperElement::getMSEBcSideFaceFlag() const { return m_SE_bcSideFaceFlag; }
    int SuperElement::getMSEC() const { return m_SE_c; }
    int SuperElement::getMSEA() const { return m_SE_a; }
    int SuperElement::getMSEB() const { return m_SE_b; }
    int SuperElement::getMSEL() const { return m_SE_l; }
    int SuperElement::getMSED() const { return m_SE_d; }
    Eigen::MatrixXi SuperElement::getMSEEZT() const { return m_SE_EZT; }
    Eigen::MatrixXi SuperElement::getMSELtg() const { return m_SE_ltg; }
    Eigen::MatrixXi SuperElement::getMSELtgC() const { return m_SE_ltgC; }
    Eigen::MatrixXi SuperElement::getMSELtgU() const { return m_SE_ltgU; }
    Eigen::MatrixXi SuperElement::getMSELtgCNew() const { return m_SE_ltgCNew; }
    Eigen::MatrixXi SuperElement::getMSELtgUNew() const { return m_SE_ltgUNew; }
    int SuperElement::getMSEDimC() const { return m_SE_dimC; }
    int SuperElement::getMSEDimU() const { return m_SE_dimU; }

    // Setters
    void SuperElement::setMSEType(const std::string& type) { this->m_SE_type = type; }
    void SuperElement::setMSESbfemDomainType(const std::string& sbfemDomainType) { this->m_SE_sbfemDomainType = sbfemDomainType; }
    void SuperElement::setMSEStructName(const std::string& structName) { this->m_SE_structName = structName; }
    void SuperElement::setMSEMeshPartition(const std::string& meshPartition) { this->m_SE_meshPartition = meshPartition; }
    void SuperElement::setMSEShapeFct(const std::string& shapeFct) { this->m_SE_shapeFct = shapeFct; }
    void SuperElement::setMSEPolyOrd(int polyOrd) { this->m_SE_polyOrd = polyOrd; }
    void SuperElement::setMSEGeom(const std::string& geom) { this->m_SE_geom = geom; }
    void SuperElement::setMSECentre(const Eigen::VectorXd& centre) { this->m_SE_centre = centre; }
    void SuperElement::setMSEIelem(int ielem) { this->m_SE_ielem = ielem; }
    void SuperElement::setMSEDim(int dim) { this->m_SE_dim = dim; }
    void SuperElement::setMSENodedim(int nodedim) { this->m_SE_nodedim = nodedim; }
    void SuperElement::setMSEEvMethod(const std::string& evMethod) { this->m_SE_evMethod = evMethod; }
    void SuperElement::setMSEBlockPartition(const std::string& blockPartition) { this->m_SE_blockPartition = blockPartition; }
    void SuperElement::setMSESchurAlgorithm(const std::string& schurAlgorithm) { this->m_SE_schurAlgorithm = schurAlgorithm; }
    void SuperElement::setMSENodes(const Eigen::MatrixXd& nodes) { this->m_SE_nodes = nodes; }
    void SuperElement::setMSENodesC(const Eigen::MatrixXd& nodesC) { this->m_SE_nodesC = nodesC; }
    void SuperElement::setMSEElementsM(const Eigen::MatrixXd& elementsM) { this->m_SE_elementsM = elementsM; }
    void SuperElement::setMSEElementsMC(const Eigen::MatrixXd& elementsMC) { this->m_SE_elementsMC = elementsMC; }
    void SuperElement::setMSEOffset(const Eigen::VectorXd& offset) { this->m_SE_offset = offset; }
    void SuperElement::setMSEBcSideFaceFlag(bool bcSideFaceFlag) { this->m_SE_bcSideFaceFlag = bcSideFaceFlag; }
    void SuperElement::setMSEC(int c) { this->m_SE_c = c; }
    void SuperElement::setMSEA(int a) { this->m_SE_a = a; }
    void SuperElement::setMSEB(int b) { this->m_SE_b = b; }
    void SuperElement::setMSEL(int l) { this->m_SE_l = l; }
    void SuperElement::setMSED(int d) { this->m_SE_d = d; }
    void SuperElement::setMSEEZT(const Eigen::MatrixXi& EZT) { this->m_SE_EZT = EZT; }
    void SuperElement::setMSELtg(const Eigen::MatrixXi& ltg) { this->m_SE_ltg = ltg; }
    void SuperElement::setMSELtgC(const Eigen::MatrixXi& ltgC) { this->m_SE_ltgC = ltgC; }
    void SuperElement::setMSELtgU(const Eigen::MatrixXi& ltgU) { this->m_SE_ltgU = ltgU; }
    void SuperElement::setMSELtgCNew(const Eigen::MatrixXi& ltgCNew) { this->m_SE_ltgCNew = ltgCNew; }
    void SuperElement::setMSELtgUNew(const Eigen::MatrixXi& ltgUNew) { this->m_SE_ltgUNew = ltgUNew; }
    void SuperElement::setMSEDimC(int dimC) { this->m_SE_dimC = dimC; }
    void SuperElement::setMSEDimU(int dimU) { this->m_SE_dimU = dimU; }

std::ostream& operator << (std::ostream& os, const SuperElement& sE)
{
    os << "SuperElement: \n" <<
        "polygonOrder: " << sE.m_SE_polyOrd << '\n' <<
        "number of elements: " << sE.m_SE_ielem << '\n'
        << "\n";
    return os;
}


