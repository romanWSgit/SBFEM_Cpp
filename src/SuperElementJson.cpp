/**
 * \file   SuperElementJson.cpp
 * \author Roman Wallner- Silberhuber
 * \date   21.05.23
 * \brief  This implements the SuperElementJson header file
 */

#include "SuperElementJson.h"

//
// SuperElementJson::SuperElementJson()
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

SuperElementJson::SuperElementJson(const nlohmann::json &json)
{
    try
    {
        m_SE_type = json.at("Type");
        m_SE_sbfemDomainType = json.at("sbfemDomainType");
        m_SE_structName = json.at("struct");
        m_SE_meshPartition = json.at("meshPartition");
        std::string shapeFct = json.at("shapeFct");
        if (shapeFct == "standard shape functions")
        {
            m_SE_shapeFct = ShapeFunctionType::STANDARD;
        }
        else if (shapeFct == "hierarchical shape functions")
        {
            m_SE_shapeFct = ShapeFunctionType::HIERARCHICAL;
        }
        else
        {
            throw std::runtime_error("Invalid shape function: " + shapeFct);
        }
        m_SE_polyOrd = json.at("polyOrd");
        m_SE_geom = json.at("geom");
        m_SE_centre = Eigen::VectorXd::Map(
            json.at("centre").get<std::vector<double>>().data(),
            json.at("centre").size());
        m_SE_ielem = json.at("ielem");
        m_SE_dim = json.at("dim");
        m_SE_nodedim = json.at("nodedim");
        m_SE_evMethod = json.at("evMethod");
        m_SE_blockPartition = json.at("blockPartition");
        m_SE_schurAlgorithm = json.at("schurAlgorithm");
        m_SE_nodes = jsonToEigenMatrix<Eigen::MatrixXd>(json.at("nodes"));
        m_SE_nodesC = jsonToEigenMatrix<Eigen::MatrixXd>(json.at("nodesC"));
        m_SE_elementsM =
            jsonToEigenMatrix<Eigen::MatrixXd>(json.at("elementsM"));
        m_SE_elementsMC =
            jsonToEigenMatrix<Eigen::MatrixXd>(json.at("elementsMC"));
        m_SE_offset = Eigen::VectorXd::Map(
            json.at("offset").get<std::vector<double>>().data(),
            json.at("offset").size());
        m_SE_bcSideFaceFlag = json.at("bcSideFaceFlag");
        m_SE_c = json.at("c");
        m_SE_a = json.at("a");
        m_SE_b = json.at("b");
        m_SE_l = json.at("l");
        m_SE_d = json.at("d");
        m_SE_EZT = jsonToEigenMatrix<Eigen::MatrixXi>(json.at("EZT"));
        m_SE_ltg = jsonToEigenMatrix<Eigen::MatrixXi>(json.at("ltg"));
        m_SE_ltgC = jsonToEigenMatrix<Eigen::MatrixXi>(json.at("ltgC"));
        m_SE_ltgU = jsonToEigenMatrix<Eigen::MatrixXi>(json.at("ltgU"));
        m_SE_ltgCNew = jsonToEigenMatrix<Eigen::MatrixXi>(json.at("ltgCNew"));
        m_SE_ltgUNew = jsonToEigenMatrix<Eigen::MatrixXi>(json.at("ltgUNew"));
        m_SE_dimC = json.at("dimC");
        m_SE_dimU = json.at("dimU");
    }
    catch (nlohmann::json::out_of_range &e)
    {
        throw std::runtime_error("JSON key not found: " +
                                 std::string(e.what()));
    }
    catch (std::exception &e)
    {
        // Catch all other exceptions
        throw std::runtime_error("An error occurred during parsing: " +
                                 std::string(e.what()));
    }
}

// Getters
std::string SuperElementJson::getMSEType() const
{
    return m_SE_type;
}
std::string SuperElementJson::getMSESbfemDomainType() const
{
    return m_SE_sbfemDomainType;
}
std::string SuperElementJson::getMSEStructName() const
{
    return m_SE_structName;
}
std::string SuperElementJson::getMSEMeshPartition() const
{
    return m_SE_meshPartition;
}
ShapeFunctionType SuperElementJson::getMSEShapeFct() const
{
    return m_SE_shapeFct;
}
int SuperElementJson::getMSEPolyOrd() const
{
    return m_SE_polyOrd;
}
std::string SuperElementJson::getMSEGeom() const
{
    return m_SE_geom;
}
Eigen::VectorXd SuperElementJson::getMSECentre() const
{
    return m_SE_centre;
}
int SuperElementJson::getMSEIelem() const
{
    return m_SE_ielem;
}
int SuperElementJson::getMSEDim() const
{
    return m_SE_dim;
}
int SuperElementJson::getMSENodedim() const
{
    return m_SE_nodedim;
}
std::string SuperElementJson::getMSEEvMethod() const
{
    return m_SE_evMethod;
}
std::string SuperElementJson::getMSEBlockPartition() const
{
    return m_SE_blockPartition;
}
std::string SuperElementJson::getMSESchurAlgorithm() const
{
    return m_SE_schurAlgorithm;
}
Eigen::MatrixXd SuperElementJson::getMSENodes() const
{
    return m_SE_nodes;
}
Eigen::MatrixXd SuperElementJson::getMSENodesC() const
{
    return m_SE_nodesC;
}
Eigen::MatrixXd SuperElementJson::getMSEElementsM() const
{
    return m_SE_elementsM;
}
Eigen::MatrixXd SuperElementJson::getMSEElementsMC() const
{
    return m_SE_elementsMC;
}
Eigen::VectorXd SuperElementJson::getMSEOffset() const
{
    return m_SE_offset;
}
bool SuperElementJson::getMSEBcSideFaceFlag() const
{
    return m_SE_bcSideFaceFlag;
}
int SuperElementJson::getMSEC() const
{
    return m_SE_c;
}
int SuperElementJson::getMSEA() const
{
    return m_SE_a;
}
int SuperElementJson::getMSEB() const
{
    return m_SE_b;
}
int SuperElementJson::getMSEL() const
{
    return m_SE_l;
}
int SuperElementJson::getMSED() const
{
    return m_SE_d;
}
Eigen::MatrixXi SuperElementJson::getMSEEZT() const
{
    return m_SE_EZT;
}
Eigen::MatrixXi SuperElementJson::getMSELtg() const
{
    return m_SE_ltg;
}
Eigen::MatrixXi SuperElementJson::getMSELtgC() const
{
    return m_SE_ltgC;
}
Eigen::MatrixXi SuperElementJson::getMSELtgU() const
{
    return m_SE_ltgU;
}
Eigen::MatrixXi SuperElementJson::getMSELtgCNew() const
{
    return m_SE_ltgCNew;
}
Eigen::MatrixXi SuperElementJson::getMSELtgUNew() const
{
    return m_SE_ltgUNew;
}
int SuperElementJson::getMSEDimC() const
{
    return m_SE_dimC;
}
int SuperElementJson::getMSEDimU() const
{
    return m_SE_dimU;
}

// Setters
void SuperElementJson::setMSEType(const std::string &type)
{
    this->m_SE_type = type;
}
void SuperElementJson::setMSESbfemDomainType(const std::string &sbfemDomainType)
{
    this->m_SE_sbfemDomainType = sbfemDomainType;
}
void SuperElementJson::setMSEStructName(const std::string &structName)
{
    this->m_SE_structName = structName;
}
void SuperElementJson::setMSEMeshPartition(const std::string &meshPartition)
{
    this->m_SE_meshPartition = meshPartition;
}
void SuperElementJson::setMSEShapeFct(const enum ShapeFunctionType &shapeFct)
{
    this->m_SE_shapeFct = shapeFct;
}
void SuperElementJson::setMSEPolyOrd(int polyOrd)
{
    this->m_SE_polyOrd = polyOrd;
}
void SuperElementJson::setMSEGeom(const std::string &geom)
{
    this->m_SE_geom = geom;
}
void SuperElementJson::setMSECentre(const Eigen::VectorXd &centre)
{
    this->m_SE_centre = centre;
}
void SuperElementJson::setMSEIelem(int ielem)
{
    this->m_SE_ielem = ielem;
}
void SuperElementJson::setMSEDim(int dim)
{
    this->m_SE_dim = dim;
}
void SuperElementJson::setMSENodedim(int nodedim)
{
    this->m_SE_nodedim = nodedim;
}
void SuperElementJson::setMSEEvMethod(const std::string &evMethod)
{
    this->m_SE_evMethod = evMethod;
}
void SuperElementJson::setMSEBlockPartition(const std::string &blockPartition)
{
    this->m_SE_blockPartition = blockPartition;
}
void SuperElementJson::setMSESchurAlgorithm(const std::string &schurAlgorithm)
{
    this->m_SE_schurAlgorithm = schurAlgorithm;
}
void SuperElementJson::setMSENodes(const Eigen::MatrixXd &nodes)
{
    this->m_SE_nodes = nodes;
}
void SuperElementJson::setMSENodesC(const Eigen::MatrixXd &nodesC)
{
    this->m_SE_nodesC = nodesC;
}
void SuperElementJson::setMSEElementsM(const Eigen::MatrixXd &elementsM)
{
    this->m_SE_elementsM = elementsM;
}
void SuperElementJson::setMSEElementsMC(const Eigen::MatrixXd &elementsMC)
{
    this->m_SE_elementsMC = elementsMC;
}
void SuperElementJson::setMSEOffset(const Eigen::VectorXd &offset)
{
    this->m_SE_offset = offset;
}
void SuperElementJson::setMSEBcSideFaceFlag(bool bcSideFaceFlag)
{
    this->m_SE_bcSideFaceFlag = bcSideFaceFlag;
}
void SuperElementJson::setMSEC(int c)
{
    this->m_SE_c = c;
}
void SuperElementJson::setMSEA(int a)
{
    this->m_SE_a = a;
}
void SuperElementJson::setMSEB(int b)
{
    this->m_SE_b = b;
}
void SuperElementJson::setMSEL(int l)
{
    this->m_SE_l = l;
}
void SuperElementJson::setMSED(int d)
{
    this->m_SE_d = d;
}
void SuperElementJson::setMSEEZT(const Eigen::MatrixXi &EZT)
{
    this->m_SE_EZT = EZT;
}
void SuperElementJson::setMSELtg(const Eigen::MatrixXi &ltg)
{
    this->m_SE_ltg = ltg;
}
void SuperElementJson::setMSELtgC(const Eigen::MatrixXi &ltgC)
{
    this->m_SE_ltgC = ltgC;
}
void SuperElementJson::setMSELtgU(const Eigen::MatrixXi &ltgU)
{
    this->m_SE_ltgU = ltgU;
}
void SuperElementJson::setMSELtgCNew(const Eigen::MatrixXi &ltgCNew)
{
    this->m_SE_ltgCNew = ltgCNew;
}
void SuperElementJson::setMSELtgUNew(const Eigen::MatrixXi &ltgUNew)
{
    this->m_SE_ltgUNew = ltgUNew;
}
void SuperElementJson::setMSEDimC(int dimC)
{
    this->m_SE_dimC = dimC;
}
void SuperElementJson::setMSEDimU(int dimU)
{
    this->m_SE_dimU = dimU;
}

std::ostream &operator<<(std::ostream &os, const SuperElementJson &sE)
{
    os << "SuperElementJson: \n"
       << "polygonOrder: " << sE.m_SE_polyOrd << '\n'
       << "number of elements: " << sE.m_SE_ielem << '\n'
       << "\n";
    return os;
}
