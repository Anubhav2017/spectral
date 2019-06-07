#include "parameterization.h"

#include "CellMesh/celledge.h"

#include <cmath>

#include <QSet>
#include <QImage>
#include <QPainter>
#include <QFile>
#include <QTextStream>

#include <QDebug>

Parameterization::Parameterization()
{
}

bool Parameterization::contains(Vertex *v)
{
    return m_vertexLocalIdMap.contains(v);
}

void Parameterization::updateParameter(Vertex *v, double newU, double newV)
{
    const int index = m_vertexLocalIdMap[v];
    m_parameterValuesU[index] = newU;
    m_parameterValuesV[index] = newV;
}

QVector2D Parameterization::getParameter(Vertex *v)
{
    const int index = m_vertexLocalIdMap[v];
    return QVector2D(m_parameterValuesU[index], m_parameterValuesV[index]);
}

QVector<double> Parameterization::getParameterizationOfBoundaryCurve(int id)
{
    if (m_borderType != Parameterization::BT_equidistantUnitSquare && m_borderType != Parameterization::BT_lengthWeightedUnitSquare) {
        qDebug() << "ERROR: Tried to access border curve parameterization for a non-quad cell! ID:" << m_cellFace->getId();
        return QVector<double>();
    }

    CellEdge *ce = (CellEdge *) m_cellFace->getEdge(id);
    const int numPLV = ce->getPolylineVertices()->size();

    QVector<double> ceParam(numPLV);
    for (int i = 0; i < numPLV; i++) {
        QVector2D plvParam = this->getParameter(ce->getPolylineVertices()->at(i));
        if (m_cellFace->edgeIsInverted(id)) {
            switch (id) {
            case 0:
                ceParam[i] = m_parameterInterval.second - plvParam.x();
                break;
            case 1:
                ceParam[i] = m_parameterInterval.second - plvParam.y();
                break;
            case 2:
                ceParam[i] = plvParam.x();
                break;
            case 3:
                ceParam[i] = plvParam.y();
                break;
            default:
                break;
            }
        } else {
            switch (id) {
            case 0:
                ceParam[i] = plvParam.x();
                break;
            case 1:
                ceParam[i] = plvParam.y();
                break;
            case 2:
                ceParam[i] = m_parameterInterval.second - plvParam.x();
                break;
            case 3:
                ceParam[i] = m_parameterInterval.second - plvParam.y();
                break;
            default:
                break;
            }
        }
    }
    return ceParam;
}

QVector2D Parameterization::getPolarCoordinates(Vertex *v)
{
    const QVector2D center(m_parameterInterval.first/2, m_parameterInterval.second/2);
    const QVector2D paramValue = this->getParameter(v) - center;
    const double x = paramValue.x();
    const double y = paramValue.y();
    const double length = paramValue.length();
    double angle = 0;
    if (y > 0) {
        double ratio = x/length;
        if (ratio > 1)
            ratio = 1;
        else if (ratio < -1)
            ratio = -1;
        angle = acos(ratio);
    } else if (y == 0 && x < 0) {
        angle = M_PI;
    } else if (y < 0) {
        double ratio = x/length;
        if (ratio > 1)
            ratio = 1;
        else if (ratio < -1)
            ratio = -1;
        angle = (2 * M_PI) - acos(ratio);
    }

    return QVector2D(angle, length);
}

double Parameterization::getPolarAngle(Vertex *v)
{
    const double centerParameterValue = m_parameterInterval.first + (m_parameterInterval.second - m_parameterInterval.first)/2;
    const QVector2D center(centerParameterValue, centerParameterValue);
    const QVector2D paramValue = this->getParameter(v) - center;
    const double x = paramValue.x();
    const double y = paramValue.y();
    const double length = paramValue.length();
    double angle = 0;
    if (y > 0) {
        double ratio = x/length;
        if (ratio > 1)
            ratio = 1;
        else if (ratio < -1)
            ratio = -1;
        angle = acos(ratio);
    } else if (y == 0 && x < 0) {
        angle = M_PI;
    } else if (y < 0) {
        double ratio = x/length;
        if (ratio > 1)
            ratio = 1;
        else if (ratio < -1)
            ratio = -1;
        angle = (2 * M_PI) - acos(ratio);
    }
    return angle;
}

void Parameterization::saveAsImage(QString filename, int size, int symbolRadius1, int symbolRadius2)
{
    QImage img(size, size, QImage::Format_ARGB32);
    img.fill(0);
    QPainter p(&img);
    p.setBrush(QBrush(Qt::white));
    p.drawRect(-1, -1, size + 2, size + 2);

    const double a = m_parameterInterval.first;
    const QVector2D aa(a, a);
    const double b = m_parameterInterval.second;
    const double scale = size * 0.9 /(b - a);
    const QVector2D origin(size * 0.05, size * 0.05);


    //Draw endpoints of cell edges
    p.setPen(Qt::black);
    for (int i = 0; i < m_cellFace->getNumberOfEdges(); i++) {
        CellEdge *ce = (CellEdge *) m_cellFace->getEdge(i);
        QVector2D rectCenter = origin + scale * (this->getParameter(ce->getPolylineVertices()->first()) - aa);
        p.drawRect(QRect(QPoint(rectCenter.x() - symbolRadius2, rectCenter.y() - symbolRadius2),
                         QPoint(rectCenter.x() + symbolRadius2, rectCenter.y() + symbolRadius2)));
        rectCenter = origin + scale * (this->getParameter(ce->getPolylineVertices()->last()) - aa);
        p.drawRect(QRect(QPoint(rectCenter.x() - symbolRadius2, rectCenter.y() - symbolRadius2),
                         QPoint(rectCenter.x() + symbolRadius2, rectCenter.y() + symbolRadius2)));
    }

    //Draw Vertices (inner and border)
    p.setPen(Qt::darkGray);
    for (int i = 0; i < m_numberOfVertices; i++)
        p.drawEllipse((origin + scale * (QVector2D(m_parameterValuesU[i], m_parameterValuesV[i]) - aa)).toPoint(), symbolRadius1, symbolRadius1);

    //Draw inner edges
    p.setPen(Qt::lightGray);
    for (int i = 0; i < m_cellFace->getMeshEdges()->size(); i++) {
        Edge *e = m_cellFace->getMeshEdges()->at(i);
        p.drawLine((origin + scale * this->getParameter(e->getVertex(0)) - aa).toPoint(), (origin + scale * this->getParameter(e->getVertex(1)) - aa).toPoint());
    }

    //Draw edges connecting border vertices to inner vertices
    for (int i = 0; i < m_cellFace->getNumberOfEdges(); i++) {
        CellEdge *ce = (CellEdge *) m_cellFace->getEdge(i);
        foreach (PolylineVertex *plv, *ce->getPolylineVertices()) {
            QVector2D pos = this->getParameter(plv) - aa;
            if (plv->isMeshVertex()) {
                foreach (int eId, *plv->getMeshVertex()->getIncidentEdgeIds()) {
                    Edge *e = m_originalMesh->getEdge(eId);
                    Vertex *v = e->getConnectedVertex(plv->getMeshVertex());
                    if (m_vertexLocalIdMap.contains(v))
                        p.drawLine((origin + scale * pos).toPoint(), (origin + scale * (this->getParameter(v) - aa)).toPoint());
                }
            } else if (plv->isCrossingVertex()) {
                for (int j = 0; j < 2; j++) {
                    Vertex *v = plv->getCrossingEdge()->getVertex(j);
                    if (m_vertexLocalIdMap.contains(v))
                        p.drawLine((origin + scale * pos).toPoint(), (origin + scale * (this->getParameter(v) - aa)).toPoint());
                }
            }
        }
    }

    //Draw edges that connect the border vertices
    p.setPen(Qt::black);
    for (int i = 0; i < m_numberOfBorderVertices; i++) {
        const int nextId = (i + 1) % m_numberOfBorderVertices;
        const QVector2D point0(m_parameterValuesU[m_numberOfInnerVertices + i], m_parameterValuesV[m_numberOfInnerVertices + i]);
        const QVector2D point1(m_parameterValuesU[m_numberOfInnerVertices + nextId], m_parameterValuesV[m_numberOfInnerVertices + nextId]);
        p.drawLine((origin + scale * (point0 - aa)).toPoint(),
                   (origin + scale * (point1 - aa)).toPoint());
    }

    img.save(filename);
}

void Parameterization::addCopy(Vertex *original, Vertex *copy)
{
    m_vertexLocalIdMap[copy] = m_vertexLocalIdMap[original];
}

Parameterization::BorderType Parameterization::getBorderType()
{
    return m_borderType;
}

void Parameterization::saveToFile(QString filename, QString optionalHeader)
{
    const int numVertices = m_cellFace->getMeshVertices()->size();
    const int numberOfBorders = m_cellFace->getNumberOfEdges();

    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream out(&file);

    out << m_info;
    if (!optionalHeader.isEmpty()) {
        out << optionalHeader;
        if (!optionalHeader.endsWith(QString("\n")))
            out << "\n";
    }
    out << "===\n";
    out << QString::number(m_borderType) << "\n";
    out << QString::number(m_parameterInterval.first) << " " << QString::number(m_parameterInterval.second) << "\n";
    out << "===\n";

    for (int i = 0; i < numVertices; i++) {
        QVector2D param = this->getParameter(m_cellFace->getMeshVertices()->at(i));
        out << QString::number(m_cellFace->getMeshVertices()->at(i)->getId()) << " " << QString::number(param.x()) << " " << QString::number(param.y()) << "\n";
    }

    out << "===\n";
    for (int iBorder = 0; iBorder < numberOfBorders; iBorder++) {
        CellEdge *ce = (CellEdge *) m_cellFace->getEdge(iBorder);
        const int numberOfPLV = ce->getPolylineVertices()->size();
        for (int i = 0; i < numberOfPLV; i++) {
            QVector2D param = this->getParameter(ce->getPolylineVertices()->at(i));
            out << QString::number(i) << " " << QString::number(param.x()) << " " << QString::number(param.y()) << "\n";
        }
        out << "---\n";
    }
    out << "===\n";
    out << "\n" << "\n";

    file.close();
}

Parameterization *Parameterization::loadFromFile(QString filename, CellFace *cf, Mesh *originalMesh)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly))
        return 0;
    QTextStream in(&file);

    Parameterization *param = new Parameterization();
    param->m_cellFace = cf;
    param->m_originalMesh = originalMesh;

    QString info;
    QString line = in.readLine();
    //skip until first seperator
    while (line != "===") {
        info.append(line).append("\n");
        line = in.readLine();
    }
    param->m_info = info;
    param->fillVertexLocalIdMap();

    //Read mesh dimensions
    line = in.readLine();   //This line contains meta information
    QStringList split = line.split(" ");
    param->m_borderType = static_cast<Parameterization::BorderType>(split[0].toInt());

    line = in.readLine();   //next separator (old version) or parameter interval
    if (line != "===") {
        split = line.split(" ");
        const double paramA = split[0].toDouble();
        const double paramB = split[1].toDouble();
        param->m_parameterInterval = QPair<double, double>(paramA, paramB);
        line = in.readLine();   //next separator
    } else {
        param->m_parameterInterval = QPair<double, double>(0, 1);   //default value for old files
    }

    //read mesh vertices
    line = in.readLine();
    while (line != "===") {
        split = line.split(" ");
        const int id = split[0].toInt();
        param->updateParameter(originalMesh->getVertex(id), split[1].toDouble(), split[2].toDouble());
        line = in.readLine();
    }

    //read polyline vertices
    line = in.readLine();
    int iBorder = 0;
    param->m_numberOfBorderVertices = 0;
    while (line != "===") {
        CellEdge *ce = (CellEdge *) cf->getEdge(iBorder);
        param->m_numberOfBorderVertices += ce->getPolylineVertices()->size();
        while (line != "---") {
            split = line.split(" ");
            const int i = split[0].toInt();
            param->updateParameter(ce->getPolylineVertices()->at(i), split[1].toDouble(), split[2].toDouble());
            line = in.readLine();
        }
        line = in.readLine();

        iBorder++;
    }

    file.close();

    param->m_numberOfInnerVertices = cf->getMeshVertices()->size();
    param->m_numberOfVertices = param->m_numberOfBorderVertices + param->m_numberOfInnerVertices;

    return param;
}

QPair<double, double> Parameterization::getParameterInterval()
{
    return m_parameterInterval;
}

QString Parameterization::getInfo()
{
    return m_info;
}

void Parameterization::addInfo(QString info)
{
    m_info.append(info);
}

double Parameterization::limitToParameterInterval()
{
    return this->limitTo(m_parameterInterval.first, m_parameterInterval.second, m_parameterInterval.first, m_parameterInterval.second);
}

double Parameterization::limitTo(const double lowerBoundU, const double upperBoundU, const double lowerBoundV, const double upperBoundV)
{
    const int numberOfParamValues = m_parameterValuesU.size();

    double error = 0;
    for (int i = 0; i < numberOfParamValues; i++) {
        double paramU = m_parameterValuesU[i];
        double paramV = m_parameterValuesV[i];
        const double origU = paramU;    //copy to compare after modifications are made
        const double origV = paramV;

        if (paramU < lowerBoundU)
            paramU = lowerBoundU;
        if (paramU > upperBoundU)
            paramU = upperBoundU;
        if (paramV < lowerBoundV)
            paramV = lowerBoundV;
        if (paramV > upperBoundV)
            paramV = upperBoundV;

        error += sqrt((origU - paramU)*(origU - paramU) + (origV - paramV) * (origV - paramV));
        if (paramU != origU)
            m_parameterValuesU[i] = paramU;
        if (paramV != origV)
            m_parameterValuesV[i] = paramV;
    }

    return error;
}

Parameterization *Parameterization::harmonicMap(CellFace *cf, Mesh *originalMesh, Parameterization::WeightType wt, Parameterization::BorderType bt, QPair<double, double> parameterInterval, const int expectedEdgesPerVertex)
{
    //Most computation time is needed for:
    //Setting up the weight matrices
    //Solving the sparse linear system

    const int numberOfCellEdges = cf->getNumberOfEdges();
    QSet<Vertex *> faceVertices;
    for (int i = 0; i < numberOfCellEdges; i++) {
        Vertex *cv = cf->getVertex(i);
        if (!faceVertices.contains(cv))
            faceVertices.insert(cv);
        else {
            qDebug() << "CellFace" << cf->getId() << "has self intersecting boundary in CellVertex" << cv->getId() << ". Parameterization not possible!";
            return 0;
        }
    }

    Parameterization *param = new Parameterization();

    param->m_borderType = bt;
    param->m_cellFace = cf;
    param->m_originalMesh = originalMesh;
    param->m_parameterInterval = parameterInterval;
    param->fillVertexLocalIdMap();

    QVector<PolylineVertex *> borderVertices;
    QVector<int> corners;
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = (CellEdge *) cf->getEdge(i);
        const int numberOfEdgeVertices = ce->getPolylineVertices()->size();
        corners.append(borderVertices.size());
        //Add border vertices of that edge to the list (with exception of the last one, which equals the first vertex of the next edge)
        if (cf->edgeIsInverted(i))
            for (int j = numberOfEdgeVertices - 1; j > 0; j--)
                borderVertices.append(ce->getPolylineVertices()->at(j));
        else
            for (int j = 0; j < numberOfEdgeVertices - 1; j++)
                borderVertices.append(ce->getPolylineVertices()->at(j));
    }

    const int numberOfInnerVertices = cf->getMeshVertices()->size();
    const int numberOfBorderVertices = borderVertices.size();
    const int numberOfVertices = numberOfInnerVertices + numberOfBorderVertices;

    //Setup Eigen coefficient matrices and solution vectors
    Eigen::SparseMatrix<double> matrixBorder(numberOfVertices, numberOfBorderVertices);
    Eigen::SparseMatrix<double> matrixInner(numberOfVertices, numberOfInnerVertices);
    matrixInner.reserve(Eigen::VectorXi::Constant(numberOfInnerVertices, expectedEdgesPerVertex));      //Average 7, Expected Min/Max: 5/11
    matrixBorder.reserve(Eigen::VectorXi::Constant(numberOfBorderVertices, expectedEdgesPerVertex));    //Average 4, Expected Min/Max: 3/4
    Eigen::VectorXd borderValuesX(numberOfBorderVertices);
    Eigen::VectorXd borderValuesY(numberOfBorderVertices);

    Parameterization::calculateBorderProjection(&borderVertices, &corners, borderValuesX, borderValuesY, bt, parameterInterval);

    //Matrix coefficients:
    //
    //        / -w_ij                ,if i != j and Edge (i,j) exists
    // A_ij = | \sum_{k != i} w_ik   ,if i == j
    //        \ 0                    , otherwise

    //Write edges between inner vertices into matrix
    const int numMeshEdges = cf->getMeshEdges()->size();
    QVector<double> diagonalElements(numberOfVertices, 0);

    for (int i = 0; i < numMeshEdges; i++) {    //This part takes about as much time as the solver itself
        Edge *e = cf->getMeshEdges()->at(i);
        Vertex *v0 = e->getVertex(0);
        Vertex *v1 = e->getVertex(1);

        const int id0 = param->m_vertexLocalIdMap[v0];
        const int id1 = param->m_vertexLocalIdMap[v1];
        const double weight = param->getWeight(e, wt);

        diagonalElements[id0] += weight;
        diagonalElements[id1] += weight;
        matrixInner.insert(id0, id1) = -weight;
        matrixInner.insert(id1, id0) = -weight;
    }

    //Write edges and partial edges between border vertices and inner vertices into matrix
    for (int i = 0; i < numberOfBorderVertices; i++) {
        PolylineVertex *plv = borderVertices[i];
        if (plv->isMeshVertex()) {  //plv is an actual vertex -> need to look at incident edges of the original mesh
            Vertex *meshVertex = plv->getMeshVertex();
            foreach (int eId, *plv->getMeshVertex()->getIncidentEdgeIds()) {
                Edge *e = originalMesh->getEdge(eId);
                Vertex *v = e->getConnectedVertex(meshVertex);

                if (param->m_vertexLocalIdMap.contains(v)) {
                    int neighborId = param->m_vertexLocalIdMap[v];
                    if (neighborId < numberOfInnerVertices) {
                        const double weight = param->getWeight(e, wt);

                        diagonalElements[numberOfInnerVertices + i] += weight;
                        diagonalElements[neighborId] += weight;
                        matrixInner.insert(numberOfInnerVertices + i, neighborId) = -weight;
                        matrixBorder.insert(neighborId, i) = -weight;
                    }
                }
            }
        } else if (plv->isCrossingVertex()) {    //plv is only an artificial vertex -> only look at a part of the crossing edge
            Edge *e = plv->getCrossingEdge();
            for (int j = 0; j < 2; j++) {
                if (param->m_vertexLocalIdMap.contains(e->getVertex(j))) {
                    const int neighborId = param->m_vertexLocalIdMap[e->getVertex(j)];
                    if (neighborId < numberOfInnerVertices) {
                        const double lengthInner = (e->getVertex(j)->getPosition() - plv->getPosition()).length();
                        const double lengthTotal = (e->getVertex(0)->getPosition() - e->getVertex(1)->getPosition()).length();
                        const double weight = param->getWeight(e, wt) * lengthTotal / lengthInner; // = 1/(lengthInner/lengthTotal) -> short edge -> large weight

                        diagonalElements[numberOfInnerVertices + i] += weight;
                        diagonalElements[neighborId] += weight;
                        matrixInner.insert(numberOfInnerVertices + i, neighborId) = -weight;
                        matrixBorder.insert(neighborId, i) = -weight;
                    }
                }
            }
        }
    }

    //Write edges between neighboring border vertices into matrix (has low impact on the result)
    for (int i = 0; i < numberOfBorderVertices; i++) {
        const int nextId = (i + 1) % numberOfBorderVertices;
        //Better version, but still buggy:
//        Edge *tmpEdge = new Edge(borderVertices[i], borderVertices[nextId], -1);
//        Vertex *outerVertex = 0;
//        QVector<Vertex *> connectedVertices0, connectedVertices1;
//        if (borderVertices[i]->isCrossingVertex())
//            connectedVertices0 << borderVertices[i]->getCrossingEdge()->getVertex(0) << borderVertices[i]->getCrossingEdge()->getVertex(1);
//        if (borderVertices[i]->isMeshVertex())
//            foreach (int eId, *borderVertices[i]->getIncidentEdgeIds())
//                connectedVertices0 << originalMesh->getEdge(eId)->getConnectedVertex(borderVertices[i]->getMeshVertex());
//        if (borderVertices[nextId]->isCrossingVertex())
//            connectedVertices1 << borderVertices[nextId]->getCrossingEdge()->getVertex(0) << borderVertices[nextId]->getCrossingEdge()->getVertex(1);
//        if (borderVertices[nextId]->isMeshVertex())
//            foreach (int eId, *borderVertices[nextId]->getIncidentEdgeIds())
//                connectedVertices1 << originalMesh->getEdge(eId)->getConnectedVertex(borderVertices[nextId]->getMeshVertex());
//        foreach (Vertex *v, connectedVertices0)
//            if (connectedVertices1.contains(v)) {
//                outerVertex = v;
//                break;
//            }
//        const double weight = param->getWeight(tmpEdge, wt, outerVertex);
//        delete tmpEdge;

        const double weight = 1;

        diagonalElements[numberOfInnerVertices + i] += weight;
        diagonalElements[numberOfInnerVertices + nextId] += weight;
        matrixBorder.insert(numberOfInnerVertices + i, nextId) = -weight;
        matrixBorder.insert(numberOfInnerVertices + nextId, i) = -weight;
    }

    //Write diagonal elements into matrix
    for (int i = 0; i < numberOfInnerVertices; i++)
        matrixInner.insert(i, i) = diagonalElements[i];

    for (int i = 0; i < numberOfBorderVertices; i++)
        matrixBorder.insert(numberOfInnerVertices + i, i) = diagonalElements[numberOfInnerVertices + i];

    //Solve matrixInner * innerValues = - matrixBorder * borderValues
//    matrixInner.makeCompressed();
//    matrixBorder.makeCompressed();
//    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver(matrixInner);
//    Eigen::VectorXd innerValuesX = solver.solve(-matrixBorder*borderValuesX);
//    Eigen::VectorXd innerValuesY = solver.solve(-matrixBorder*borderValuesY);

    //Solve matrixInner^T * matrixInner * innerValues = - matrixInner^T * matrixBorder * borderValues:
    Eigen::SparseMatrix<double, Eigen::ColMajor> mm = matrixInner.transpose() * matrixInner;
    Eigen::SparseMatrix<double> mb = matrixInner.transpose() * matrixBorder;
    mm.makeCompressed();
    mb.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solverSquaredSystem(mm);
    Eigen::VectorXd innerValuesX = solverSquaredSystem.solve(-mb*borderValuesX);
    Eigen::VectorXd innerValuesY = solverSquaredSystem.solve(-mb*borderValuesY);

    //Write parameter values into Parameterization-object
    param->m_parameterValuesU = Eigen::VectorXd(numberOfVertices);
    param->m_parameterValuesV = Eigen::VectorXd(numberOfVertices);
    param->m_parameterValuesU << innerValuesX, borderValuesX;
    param->m_parameterValuesV << innerValuesY, borderValuesY;

    param->m_numberOfVertices = numberOfVertices;
    param->m_numberOfBorderVertices = numberOfBorderVertices;
    param->m_numberOfInnerVertices = numberOfInnerVertices;
    param->m_info = QString("Harmonic Map with weight type %1 and border type %2, parameter interval: [%3,%4]\n")
            .arg(wt).arg(bt).arg(parameterInterval.first).arg(parameterInterval.second);

    return param;
}

double Parameterization::getWeight(Edge *e, Parameterization::WeightType wt, Vertex *optionalOuterTriangleVertex)
{
    switch (wt) {
    case Parameterization::WT_constant:
        return 1;
    case Parameterization::WT_length:
        return (e->getVertex(0)->getPosition() - e->getVertex(1)->getPosition()).length();
    case Parameterization::WT_cotangent: {
        double result = 0;
        for (int i = 0; i < 2; i++) {
            const int fId = e->getIncidentFaceId(i);
            if (fId != -1 && !optionalOuterTriangleVertex) {
                Face *f = m_originalMesh->getFace(fId);
                if (f->getNumberOfEdges() != 3)
                    return 1;

                Vertex *v0 = e->getVertex(0);
                Vertex *v1 = e->getVertex(1);
                Vertex *vOut = 0;
                for (int k = 0; !vOut && k < 3; k++) {
                    vOut = f->getVertex(k);
                    if (vOut == v0 || vOut == v1)
                        vOut = 0;
                }
                QVector3D vOut_v0 = v0->getPosition() - vOut->getPosition();
                QVector3D vOut_v1 = v1->getPosition() - vOut->getPosition();
                const double productLength = vOut_v0.length() * vOut_v1.length();
                const double cos = QVector3D::dotProduct(vOut_v0, vOut_v1)/productLength;
                const double sin = QVector3D::crossProduct(vOut_v0, vOut_v1).length()/productLength;
                result += 0.5 * cos/sin; // = cotangent
            } else if (optionalOuterTriangleVertex) {
                Vertex *v0 = e->getVertex(0);
                Vertex *v1 = e->getVertex(1);
                Vertex *vOut = optionalOuterTriangleVertex;

                QVector3D vOut_v0 = v0->getPosition() - vOut->getPosition();
                QVector3D vOut_v1 = v1->getPosition() - vOut->getPosition();
                const double productLength = vOut_v0.length() * vOut_v1.length();
                const double cos = QVector3D::dotProduct(vOut_v0, vOut_v1)/productLength;
                const double sin = QVector3D::crossProduct(vOut_v0, vOut_v1).length()/productLength;
                result += 0.5 * cos/sin; // = cotangent
            }
        }
        return result;
    }
    default:
        return 1;
    }
}

void Parameterization::calculateBorderProjection(QVector<PolylineVertex *> *borderVertices, QVector<int> *corners, Eigen::VectorXd &projectedX, Eigen::VectorXd &projectedY, Parameterization::BorderType bt, QPair<double, double> parameterInterval)
{
    const int numberOfBorderVertices = borderVertices->size();
    const double a = parameterInterval.first;
    const double b = parameterInterval.second;
    switch (bt) {
    case Parameterization::BT_equidistantCircle: {
        //Project border onto a circle
        for (int i = 0; i < numberOfBorderVertices; i++) {
            projectedX[i] = a + (b - a) * cos(2*M_PI*(double) i / numberOfBorderVertices);
            projectedY[i] = a + (b - a) * sin(2*M_PI*(double) i / numberOfBorderVertices);
        }
        break;
    }
    case BT_lengthWeightedCircle: {
        const double minValue = 0.001;
        double totalLength = (borderVertices->first()->getPosition() - borderVertices->last()->getPosition()).length();
        if (totalLength < minValue)
            totalLength = minValue;
        for (int i = 1; i < numberOfBorderVertices; i++) {
            double dist = (borderVertices->at(i)->getPosition() - borderVertices->at(i-1)->getPosition()).length();
            if (dist < minValue)
                dist = minValue;
            totalLength += dist;
        }

        double currentLength = 0;
        projectedX[0] = a + (b - a) * (1 + cos(0))/2;
        projectedY[0] = a + (b - a) * (1 + sin(0))/2;
        for (int i = 1; i < numberOfBorderVertices; i++) {
            double dist = (borderVertices->at(i)->getPosition() - borderVertices->at(i-1)->getPosition()).length();
            if (dist < minValue)
                dist = minValue;
            currentLength += dist;
            projectedX[i] = a + (b - a) * (1 + cos(2*M_PI* currentLength / totalLength))/2;
            projectedY[i] = a + (b - a) * (1 + sin(2*M_PI* currentLength / totalLength))/2;
        }
        break;
    }
    case Parameterization::BT_equidistantPolygon: {
        const int numberOfCorners = corners->size();
        //Project corners onto equidistant circle
        QVector<QVector2D> cornerProjections(numberOfCorners);
        for (int i = 0; i < numberOfCorners; i++)
            cornerProjections[i] = QVector2D((1 + cos(2*M_PI*(double) i / numberOfCorners))/2, (1 + sin(2*M_PI*(double) i / numberOfCorners))/2);

        //Project points onto lines between two corners
        for (int i = 1; i < numberOfCorners; i++) {
            const int numVertices = corners->at(i) - corners->at(i-1);
            for (int j = 0; j < numVertices; j++) {
                const double lambda = (double) j / (double) numVertices;
                QVector2D pos = (1 - lambda) * cornerProjections[i - 1] + lambda * cornerProjections[i];
                projectedX[corners->at(i-1) + j] = a + (b - a) * pos.x();
                projectedY[corners->at(i-1) + j] = a + (b - a) * pos.y();
            }
        }
        const int numLastVertices = numberOfBorderVertices - corners->last();
        for (int j = 0; j < numLastVertices; j++) {
            const double lambda = (double) j / (double) numLastVertices;
            QVector2D pos = (1 - lambda) * cornerProjections.last() + lambda * cornerProjections.first();
            projectedX[corners->last() + j] = a + (b - a) * pos.x();
            projectedY[corners->last() + j] = a + (b - a) * pos.y();
        }

        break;
    }
    case Parameterization::BT_equidistantUnitSquare: {
        const int numberOfCorners = corners->size();
        if (numberOfCorners != 4)
            calculateBorderProjection(borderVertices, corners, projectedX, projectedY, BT_equidistantPolygon, parameterInterval);

        QVector<QVector2D> cornerProjections(4);
        cornerProjections[0] = QVector2D(a, a);
        cornerProjections[1] = QVector2D(b, a);
        cornerProjections[2] = QVector2D(b, b);
        cornerProjections[3] = QVector2D(a, b);

        //Project points onto lines between two corners
        for (int i = 1; i < numberOfCorners; i++) {
            const int numVertices = corners->at(i) - corners->at(i-1);
            for (int j = 0; j < numVertices; j++) {
                const double lambda = (double) j / (double) numVertices;
                QVector2D pos = (1 - lambda) * cornerProjections[i - 1] + lambda * cornerProjections[i];
                projectedX[corners->at(i-1) + j] = pos.x();
                projectedY[corners->at(i-1) + j] = pos.y();
            }
        }
        const int numLastVertices = numberOfBorderVertices - corners->last();
        for (int j = 0; j < numLastVertices; j++) {
            const double lambda = (double) j / (double) numLastVertices;
            QVector2D pos = (1 - lambda) * cornerProjections.last() + lambda * cornerProjections.first();
            projectedX[corners->last() + j] = pos.x();
            projectedY[corners->last() + j] = pos.y();
        }
        break;
    }
    case Parameterization::BT_lengthWeightedUnitSquare: {
        const int numberOfCorners = corners->size();
        if (numberOfCorners != 4)
            calculateBorderProjection(borderVertices, corners, projectedX, projectedY, BT_equidistantPolygon, parameterInterval);

        QVector<QVector2D> cornerProjections(4);
        cornerProjections[0] = QVector2D(a, a);
        cornerProjections[1] = QVector2D(b, a);
        cornerProjections[2] = QVector2D(b, b);
        cornerProjections[3] = QVector2D(a, b);

        //Project points onto lines between two corners
        for (int i = 1; i < numberOfCorners; i++) {
            double totalLength = 0;
            for (int j = corners->at(i-1) + 1; j < corners->at(i) + 1; j++)
                totalLength += (borderVertices->at(j-1)->getPosition() - borderVertices->at(j)->getPosition()).length();

            double currentLength = 0;
            projectedX[corners->at(i-1)] = cornerProjections[i-1].x();
            projectedY[corners->at(i-1)] = cornerProjections[i-1].y();
            for (int j = corners->at(i-1) + 1; j < corners->at(i); j++) {
                currentLength += (borderVertices->at(j-1)->getPosition() - borderVertices->at(j)->getPosition()).length();
                const double lambda = (double) currentLength / (double) totalLength;
                QVector2D pos = (1 - lambda) * cornerProjections[i - 1] + lambda * cornerProjections[i];
                projectedX[j] = pos.x();
                projectedY[j] = pos.y();
            }
        }

        //Handle last line seperately
        double totalLength = 0;
        for (int j = corners->last(); j < numberOfBorderVertices; j++)
            totalLength += (borderVertices->at(j-1)->getPosition() - borderVertices->at(j)->getPosition()).length();
        totalLength += (borderVertices->last()->getPosition() - borderVertices->first()->getPosition()).length();

        double currentLength = 0;
        projectedX[corners->last()] = cornerProjections.last().x();
        projectedY[corners->last()] = cornerProjections.last().y();
        for (int j = corners->last() + 1; j < numberOfBorderVertices; j++) {
            currentLength += (borderVertices->at(j-1)->getPosition() - borderVertices->at(j)->getPosition()).length();
            const double lambda = (double) currentLength / (double) totalLength;
            QVector2D pos = (1 - lambda) * cornerProjections.last() + lambda * cornerProjections.first();
            projectedX[j] = pos.x();
            projectedY[j] = pos.y();
        }
        break;
    }
    default:
        for (int i = 0; i < numberOfBorderVertices; i++) {
            projectedX[i] = 0;
            projectedY[i] = 0;
        }
    }
}

void Parameterization::fillVertexLocalIdMap()
{
    //Take inner vertices
    QVector<Vertex *> *meshVertices = m_cellFace->getMeshVertices();
    const int numberOfInnerVertices = meshVertices->size();

    //Assign local ids to inner vertices
    for (int i = 0; i < numberOfInnerVertices; i++)
        m_vertexLocalIdMap[meshVertices->at(i)] = i;

    int nextLocalId = numberOfInnerVertices;

    //Get border vertices
    const int numberOfCellEdges = m_cellFace->getNumberOfEdges();
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = (CellEdge *) m_cellFace->getEdge(i);
        const int numberOfEdgeVertices = ce->getPolylineVertices()->size();
        //For later queries, add also the cell vertices to the hash map
        m_vertexLocalIdMap[m_cellFace->getVertex(i)] = nextLocalId;
        //Add border vertices of that edge to the list (with exception of the last one, which equals the first vertex of the next edge)
        if (m_cellFace->edgeIsInverted(i)) {
            for (int j = numberOfEdgeVertices - 1; j > 0; j--) {
                PolylineVertex *plv = ce->getPolylineVertices()->at(j);
                m_vertexLocalIdMap[plv] = nextLocalId;
                if (plv->isMeshVertex())
                    m_vertexLocalIdMap[plv->getMeshVertex()] = nextLocalId;
                nextLocalId++;
            }
            if (i < numberOfCellEdges - 1) {
                m_vertexLocalIdMap[ce->getPolylineVertices()->first()] = nextLocalId;
                if (ce->getPolylineVertices()->first()->isMeshVertex())
                    m_vertexLocalIdMap[ce->getPolylineVertices()->first()] = nextLocalId;
            } else {
                m_vertexLocalIdMap[ce->getPolylineVertices()->first()] = numberOfInnerVertices;
                if (ce->getPolylineVertices()->first()->isMeshVertex())
                    m_vertexLocalIdMap[ce->getPolylineVertices()->first()] = numberOfInnerVertices;
            }
        } else {
            for (int j = 0; j < numberOfEdgeVertices - 1; j++) {
                PolylineVertex *plv = ce->getPolylineVertices()->at(j);
                m_vertexLocalIdMap[plv] = nextLocalId;
                if (plv->isMeshVertex())
                    m_vertexLocalIdMap[plv->getMeshVertex()] = nextLocalId;
                nextLocalId++;
            }
            if (i < numberOfCellEdges - 1) {
                m_vertexLocalIdMap[ce->getPolylineVertices()->last()] = nextLocalId;
                if (ce->getPolylineVertices()->last()->isMeshVertex())
                    m_vertexLocalIdMap[ce->getPolylineVertices()->last()] = nextLocalId;
            } else {
                m_vertexLocalIdMap[ce->getPolylineVertices()->last()] = numberOfInnerVertices;
                if (ce->getPolylineVertices()->last()->isMeshVertex())
                    m_vertexLocalIdMap[ce->getPolylineVertices()->last()] = numberOfInnerVertices;
            }
        }
    }

    m_parameterValuesU = Eigen::VectorXd(m_vertexLocalIdMap.keys().size());
    m_parameterValuesV = Eigen::VectorXd(m_vertexLocalIdMap.keys().size());
}
