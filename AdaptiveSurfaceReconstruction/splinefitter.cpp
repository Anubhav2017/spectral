#include "splinefitter.h"

#include "parameteroptimizer.h"

//c/c++ includes
#include <iostream>

//GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

//NLOpt
#include <nlopt.h>

//Qt
#include <QTime>
#include <QDebug>

Splinefitter::Splinefitter(CellMesh *cm, BSplineProperties *properties):
    m_cellMesh(cm), m_initialGuess(0),
    m_properties(properties),
    m_energyMatrix(Eigen::SparseMatrix<double>(36, 36)), m_weightLSTerm(0.5), m_weightEnergyTerm(0.5), m_useFixedWeights(true),
    m_numberOfG1VertexPoints(4), m_sizeOfG1VertexDomain(0.25), m_numberOfG1EdgePoints(2),
    m_g1VertexConstraintsAsOpt(true), m_g1VertexConstraintsOptFactor(1000), m_g1EdgeConstraintsAsOpt(true), m_g1EdgeConstraintsOptFactor(1000),
    m_nloMaxIterations(1000), m_nloFPrec(-1), m_nloXPrec(-1), m_nloG0constPrec(-1), m_nloEqualConstPrec(-1), m_nloG1constPrec(-1),
    m_detailedOutput(false)
{

}

struct TangentVectorCoeffientsNormal {
    Eigen::SparseVector<double> Nu0;
    Eigen::SparseVector<double> Nv0;
    Eigen::SparseVector<double> Nu1;
    Eigen::SparseVector<double> Nv1;
};

struct NormalConstraintParameters {
    int numberOfControlPoints;
    int id;
    TangentVectorCoeffientsNormal tangents;
};

double normalConstraint_normalized(unsigned int n, const double *x, double *grad, void *params) {
    NormalConstraintParameters *parameters = (NormalConstraintParameters *) params;
    const int numControlPoints = parameters->numberOfControlPoints;

    Eigen::VectorXd ctrPointsX(numControlPoints);
    Eigen::VectorXd ctrPointsY(numControlPoints);
    Eigen::VectorXd ctrPointsZ(numControlPoints);
    for (int i = 0; i < numControlPoints; i++) {
        ctrPointsX[i] = x[i                       ];
        ctrPointsY[i] = x[i +     numControlPoints];
        ctrPointsZ[i] = x[i + 2 * numControlPoints];
    }

    const QVector3D tu0(ctrPointsX.transpose() * parameters->tangents.Nu0,
                        ctrPointsY.transpose() * parameters->tangents.Nu0,
                        ctrPointsZ.transpose() * parameters->tangents.Nu0);
    const QVector3D tv0(ctrPointsX.transpose() * parameters->tangents.Nv0,
                        ctrPointsY.transpose() * parameters->tangents.Nv0,
                        ctrPointsZ.transpose() * parameters->tangents.Nv0);
    const QVector3D tu1(ctrPointsX.transpose() * parameters->tangents.Nu1,
                        ctrPointsY.transpose() * parameters->tangents.Nu1,
                        ctrPointsZ.transpose() * parameters->tangents.Nu1);
    const QVector3D tv1(ctrPointsX.transpose() * parameters->tangents.Nv1,
                        ctrPointsY.transpose() * parameters->tangents.Nv1,
                        ctrPointsZ.transpose() * parameters->tangents.Nv1);

    const QVector3D normal0 = QVector3D::crossProduct(tu0, tv0);
    const QVector3D normal1 = QVector3D::crossProduct(tu1, tv1);
    const QVector3D normalizedNormal0 = normal0.normalized();
    const QVector3D normalizedNormal1 = normal1.normalized();
    const QVector3D normalizedNormalDif = normalizedNormal0 - normalizedNormal1;
    const double normal0Length = normal0.length();
    const double normal0LengthSquared = normal0Length * normal0Length;
    const double normal1Length = normal1.length();
    const double normal1LengthSquared = normal1Length * normal1Length;
    if (grad) {
        for (int i = 0; i < numControlPoints; i++) {
            const double nu0_i = parameters->tangents.Nu0.coeff(i);
            const double nv0_i = parameters->tangents.Nv0.coeff(i);
            const double nu1_i = parameters->tangents.Nu1.coeff(i);
            const double nv1_i = parameters->tangents.Nv1.coeff(i);

            if (nu0_i != 0 || nv0_i != 0 || nu1_i != 0 || nv1_i != 0) {

                const QVector3D normal0Dxi(0, tu0.z() * nv0_i - nu0_i * tv0.z(), nu0_i * tv0.y() - tu0.y() * nv0_i);
                const QVector3D normal0Dyi(nu0_i * tv0.z() - tu0.z() * nv0_i, 0, tu0.x() * nv0_i - nu0_i * tv0.x());
                const QVector3D normal0Dzi(tu0.y() * nv0_i - nu0_i * tv0.y(), nu0_i * tv0.x() - tu0.x() * nv0_i, 0);
                const QVector3D normal1Dxi(0, tu1.z() * nv1_i - nu1_i * tv1.z(), nu1_i * tv1.y() - tu1.y() * nv1_i);
                const QVector3D normal1Dyi(nu1_i * tv1.z() - tu1.z() * nv1_i, 0, tu1.x() * nv1_i - nu1_i * tv1.x());
                const QVector3D normal1Dzi(tu1.y() * nv1_i - nu1_i * tv1.y(), nu1_i * tv1.x() - tu1.x() * nv1_i, 0);

                const QVector3D normalizedNormal0Dxi(normal0Dxi/normal0Length -
                                                     normalizedNormal0 * QVector3D::dotProduct(normal0Dxi, normal0)/normal0LengthSquared);
                const QVector3D normalizedNormal0Dyi(normal0Dyi/normal0Length -
                                                     normalizedNormal0 * QVector3D::dotProduct(normal0Dyi, normal0)/normal0LengthSquared);
                const QVector3D normalizedNormal0Dzi(normal0Dzi/normal0Length -
                                                     normalizedNormal0 * QVector3D::dotProduct(normal0Dzi, normal0)/normal0LengthSquared);

                const QVector3D normalizedNormal1Dxi(normal1Dxi/normal1Length -
                                                     normalizedNormal1 * QVector3D::dotProduct(normal1Dxi, normal1)/normal1LengthSquared);
                const QVector3D normalizedNormal1Dyi(normal1Dyi/normal1Length -
                                                     normalizedNormal1 * QVector3D::dotProduct(normal1Dyi, normal1)/normal1LengthSquared);
                const QVector3D normalizedNormal1Dzi(normal1Dzi/normal1Length -
                                                     normalizedNormal1 * QVector3D::dotProduct(normal1Dzi, normal1)/normal1LengthSquared);


                grad[i                       ] = QVector3D::dotProduct(normalizedNormal0Dxi - normalizedNormal1Dxi, normalizedNormalDif);
                grad[i +     numControlPoints] = QVector3D::dotProduct(normalizedNormal0Dyi - normalizedNormal1Dyi, normalizedNormalDif);
                grad[i + 2 * numControlPoints] = QVector3D::dotProduct(normalizedNormal0Dzi - normalizedNormal1Dzi, normalizedNormalDif);
            } else {
                grad[i                       ] = 0;
                grad[i +     numControlPoints] = 0;
                grad[i + 2 * numControlPoints] = 0;
            }
        }

        for (unsigned int i = 3 * numControlPoints; i < n; i++)
            grad[i] = 0;
    }

    return 0.5 * normalizedNormalDif.lengthSquared();
}

struct NonlinearOptParameters {
    Eigen::SparseMatrix<double> *M;
    Eigen::SparseMatrix<double, Eigen::RowMajor> *N;
    Eigen::SparseMatrix<double> *matGradientFirstHalf;
    Eigen::VectorXd *dataX;
    Eigen::VectorXd *dataY;
    Eigen::VectorXd *dataZ;
    Eigen::VectorXd *NDataX;
    Eigen::VectorXd *NDataY;
    Eigen::VectorXd *NDataZ;
    double weightLS;
    double weightEnergy;
    QVector<NormalConstraintParameters *> normalConstraints;
    double weightNormalConstraints;
};

double nonlinearOptFunction(unsigned int n, const double *x, double *grad, void *params) {
    NonlinearOptParameters *parameters = (NonlinearOptParameters *) params;

    const int numControlPoints = parameters->N->cols();
    const int numControlPoints3 = 3 * numControlPoints;

    Eigen::VectorXd controlPoints(numControlPoints3);
    for (int i = 0; i < numControlPoints3; i++)
        controlPoints[i] = x[i];

    const Eigen::VectorXd &ctrPointsX = controlPoints.segment(0                   , numControlPoints);
    const Eigen::VectorXd &ctrPointsY = controlPoints.segment(numControlPoints    , numControlPoints);
    const Eigen::VectorXd &ctrPointsZ = controlPoints.segment(numControlPoints * 2, numControlPoints);
//    for (int i = 0; i < numControlPoints; i++) {
//        ctrPointsX[i] = controlPoints[i                       ];
//        ctrPointsY[i] = controlPoints[i +     numControlPoints];
//        ctrPointsZ[i] = controlPoints[i + 2 * numControlPoints];
//    }

    const Eigen::VectorXd errorX = (*parameters->N) * ctrPointsX - (*parameters->dataX);
    const Eigen::VectorXd errorY = (*parameters->N) * ctrPointsY - (*parameters->dataY);
    const Eigen::VectorXd errorZ = (*parameters->N) * ctrPointsZ - (*parameters->dataZ);
    const double lsError =  0.5 * parameters->weightLS * (errorX.dot(errorX) + errorY.dot(errorY) + errorZ.dot(errorZ));

    double energy = ctrPointsX.transpose() * ((*parameters->M) * ctrPointsX);
    energy += ctrPointsY.transpose() * ((*parameters->M) * ctrPointsY);
    energy += ctrPointsZ.transpose() * ((*parameters->M) * ctrPointsZ);
    energy *= 0.5 * parameters->weightEnergy;

    if (grad) {
        //First version
//        const Eigen::VectorXd gradX = (parameters->weightLS * (*parameters->N).transpose() * (*parameters->N) + parameters->weightEnergy * (*parameters->M)) * ctrPointsX
//                                    - parameters->weightLS * (*parameters->N).transpose() * (*parameters->dataX);
//        const Eigen::VectorXd gradY = (parameters->weightLS * (*parameters->N).transpose() * (*parameters->N) + parameters->weightEnergy * (*parameters->M)) * ctrPointsY
//                                    - parameters->weightLS * (*parameters->N).transpose() * (*parameters->dataY);
//        const Eigen::VectorXd gradZ = (parameters->weightLS * (*parameters->N).transpose() * (*parameters->N) + parameters->weightEnergy * (*parameters->M)) * ctrPointsZ
//                                    - parameters->weightLS * (*parameters->N).transpose() * (*parameters->dataZ);
        //Second Version - calculate matrix product once per iteration
//        Eigen::SparseMatrix<double> mat = (parameters->weightLS * (*parameters->N).transpose() * (*parameters->N) + parameters->weightEnergy * (*parameters->M));
        //Third Version - calculate matrix product only once before optimization
//        const Eigen::VectorXd gradX = *parameters->matGradientFirstHalf * ctrPointsX - parameters->weightLS * (*parameters->N).transpose() * (*parameters->dataX);
//        const Eigen::VectorXd gradY = *parameters->matGradientFirstHalf * ctrPointsY - parameters->weightLS * (*parameters->N).transpose() * (*parameters->dataY);
//        const Eigen::VectorXd gradZ = *parameters->matGradientFirstHalf * ctrPointsZ - parameters->weightLS * (*parameters->N).transpose() * (*parameters->dataZ);
        //Fourth Version - calculate all fixed products only once
        const Eigen::VectorXd gradX = *parameters->matGradientFirstHalf * ctrPointsX - *parameters->NDataX;
        const Eigen::VectorXd gradY = *parameters->matGradientFirstHalf * ctrPointsY - *parameters->NDataY;
        const Eigen::VectorXd gradZ = *parameters->matGradientFirstHalf * ctrPointsZ - *parameters->NDataZ;

        for (int i = 0; i < numControlPoints; i++) {
            grad[i]                        = gradX[i];
            grad[i + numControlPoints]     = gradY[i];
            grad[i + 2 * numControlPoints] = gradZ[i];
        }

        for (unsigned int i = numControlPoints3; i < n; i++)
             grad[i] = 0;
    }

    const int numNormalConstraints = parameters->normalConstraints.size();
    const double normalFactor = parameters->weightNormalConstraints;
    double normalError = 0;
    double *gradConstraint = 0;
    if (grad)
        gradConstraint = new double[n];
    for (int i = 0; i < numNormalConstraints; i++) {
        normalError += normalConstraint_normalized(n, x, gradConstraint ,parameters->normalConstraints.at(i));
        if (grad)
            for (unsigned int j = 0; j < n; j++)
                grad[j] += normalFactor * gradConstraint[j];
    }
    if (grad)
        delete[] gradConstraint;

    normalError *= normalFactor;

    return lsError + energy + normalError;
}

struct G0ConstraintParameters {
    int index0;
    int coef0;
    int index1;
    int coef1;
};

double g0Constraint(unsigned int n, const double *x, double *grad, void *params) {
    G0ConstraintParameters *parameters = (G0ConstraintParameters *) params;

    const int index0 = parameters->index0;
    const int index1 = parameters->index1;
    const int coef0 = parameters->coef0;
    const int coef1 = parameters->coef1;
    if (grad) {
        for (unsigned int i = 0; i < n; i++)
            grad[i] = 0;

        grad[index0] = coef0;
        grad[index1] = coef1;
    }

    return coef0 * x[index0] + coef1 * x[index1];
}

double g0QuadraticConstraint(unsigned int n, const double *x, double *grad, void *params) {
    G0ConstraintParameters *parameters = (G0ConstraintParameters *) params;

    const int index0 = parameters->index0;
    const int index1 = parameters->index1;
    const int coef0 = parameters->coef0;
    const int coef1 = parameters->coef1;
    const double sum = coef0 * x[index0] + coef1 * x[index1];
    if (grad) {
        for (unsigned int i = 0; i < n; i++)
            grad[i] = 0;

        grad[index0] = coef0 * sum;
        grad[index1] = coef1 * sum;
    }


    return 0.5 * sum * sum;
}




double normalConstraint_alphaNormal(unsigned int n, const double *x, double *grad, void *params) {
    NormalConstraintParameters *parameters = (NormalConstraintParameters *) params;
    const int numControlPoints = parameters->numberOfControlPoints;
    const int constraintId = parameters->id;

    Eigen::VectorXd ctrPointsX(numControlPoints);
    Eigen::VectorXd ctrPointsY(numControlPoints);
    Eigen::VectorXd ctrPointsZ(numControlPoints);
    for (int i = 0; i < numControlPoints; i++) {
        ctrPointsX[i] = x[i                       ];
        ctrPointsY[i] = x[i +     numControlPoints];
        ctrPointsZ[i] = x[i + 2 * numControlPoints];
    }

    const QVector3D tu0(ctrPointsX.transpose() * parameters->tangents.Nu0,
                        ctrPointsY.transpose() * parameters->tangents.Nu0,
                        ctrPointsZ.transpose() * parameters->tangents.Nu0);
    const QVector3D tv0(ctrPointsX.transpose() * parameters->tangents.Nv0,
                        ctrPointsY.transpose() * parameters->tangents.Nv0,
                        ctrPointsZ.transpose() * parameters->tangents.Nv0);
    const QVector3D tu1(ctrPointsX.transpose() * parameters->tangents.Nu1,
                        ctrPointsY.transpose() * parameters->tangents.Nu1,
                        ctrPointsZ.transpose() * parameters->tangents.Nu1);
    const QVector3D tv1(ctrPointsX.transpose() * parameters->tangents.Nv1,
                        ctrPointsY.transpose() * parameters->tangents.Nv1,
                        ctrPointsZ.transpose() * parameters->tangents.Nv1);

    const double alpha = x[3 * numControlPoints + constraintId];
    const double alphaSquared = alpha * alpha;
    const QVector3D normal0 = QVector3D::crossProduct(tu0, tv0);
    const QVector3D normal1 = QVector3D::crossProduct(tu1, tv1);
    const QVector3D alphaNormalDist = normal0 - alphaSquared * normal1;

    if (grad) {
        for (int i = 0; i < numControlPoints; i++) {
            const double nu0_i = parameters->tangents.Nu0.coeff(i);
            const double nv0_i = parameters->tangents.Nv0.coeff(i);
            const double nu1_i = parameters->tangents.Nu1.coeff(i);
            const double nv1_i = parameters->tangents.Nv1.coeff(i);

            const QVector3D normal0Dxi(0, tu0.z() * nv0_i - nu0_i * tv0.z(), nu0_i * tv0.y() - tu0.y() * nv0_i);
            const QVector3D normal0Dyi(nu0_i * tv0.z() - tu0.z() * nv0_i, 0, tu0.x() * nv0_i - nu0_i * tv0.x());
            const QVector3D normal0Dzi(tu0.y() * nv0_i - nu0_i * tv0.y(), nu0_i * tv0.x() - tu0.x() * nv0_i, 0);
            const QVector3D normal1Dxi(0, tu1.z() * nv1_i - nu1_i * tv1.z(), nu1_i * tv1.y() - tu1.y() * nv1_i);
            const QVector3D normal1Dyi(nu1_i * tv1.z() - tu1.z() * nv1_i, 0, tu1.x() * nv1_i - nu1_i * tv1.x());
            const QVector3D normal1Dzi(tu1.y() * nv1_i - nu1_i * tv1.y(), nu1_i * tv1.x() - tu1.x() * nv1_i, 0);

            grad[i                       ] = QVector3D::dotProduct(normal0Dxi - alphaSquared * normal1Dxi, alphaNormalDist);
            grad[i +     numControlPoints] = QVector3D::dotProduct(normal0Dyi - alphaSquared * normal1Dyi, alphaNormalDist);
            grad[i + 2 * numControlPoints] = QVector3D::dotProduct(normal0Dzi - alphaSquared * normal1Dzi, alphaNormalDist);
        }

        for (unsigned int i = 3 * numControlPoints; i < n; i++)
            grad[i] = 0;

        grad[3 * numControlPoints + constraintId] = QVector3D::dotProduct(- 2 * alpha * normal1, alphaNormalDist);
    }

    return 0.5 * alphaNormalDist.lengthSquared();
}

double normalConstraint_unnormalized(unsigned int n, const double *x, double *grad, void *params) {
    NormalConstraintParameters *parameters = (NormalConstraintParameters *) params;
    const int numControlPoints = parameters->numberOfControlPoints;

    Eigen::VectorXd ctrPointsX(numControlPoints);
    Eigen::VectorXd ctrPointsY(numControlPoints);
    Eigen::VectorXd ctrPointsZ(numControlPoints);
    for (int i = 0; i < numControlPoints; i++) {
        ctrPointsX[i] = x[i                       ];
        ctrPointsY[i] = x[i +     numControlPoints];
        ctrPointsZ[i] = x[i + 2 * numControlPoints];
    }

    const QVector3D tu0(ctrPointsX.transpose() * parameters->tangents.Nu0,
                        ctrPointsY.transpose() * parameters->tangents.Nu0,
                        ctrPointsZ.transpose() * parameters->tangents.Nu0);
    const QVector3D tv0(ctrPointsX.transpose() * parameters->tangents.Nv0,
                        ctrPointsY.transpose() * parameters->tangents.Nv0,
                        ctrPointsZ.transpose() * parameters->tangents.Nv0);
    const QVector3D tu1(ctrPointsX.transpose() * parameters->tangents.Nu1,
                        ctrPointsY.transpose() * parameters->tangents.Nu1,
                        ctrPointsZ.transpose() * parameters->tangents.Nu1);
    const QVector3D tv1(ctrPointsX.transpose() * parameters->tangents.Nv1,
                        ctrPointsY.transpose() * parameters->tangents.Nv1,
                        ctrPointsZ.transpose() * parameters->tangents.Nv1);

    const QVector3D normal0 = QVector3D::crossProduct(tu0, tv0);
    const QVector3D normal1 = QVector3D::crossProduct(tu1, tv1);
    const QVector3D normalDif = normal0 - normal1;;
    if (grad) {
        for (int i = 0; i < numControlPoints; i++) {
            const double nu0_i = parameters->tangents.Nu0.coeff(i);
            const double nv0_i = parameters->tangents.Nv0.coeff(i);
            const double nu1_i = parameters->tangents.Nu1.coeff(i);
            const double nv1_i = parameters->tangents.Nv1.coeff(i);

            const QVector3D normal0Dxi(0, tu0.z() * nv0_i - nu0_i * tv0.z(), nu0_i * tv0.y() - tu0.y() * nv0_i);
            const QVector3D normal0Dyi(nu0_i * tv0.z() - tu0.z() * nv0_i, 0, tu0.x() * nv0_i - nu0_i * tv0.x());
            const QVector3D normal0Dzi(tu0.y() * nv0_i - nu0_i * tv0.y(), nu0_i * tv0.x() - tu0.x() * nv0_i, 0);
            const QVector3D normal1Dxi(0, tu1.z() * nv1_i - nu1_i * tv1.z(), nu1_i * tv1.y() - tu1.y() * nv1_i);
            const QVector3D normal1Dyi(nu1_i * tv1.z() - tu1.z() * nv1_i, 0, tu1.x() * nv1_i - nu1_i * tv1.x());
            const QVector3D normal1Dzi(tu1.y() * nv1_i - nu1_i * tv1.y(), nu1_i * tv1.x() - tu1.x() * nv1_i, 0);

            grad[i                       ] = QVector3D::dotProduct(normal0Dxi - normal1Dxi, normalDif);
            grad[i +     numControlPoints] = QVector3D::dotProduct(normal0Dyi - normal1Dyi, normalDif);
            grad[i + 2 * numControlPoints] = QVector3D::dotProduct(normal0Dzi - normal1Dzi, normalDif);
        }

        for (unsigned int i = 3 * numControlPoints; i < n; i++)
            grad[i] = 0;
    }

    return 0.5 * normalDif.lengthSquared();
}

struct EqualConstraintParameters {
    int index;
    double value;
};

double equalConstraint(unsigned int n, const double *x, double *grad, void *params) {
    EqualConstraintParameters *parameters = (EqualConstraintParameters *) params;
    const int index = parameters->index;
    const double value = parameters->value;

    if (grad) {
        for (unsigned int i = 0; i < n; i++)
            grad[i] = 0;

        grad[index] = 1;
    }

    return x[index] - value;
}


BSplinePatchNetwork *Splinefitter::doApproximationLocalG0ConstraintG1Sample()
{
    QVector<BSplineSurface *> bSplineSurfaces;

    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();
    const int numberOfControlPointsPerPatch = numberOfControlPointsPerRow * numberOfControlPointsPerRow;
    //Control points around a cell vertex (only need to define this once)
    QVector<QPair<int, int> > CP_around_V[4] = {QVector<QPair<int, int> >(4), QVector<QPair<int, int> >(4), QVector<QPair<int, int> >(4), QVector<QPair<int, int> >(4)};
    CP_around_V[0][0] = QPair<int, int>(0, 0);
    CP_around_V[0][1] = QPair<int, int>(0, 1);
    CP_around_V[0][2] = QPair<int, int>(1, 0);
    CP_around_V[0][3] = QPair<int, int>(1, 1);
    CP_around_V[1][0] = QPair<int, int>(numberOfControlPointsPerRow - 1, 0);
    CP_around_V[1][1] = QPair<int, int>(numberOfControlPointsPerRow - 1, 1);
    CP_around_V[1][2] = QPair<int, int>(numberOfControlPointsPerRow - 2, 0);
    CP_around_V[1][3] = QPair<int, int>(numberOfControlPointsPerRow - 2, 1);
    CP_around_V[2][0] = QPair<int, int>(numberOfControlPointsPerRow - 1, numberOfControlPointsPerRow - 1);
    CP_around_V[2][1] = QPair<int, int>(numberOfControlPointsPerRow - 1, numberOfControlPointsPerRow - 2);
    CP_around_V[2][2] = QPair<int, int>(numberOfControlPointsPerRow - 2, numberOfControlPointsPerRow - 1);
    CP_around_V[2][3] = QPair<int, int>(numberOfControlPointsPerRow - 2, numberOfControlPointsPerRow - 2);
    CP_around_V[3][0] = QPair<int, int>(0, numberOfControlPointsPerRow - 1);
    CP_around_V[3][1] = QPair<int, int>(0, numberOfControlPointsPerRow - 2);
    CP_around_V[3][2] = QPair<int, int>(1, numberOfControlPointsPerRow - 1);
    CP_around_V[3][3] = QPair<int, int>(1, numberOfControlPointsPerRow - 2);


    const int numberOfCellFaces = m_cellMesh->getNumberOfFaces();
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    const int numberOfCellVertices = m_cellMesh->getNumberOfVertices();

    QVector<QVector<QVector<QVector3D> > >globalControlPoints(numberOfCellFaces, QVector<QVector<QVector3D> >
                                                              (numberOfControlPointsPerRow, QVector<QVector3D>
                                                               (numberOfControlPointsPerRow)));
    QVector<QVector<QVector<bool> > > computedControlPoints(numberOfCellFaces, QVector<QVector<bool> >
                                                            (numberOfControlPointsPerRow, QVector<bool>
                                                             (numberOfControlPointsPerRow, false)));

    for (int cfId = 0; cfId < numberOfCellFaces; cfId++) {
        CellFace *cf = m_cellMesh->getFace(cfId);
        globalControlPoints[cfId][0][0] = cf->getVertex(0)->getPosition();
        globalControlPoints[cfId][numberOfControlPointsPerRow - 1][0] = cf->getVertex(1)->getPosition();
        globalControlPoints[cfId][numberOfControlPointsPerRow - 1][numberOfControlPointsPerRow - 1] = cf->getVertex(2)->getPosition();
        globalControlPoints[cfId][0][numberOfControlPointsPerRow - 1] = cf->getVertex(3)->getPosition();
    }

    QVector<QVector3D> debug_vertexFitPoints;
    QVector<QVector3D> debug_edgeFitPoints;
    QVector<QVector3D> debug_innerFitPoints;

    QVector<bool> vertexFittingDone(numberOfCellVertices, false);
    QVector<bool> edgeFittingDone(numberOfCellEdges, false);
    QVector<bool> faceFittingDone(numberOfCellFaces, false);

    //Step one:1 Do B-Spline fitting around each Vertex
    #pragma omp parallel for
    for (int cvId = 0; cvId < numberOfCellVertices; cvId++) {
        QTime time;
        time.start();
        qDebug() << "CellVertex" << cvId;
        CellVertex *cv = m_cellMesh->getVertex(cvId);

        QList<CellEdge *> incidentEdges;
        foreach (int ceId, *cv->getIncidentEdgeIds())
            incidentEdges.append(m_cellMesh->getEdge(ceId));

        const int numberOfIncidentFaces = incidentEdges.size();
        const int numberOfLocalControlPoints = numberOfIncidentFaces * numberOfControlPointsPerPatch;

        //Create (sorted) list of incident face ids
        QVector<int> incidentFaceIds;
        QList<CellEdge *> incidentEdgesTmp(incidentEdges);
        CellEdge *tmpCe = incidentEdgesTmp.takeFirst();
        incidentFaceIds << tmpCe->getIncidentFaceId(0) << tmpCe->getIncidentFaceId(1);
        QVector<CellEdge *> sortedCellEdges;
        sortedCellEdges << tmpCe;

        while (!incidentEdgesTmp.isEmpty()) {
            int lastFaceId = incidentFaceIds.last();
            for (int i = 0; i < incidentEdgesTmp.size(); i++) {
                tmpCe = incidentEdgesTmp[i];
                for (int j = 0; j < 2 && tmpCe; j++) {
                    if (tmpCe->getIncidentFaceId(j) == lastFaceId) {
                        incidentEdgesTmp.removeAt(i);
                        lastFaceId = tmpCe->getIncidentFaceId((j + 1) % 2);
                        incidentFaceIds.append(lastFaceId);
                        sortedCellEdges.append(tmpCe);
                        tmpCe = 0;
                    }
                }
            }
        }

        //Cannot skip if two subsequent cell edges have to be g1
        bool skip = !(sortedCellEdges.first()->isG1Edge() && sortedCellEdges.last()->isG1Edge());
        for (int i = 1; i < numberOfIncidentFaces; i++)
            if (sortedCellEdges[i-1]->isG1Edge() && sortedCellEdges[i]->isG1Edge())
                skip = false;

        if (!skip) {
            //count the total number of datapoints in cellfaces around the cellvertex
            int numberOfLocalDataPoints = 0;
            QVector<CellFace *> incidentFaces(numberOfIncidentFaces);
            for (int i = 0; i < numberOfIncidentFaces; i++) {
                CellFace* f = m_cellMesh->getFace(incidentFaceIds[i]);
                incidentFaces[i] = f;
                numberOfLocalDataPoints += f->getMeshVertices()->size();
            }

            //Construct coefficient matrix and data point vector for approximation (least square system)
            Eigen::SparseMatrix<double, Eigen::RowMajor> matrixLS(numberOfLocalDataPoints, numberOfLocalControlPoints);
            const int nonZeroBasisFunctions = m_properties->getOrder() + 1;
            matrixLS.reserve(Eigen::VectorXi::Constant(numberOfLocalDataPoints, nonZeroBasisFunctions * nonZeroBasisFunctions));
            Eigen::VectorXd dataPointsX(numberOfLocalDataPoints);
            Eigen::VectorXd dataPointsY(numberOfLocalDataPoints);
            Eigen::VectorXd dataPointsZ(numberOfLocalDataPoints);

            int lineOffset = 0;
            for (int cfNumber = 0; cfNumber < numberOfIncidentFaces; cfNumber++) {
                CellFace *cf = incidentFaces[cfNumber];
                Parameterization *cfParam = m_cellMesh->getQuadCellParameterization(cf->getId());
                QVector<Vertex *> *data = cf->getMeshVertices();

                const int numberOfPatchDataPoints = data->size();
                const int columnOffset = numberOfControlPointsPerPatch * cfNumber;

                for (int k = 0; k < numberOfPatchDataPoints; k++) {
                    Vertex *v = data->at(k);
                    QVector2D vParam = cfParam->getParameter(v);

                    dataPointsX[lineOffset + k] = v->getPosition().x();
                    dataPointsY[lineOffset + k] = v->getPosition().y();
                    dataPointsZ[lineOffset + k] = v->getPosition().z();

                    for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                        const double nv = m_properties->N(vParam.y(), j);
                        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                            const double nu = m_properties->N(vParam.x(), i);
                            const double product = nu * nv;
                            if (product != 0)
                                matrixLS.insert(lineOffset + k, columnOffset + BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow)) = product;
                        }
                    }
                }
                lineOffset += numberOfPatchDataPoints;
            }

            //construct matrix for G0 continuity constraints
            Eigen::SparseMatrix<double, Eigen::RowMajor> matrixG0(numberOfControlPointsPerRow * numberOfIncidentFaces, numberOfLocalControlPoints);
            matrixG0.reserve(Eigen::VectorXi::Constant(numberOfControlPointsPerRow * numberOfIncidentFaces, 2));
            Eigen::SparseMatrix<double, Eigen::RowMajor> matrixG1((numberOfControlPointsPerRow - 2) * numberOfIncidentFaces, numberOfLocalControlPoints);
            matrixG1.reserve(Eigen::VectorXi::Constant((numberOfControlPointsPerRow - 2) * numberOfIncidentFaces, 4));
            for (int iFace = 0; iFace < numberOfIncidentFaces; iFace++) {
                CellFace *cf[2] = {incidentFaces.at(iFace), incidentFaces.at((iFace+1)%numberOfIncidentFaces)};
                CellEdge *ce = (CellEdge *) cf[0]->getSharedEdge(cf[1]);

                const int internalCeId[2] = {cf[0]->getInternalEdgeId(ce), cf[1]->getInternalEdgeId(ce)};
                const int offset[2] = {iFace * numberOfControlPointsPerPatch, ((iFace+1)%numberOfIncidentFaces) * numberOfControlPointsPerPatch};

                int i[2];
                int j[2];
                int iInner[2];
                int jInner[2];
                int *var[2];
                int step[2];

                for (int k = 0; k < 2; k++) {
                    switch (internalCeId[k]) {
                    case 0:
                        i[k] = 0;
                        j[k] = 0;
                        var[k] = &i[k];
                        step[k] = 1;
                        iInner[k] = 0;
                        jInner[k] = 1;
                        break;
                    case 1:
                        i[k] = numberOfControlPointsPerRow - 1;
                        j[k] = 0;
                        var[k] = &j[k];
                        step[k] = 1;
                        iInner[k] = -1;
                        jInner[k] = 0;
                        break;
                    case 2:
                        i[k] = numberOfControlPointsPerRow - 1;
                        j[k] = numberOfControlPointsPerRow - 1;
                        var[k] = &i[k];
                        step[k] = -1;
                        iInner[k] = 0;
                        jInner[k] = -1;
                        break;
                    case 3:
                        i[k] = 0;
                        j[k] = numberOfControlPointsPerRow - 1;
                        var[k] = &j[k];
                        step[k] = -1;
                        iInner[k] = 1;
                        jInner[k] = 0;
                        break;
                    default: qDebug() << "ERROR in index computation along celledge: Internal Index of Edge is not in [0,3] ->" << internalCeId;
                    }
                    if (cf[k]->edgeIsInverted(internalCeId[k])) {
                        *(var[k]) = numberOfControlPointsPerRow - 1 - *(var[k]);
                        step[k] *= -1;
                    }
                }

                for (int iCPCounter = 0; iCPCounter < numberOfControlPointsPerRow; iCPCounter++) {
                    int globalCPIndex[2];
                    int globalCPInnerIndex[2];
                    for (int k = 0; k < 2; k++) {
                        globalCPIndex[k] = offset[k] + BSplineSurface::matToVecIndexLocal(i[k], j[k], numberOfControlPointsPerRow);
                        globalCPInnerIndex[k] = offset[k] + BSplineSurface::matToVecIndexLocal(i[k] + iInner[k], j[k] + jInner[k], numberOfControlPointsPerRow);
                        *(var[k]) += step[k];
                    }
                    //qDebug() << globalCPIndex[0] << "," << i[0] << j[0] << " | " << globalCPIndex[1] << "," << i[1] << j[1];

                    const int indexRow = iFace * numberOfControlPointsPerRow + iCPCounter;
                    matrixG0.insert(indexRow, globalCPIndex[0]) =  1;
                    matrixG0.insert(indexRow, globalCPIndex[1]) = -1;

                    if (iCPCounter > 0 && iCPCounter < numberOfControlPointsPerRow - 1) {
                        const int indexRowG1 = iFace * (numberOfControlPointsPerRow - 2) + iCPCounter - 1;
                        matrixG1.insert(indexRowG1, globalCPIndex[0]) = -1;
                        matrixG1.insert(indexRowG1, globalCPIndex[1]) = -1;
                        matrixG1.insert(indexRowG1, globalCPInnerIndex[0]) = 1;
                        matrixG1.insert(indexRowG1, globalCPInnerIndex[1]) = 1;
                    }
                }

            }

            //Extract important indices
            QVector<QVector<int> > importantCtrPointIndices(numberOfIncidentFaces, QVector<int>(4));
            for (int i = 0; i < numberOfIncidentFaces; i++) {
                CellFace *cf = incidentFaces[i];

                const int internalCVId = cf->getInternalVertexId(cv);   //value is in {0,1,2,3}
                const int offset = i * numberOfControlPointsPerPatch;

                for (int k = 0; k < 4; k++)
                    importantCtrPointIndices[i][k] = offset + BSplineSurface::matToVecIndexLocal(CP_around_V[internalCVId][k], numberOfControlPointsPerRow);
            }

            //Very ugly hack to avoid linear dependencies between rows of matrixG0
            bool removed = false;
            Eigen::SparseMatrix<double, Eigen::RowMajor> matrixG0new(matrixG0.rows() - 1, matrixG0.cols());
            matrixG0new.reserve(Eigen::VectorXi::Constant(matrixG0.rows() - 1, 2));

            for (int i = 0; i < matrixG0.rows(); i++) {
                if (!removed) {
                    if (matrixG0.coeff(i, importantCtrPointIndices[0][0]) != 0)
                        removed = true;
                    else
                        for (int j = 0; j < matrixG0.cols(); j++)
                            if (matrixG0.coeff(i, j))
                                matrixG0new.insert(i, j) = matrixG0.coeff(i, j);

                } else {
                    for (int j = 0; j < matrixG0.cols(); j++)
                        if (matrixG0.coeff(i, j))
                            matrixG0new.insert(i - 1, j) = matrixG0.coeff(i, j);
                }
            }
            matrixG0 = matrixG0new;

            QVector<TangentVectorCoeffientsNormal> tangentCoefficients;
            //Extract tangent vectors for the normal constraints
//            for (int iCe = 0; iCe < numberOfIncidentFaces - 1; iCe++) {
//                CellEdge *ce = incidentEdges[iCe];

//                Eigen::SparseVector<double> Nu[2];
//                Eigen::SparseVector<double> Nv[2];

//                Eigen::SparseVector<double> NuTwist[2];
//                Eigen::SparseVector<double> NvTwist[2];

//                for (int k = 0; k < 2; k++) {
//                    CellFace *cf = m_cellMesh->getFace(ce->getIncidentFaceId(k));
//                    const int internalCVId = cf->getInternalVertexId(cv);
//                    const int cfNumber = incidentFaces.indexOf(cf); //TODO that's not nice, but ok as long as the number of incident faces is small

//                    Nu[k].resize(numberOfLocalControlPoints);
//                    Nv[k].resize(numberOfLocalControlPoints);

//                    const int internalCeId = cf->getInternalEdgeId(ce);
//                    NuTwist[k].resize(numberOfLocalControlPoints);
//                    NvTwist[k].resize(numberOfLocalControlPoints);
//                    switch (internalCVId) {
//                    case 0:
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][2]) = 1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][1]) = 1;;
//                        if (internalCeId == 0) {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = 1;
//                        } else {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = 1;
//                        }
//                        break;
//                    case 1:
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][2]) = -1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][1]) = 1;
//                        if (internalCeId == 0) {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = 1;
//                        } else {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = 1;
//                        }
//                        break;
//                    case 2:
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][2]) = -1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][1]) = -1;
//                        if (internalCeId == 2) {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = 1;
//                        } else {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                        }
//                        break;
//                    case 3:
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                        Nu[k].insert(importantCtrPointIndices[cfNumber][2]) = 1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                        Nv[k].insert(importantCtrPointIndices[cfNumber][1]) = -1;
//                        if (internalCeId == 2) {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][2]) = 1;
//                        } else {
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = -1;
//                            NuTwist[k].insert(importantCtrPointIndices[cfNumber][3]) = 1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][1]) = -1;
//                            NvTwist[k].insert(importantCtrPointIndices[cfNumber][0]) = 1;
//                        }
//                        break;
//                    default:
//                        qDebug() << "ERROR in B-spline approximation: Internal index of vertex is not in [0, 3] ->" << internalCVId;
//                    }
//                }

//                TangentVectorCoeffientsNormal tanCoefficients;
//                tanCoefficients.Nu0 = Nu[0];
//                tanCoefficients.Nv0 = Nv[0];
//                tanCoefficients.Nu1 = Nu[1];
//                tanCoefficients.Nv1 = Nv[1];
//                tangentCoefficients.append(tanCoefficients);

//                TangentVectorCoeffientsNormal tanCoefficientsTwist;
//                tanCoefficientsTwist.Nu0 = NuTwist[0];
//                tanCoefficientsTwist.Nv0 = NvTwist[0];
//                tanCoefficientsTwist.Nu1 = NuTwist[1];
//                tanCoefficientsTwist.Nv1 = NvTwist[1];
//                tangentCoefficients.append(tanCoefficientsTwist);
//            }

            //Construct G1 constraints in a small neighborhood around CV
            const double varStep = m_sizeOfG1VertexDomain / (double) (m_numberOfG1VertexPoints - 1);
            for (int iCe = 0; iCe < numberOfIncidentFaces; iCe++) {
                CellEdge *ce = incidentEdges[iCe];
                if (ce->isG1Edge()) {
                    const int cfNumbers[2] = {incidentFaceIds.indexOf(ce->getIncidentFaceId(0)), incidentFaceIds.indexOf(ce->getIncidentFaceId(1))};
                    const int offset[2] = {cfNumbers[0] * numberOfControlPointsPerPatch, cfNumbers[1] * numberOfControlPointsPerPatch};

                    CellFace *incidentCEFaces[2] = {incidentFaces[cfNumbers[0]], incidentFaces[cfNumbers[1]]};
                    const int internalCVIds[2] = {incidentCEFaces[0]->getInternalVertexId(cv), incidentCEFaces[1]->getInternalVertexId(cv)};
                    const int internalCEIds[2] = {incidentCEFaces[0]->getInternalEdgeId(ce), incidentCEFaces[1]->getInternalEdgeId(ce)};

                    double u[2];
                    double v[2];

                    double *varParam[2];
                    double varDirection[2];


                    for (int k = 0; k < 2; k++) {
                        switch (internalCVIds[k]) {
                        case 0:
                            u[k] = 0;
                            v[k] = 0;
                            break;
                        case 1:
                            u[k] = 1;
                            v[k] = 0;
                            break;
                        case 2:
                            u[k] = 1;
                            v[k] = 1;
                            break;
                        case 3:
                            u[k] = 0;
                            v[k] = 1;
                            break;
                        default:
                            qDebug() << "ERROR: Internal CV ID not in [0,3]!";
                            exit(-1);
                        }

                        switch (internalCEIds[k]) {
                        case 0:
                            varParam[k] = &u[k];
                            break;
                        case 1:
                            varParam[k] = &v[k];
                            break;
                        case 2:
                            varParam[k] = &u[k];
                            break;
                        case 3:
                            varParam[k] = &v[k];
                            break;
                        default:
                            qDebug() << "ERROR: Internal CE ID not in [0,3]!";
                            exit(-1);
                        }

                        if (*(varParam[k]) == 0)
                            varDirection[k] = 1;
                        else
                            varDirection[k] = -1;
                    }

                    for (int iEval = 0; iEval < m_numberOfG1VertexPoints; iEval++) {
                        Eigen::SparseVector<double> NTu0(numberOfLocalControlPoints);
                        Eigen::SparseVector<double> NTu1(numberOfLocalControlPoints);
                        Eigen::SparseVector<double> NTv0(numberOfLocalControlPoints);
                        Eigen::SparseVector<double> NTv1(numberOfLocalControlPoints);
                        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                            const double nu0 = m_properties->N(u[0], i);
                            const double nu0d = m_properties->NDerivative(u[0], i, 1);
                            const double nu1d = m_properties->NDerivative(u[1], i, 1);
                            const double nu1 = m_properties->N(u[1], i);
                            for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                                const int ctrPointId = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);

                                NTu0.insert(offset[0] + ctrPointId) = nu0d * m_properties->N(v[0], j);
                                NTv0.insert(offset[0] + ctrPointId) = nu0 * m_properties->NDerivative(v[0], j, 1);
                                NTu1.insert(offset[1] + ctrPointId) = nu1d * m_properties->N(v[1], j);
                                NTv1.insert(offset[1] + ctrPointId) = nu1 * m_properties->NDerivative(v[1], j, 1);
                            }
                        }

                        TangentVectorCoeffientsNormal tanCoefficients;
                        tanCoefficients.Nu0 = NTu0;
                        tanCoefficients.Nv0 = NTv0;
                        tanCoefficients.Nu1 = NTu1;
                        tanCoefficients.Nv1 = NTv1;

                        tangentCoefficients.append(tanCoefficients);

                        for (int k = 0; k < 2; k++) {
                            *(varParam[k]) += varDirection[k] * varStep;
                        }
                    }
                }
            }
            const int numberOfNormalConstraints = tangentCoefficients.size();

            //Calculate or use fixed weights
            double weightLS = m_weightLSTerm;
            double weightEnergy = m_weightEnergyTerm;
            if (!m_useFixedWeights) {
                double totalSurfaceArea = 0;
                for (int i = 0; i < numberOfIncidentFaces; i++)
                    totalSurfaceArea += incidentFaces.at(i)->calculateAreaByCorners();
                weightLS = (double) 1 / numberOfLocalDataPoints;
                weightEnergy = 1 / totalSurfaceArea;
            }

            //Create energy matrix for full system (contains individual energy matrices on its diagonals)
            Eigen::SparseMatrix<double> matrixEnergyAllPatches(numberOfLocalControlPoints, numberOfLocalControlPoints);
            matrixG0new.reserve(Eigen::VectorXi::Constant(numberOfLocalControlPoints, numberOfControlPointsPerPatch));
            for (int k = 0; k < numberOfIncidentFaces; k++)
                for (int i = 0; i < numberOfControlPointsPerPatch; i++)
                    for (int j = 0; j < numberOfControlPointsPerPatch; j++) {
                        const double coef = m_energyMatrix.coeff(i, j);
                        if (coef)
                            matrixEnergyAllPatches.insert(k * numberOfControlPointsPerPatch + i, k * numberOfControlPointsPerPatch + j) = m_energyMatrix.coeff(i, j);
                    }

//            Eigen::EigenSolver<Eigen::MatrixXd> eSolver1;
//            eSolver1.compute(matrixEnergyAllPatches.toDense(), false);
//            Eigen::VectorXcd eigenvalues1 = eSolver1.eigenvalues();
//            std::cout << "Eigenvalues" << std::endl << eigenvalues1 << std::endl;


            //Calculate A^T * A + diag(M, ... ,M)
            Eigen::SparseMatrix<double> matrixLStLS = weightLS * matrixLS.transpose() * matrixLS + weightEnergy * matrixEnergyAllPatches;

            Eigen::VectorXd rightSideX = weightLS * matrixLS.transpose() * dataPointsX;
            Eigen::VectorXd rightSideY = weightLS * matrixLS.transpose() * dataPointsY;
            Eigen::VectorXd rightSideZ = weightLS * matrixLS.transpose() * dataPointsZ;

            const int matG0Rows = matrixG0.rows();
            const int matG1Rows = matrixG1.rows();

            Eigen::VectorXd initialGuessLSG0G1X(numberOfLocalControlPoints);    //Allocation only needed when initial guess is used
            Eigen::VectorXd initialGuessLSG0G1Y(numberOfLocalControlPoints);
            Eigen::VectorXd initialGuessLSG0G1Z(numberOfLocalControlPoints);
            if (!m_initialGuess) {
                //calculate initial guess as solution of LS-fitting with G0 and simplified G1 constraints
                Eigen::SparseMatrix<double> matLSSystemG0G1(numberOfLocalControlPoints + matG0Rows + matG1Rows, numberOfLocalControlPoints + matG0Rows + matG1Rows);
                for (int i = 0; i < numberOfLocalControlPoints; i++)
                    for (int j = 0; j < numberOfLocalControlPoints; j++)
                        matLSSystemG0G1.insert(i, j) = matrixLStLS.coeff(i, j);

                for (int i = 0; i < matG0Rows; i++) {
                    for (int j = 0; j < numberOfLocalControlPoints; j++) {
                        matLSSystemG0G1.insert(numberOfLocalControlPoints + i, j) = matrixG0.coeff(i, j);
                        matLSSystemG0G1.insert(j, numberOfLocalControlPoints + i) = matrixG0.coeff(i, j);
                    }
                }

                for (int i = 0; i < matG1Rows; i++) {
                    for (int j = 0; j < numberOfLocalControlPoints; j++) {
                        matLSSystemG0G1.insert(numberOfLocalControlPoints + matG0Rows + i, j) = matrixG1.coeff(i, j);
                        matLSSystemG0G1.insert(j, numberOfLocalControlPoints + matG0Rows + i) = matrixG1.coeff(i, j);
                    }
                }

                Eigen::VectorXd rightSideLSSystemG0X(numberOfLocalControlPoints + matG0Rows + matG1Rows);
                Eigen::VectorXd rightSideLSSystemG0Y(numberOfLocalControlPoints + matG0Rows + matG1Rows);
                Eigen::VectorXd rightSideLSSystemG0Z(numberOfLocalControlPoints + matG0Rows + matG1Rows);
                for (int i = 0; i < numberOfLocalControlPoints; i++) {
                    rightSideLSSystemG0X[i] = rightSideX[i];
                    rightSideLSSystemG0Y[i] = rightSideY[i];
                    rightSideLSSystemG0Z[i] = rightSideZ[i];
                }
                for (int i = numberOfLocalControlPoints; i < numberOfLocalControlPoints + matG0Rows + matG1Rows; i++) {
                    rightSideLSSystemG0X[i] = 0;
                    rightSideLSSystemG0Y[i] = 0;
                    rightSideLSSystemG0Z[i] = 0;
                }

                Eigen::SparseLU<Eigen::SparseMatrix<double> > initialSolver;
                initialSolver.compute(matLSSystemG0G1);
                initialGuessLSG0G1X = initialSolver.solve(rightSideLSSystemG0X);
                initialGuessLSG0G1Y = initialSolver.solve(rightSideLSSystemG0Y);
                initialGuessLSG0G1Z = initialSolver.solve(rightSideLSSystemG0Z);
            } else {
                for (int cfNumber = 0; cfNumber < numberOfIncidentFaces; cfNumber++) {
                    for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                            const int ctrPointId = cfNumber * numberOfControlPointsPerPatch + BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
                            const QVector3D ctrPointInitialGuess = m_initialGuess->getBSplineSurface(incidentFaceIds[cfNumber])->getControlPointsP()->at(i).at(j);
                            initialGuessLSG0G1X[ctrPointId] = ctrPointInitialGuess.x();
                            initialGuessLSG0G1Y[ctrPointId] = ctrPointInitialGuess.y();
                            initialGuessLSG0G1Z[ctrPointId] = ctrPointInitialGuess.z();
                        }
                    }
                }
            }

            //Setup parameters for nlopt
            NonlinearOptParameters optParameters;
            optParameters.M = &matrixEnergyAllPatches;
            optParameters.N = &matrixLS;
            optParameters.matGradientFirstHalf = &matrixLStLS;
            optParameters.dataX = &dataPointsX;
            optParameters.dataY = &dataPointsY;
            optParameters.dataZ = &dataPointsZ;
            optParameters.NDataX = &rightSideX;
            optParameters.NDataY = &rightSideY;
            optParameters.NDataZ = &rightSideZ;
            optParameters.weightLS = weightLS;
            optParameters.weightEnergy = weightEnergy;
            optParameters.weightNormalConstraints = m_g1VertexConstraintsOptFactor; //possibly unused but should not slow it down

            const int numOptVariables = 3 * numberOfLocalControlPoints;

            //setup optimizer
            nlopt_opt optimizer = nlopt_create(NLOPT_LD_SLSQP, numOptVariables);
            nlopt_set_min_objective(optimizer, nonlinearOptFunction, &optParameters);
            nlopt_set_xtol_rel(optimizer, m_nloXPrec);
            nlopt_set_ftol_rel(optimizer, m_nloFPrec);
            nlopt_set_maxeval(optimizer, m_nloMaxIterations);

            //setup G0 constraints
            QVector<G0ConstraintParameters *> g0Constraints(3 * matG0Rows);
            for (int i = 0; i < matG0Rows; i++) {
                int index0 = 0;
                while (matrixG0.coeff(i, index0) == 0)
                    index0++;
                int index1 = index0 + 1;
                while (matrixG0.coeff(i, index1) == 0)
                    index1++;

                for (int k = 0; k < 3; k++) {
                    G0ConstraintParameters *g0Params = new G0ConstraintParameters;
                    g0Params->index0 = index0 + (k * numberOfLocalControlPoints);
                    g0Params->index1 = index1 + (k * numberOfLocalControlPoints);
                    g0Params->coef0 = 1;
                    g0Params->coef1 = -1;
                    nlopt_add_equality_constraint(optimizer, g0Constraint, g0Params, m_nloG0constPrec);
                    g0Constraints[3 * i + k] = g0Params;
                }
            }

            //setup normal constraints (either as real constraints or as optimization values)
            QVector<NormalConstraintParameters *> normalConstraints(numberOfNormalConstraints);
            for (int i = 0; i < numberOfNormalConstraints; i++) {
                NormalConstraintParameters *normalParams = new NormalConstraintParameters;
                normalParams->numberOfControlPoints = numberOfLocalControlPoints;
                normalParams->tangents = tangentCoefficients[i];
                normalParams->id = i;

                normalConstraints[i] = normalParams;
                if (!m_g1VertexConstraintsAsOpt)
                    nlopt_add_equality_constraint(optimizer, normalConstraint_normalized, normalParams, m_nloG1constPrec);
            }
            if (m_g1VertexConstraintsAsOpt)
                optParameters.normalConstraints = normalConstraints;

            //set initial guess
            double optVariables[numOptVariables];
            for (int i = 0; i < numberOfLocalControlPoints; i++) {
                optVariables[i                                 ] = initialGuessLSG0G1X[i];
                optVariables[i +     numberOfLocalControlPoints] = initialGuessLSG0G1Y[i];
                optVariables[i + 2 * numberOfLocalControlPoints] = initialGuessLSG0G1Z[i];
            }
//            for (int i = 3 * numberOfLocalControlPoints; i < 3 * numberOfLocalControlPoints + tangentCoefficients.size(); i++)
//                optVariables[i] = 1;

            //do optimization
            double optValue;
            int status = nlopt_optimize(optimizer, optVariables, &optValue);
            if (status < 0) {
                qDebug() << "ERROR" << status << "Value:" << optValue;
            } else {
                qDebug() << "Optimization done! Vertex:" << cvId << "Status:" << status << "Value:" << optValue;
            }

            //clean up optimizer
            nlopt_destroy(optimizer);
            for (int i = 0; i < 3 * matG0Rows; i++)
                delete g0Constraints[i];
            for (int i = 0; i < numberOfNormalConstraints; i++)
                delete normalConstraints[i];

            //extract the result (only the control points)
            Eigen::VectorXd ctrX(numberOfLocalControlPoints);
            Eigen::VectorXd ctrY(numberOfLocalControlPoints);
            Eigen::VectorXd ctrZ(numberOfLocalControlPoints);

            for (int i = 0; i < numberOfLocalControlPoints; i++) {
                ctrX[i] = optVariables[i];
                ctrY[i] = optVariables[i + numberOfLocalControlPoints];
                ctrZ[i] = optVariables[i + 2 * numberOfLocalControlPoints];
            }

            //Remove numerical errors in G0 constraints by using the average of the G0 points
            //First the center vertex, because more than two vertices need to be equal here
            //TODO save the G0 errors before overwriting
            QVector3D centerPoint;
            for (int i = 0; i < numberOfIncidentFaces; i++) {
                const int ctrPointIndex = importantCtrPointIndices[i][0];
                centerPoint += QVector3D(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);
            }
            centerPoint /= numberOfIncidentFaces;
            for (int i = 0; i < numberOfIncidentFaces; i++) {
                const int ctrPointIndex = importantCtrPointIndices[i][0];
                ctrX[ctrPointIndex] = centerPoint.x();
                ctrY[ctrPointIndex] = centerPoint.y();
                ctrZ[ctrPointIndex] = centerPoint.z();
            }

            //Now the vertices along the edges
            for (int i = 0; i < matrixG0.rows(); i++) {
                int index0 = 0;
                while (matrixG0.coeff(i, index0) == 0)
                    index0++;
                int index1 = index0 + 1;
                while (matrixG0.coeff(i, index1) == 0)
                    index1++;
                const double averageX = (ctrX[index0] + ctrX[index1])/2;
                ctrX[index0] = averageX;
                ctrX[index1] = averageX;
                const double averageY = (ctrY[index0] + ctrY[index1])/2;
                ctrY[index0] = averageY;
                ctrY[index1] = averageY;
                const double averageZ = (ctrZ[index0] + ctrZ[index1])/2;
                ctrZ[index0] = averageZ;
                ctrZ[index1] = averageZ;
            }

            //Store the G1-vertex-relevant control points in global control point matrix
            for (int cfNumber = 0; cfNumber < numberOfIncidentFaces; cfNumber++) {
                for (int i = 0; i < 4; i++) {
                    const int ctrPointIndex = importantCtrPointIndices[cfNumber][i];
                    const int localIndex = ctrPointIndex - (cfNumber * numberOfControlPointsPerPatch);
                    int index_i, index_j;
                    BSplineSurface::vecToMatIndexLocal(localIndex, index_i, index_j, numberOfControlPointsPerRow);
                    const QVector3D ctrPoint = QVector3D(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);

                    globalControlPoints[incidentFaceIds[cfNumber]][index_i][index_j] = ctrPoint;
                    computedControlPoints[incidentFaceIds[cfNumber]][index_i][index_j] = true;
                    debug_vertexFitPoints << ctrPoint;
                }
            }

            qDebug() << "time:" << time.elapsed();

            //==============//
            // OUTPUT Start //
            //==============//
            const double errorX = (matrixLS * ctrX - dataPointsX).norm();
            const double errorY = (matrixLS * ctrY - dataPointsY).norm();
            const double errorZ = (matrixLS * ctrZ - dataPointsZ).norm();
            const double lsError = errorX * errorX + errorY * errorY + errorZ * errorZ;
            const double energyX = ctrX.transpose() * matrixEnergyAllPatches * ctrX;
            const double energyY = ctrY.transpose() * matrixEnergyAllPatches * ctrY;
            const double energyZ = ctrZ.transpose() * matrixEnergyAllPatches * ctrZ;
            std::cout << "Squared error: " << lsError << " | per datapoint: " << lsError / numberOfLocalDataPoints << std::endl;
            std::cout << "Energy: " << energyX + energyY + energyZ << " x, y, z: " << energyX << " " << energyY << " " << energyZ << std::endl;

            //G0 errors will be 0 (due to postprocessing)
            //std::cout << "G0 errors [x]: " << (matrixG0 * ctrX).transpose() << " [y]: " << (matrixG0 * ctrY).transpose() << " [z]: " << (matrixG0 * ctrZ).transpose() << std::endl;

            QVector<double> normalErrors(numberOfNormalConstraints);
            double maxNormalError = -1;
            double averageNormalError = 0;
            for (int i = 0; i < numberOfNormalConstraints; i++) {
                const QVector3D Tu0(tangentCoefficients.at(i).Nu0.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nu0.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nu0.transpose() * ctrZ);
                const QVector3D Tv0(tangentCoefficients.at(i).Nv0.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nv0.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nv0.transpose() * ctrZ);
                const QVector3D Tu1(tangentCoefficients.at(i).Nu1.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nu1.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nu1.transpose() * ctrZ);
                const QVector3D Tv1(tangentCoefficients.at(i).Nv1.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nv1.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nv1.transpose() * ctrZ);

                const double normalError = (QVector3D::crossProduct(Tu0, Tv0).normalized() - QVector3D::crossProduct(Tu1, Tv1).normalized()).length();

                normalErrors[i] = normalError;
                averageNormalError += normalError;
                if (normalError > maxNormalError)
                    maxNormalError = normalError;
            }
            averageNormalError /= (double) numberOfNormalConstraints;
            std::cout << "Normal errors (" << numberOfNormalConstraints << " constraints): ";
            if (m_detailedOutput) {
                for (int i = 0; i < numberOfNormalConstraints; i++) {
                    //std::cout << "[" << normalError.x() << " " << normalError.y() << " " << normalError.z() << "](" << normalError.length() << ") ";
                    std::cout << "(" << normalErrors[i] << ") ";
                }
            }
            std::cout << "Average: " << averageNormalError << " Max: " << maxNormalError << std::endl;

            //============//
            // OUTPUT END //
            //============//

            ///=============///
            /// DEBUG START ///
            ///=============///
//            bSplineSurfaces.clear();
//            for (int cfNumber = 0; cfNumber < numberOfIncidentFaces; cfNumber++) {
//                QVector<QVector<QVector3D> > patchControlPoints(m_numberOfControlPointsPerRow, QVector<QVector3D>(m_numberOfControlPointsPerRow));
//                for (int i = 0; i < m_numberOfControlPointsPerRow; i++) {
//                    for (int j = 0; j < m_numberOfControlPointsPerRow; j++) {
//                        const int ctrPointIndex = cfNumber * m_numberOfControlPointsPerPatch + matToVecIndexLocal(i, j, m_numberOfControlPointsPerRow);
//                        patchControlPoints[i][j] = QVector3D(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);
//                    }
//                }

//                BSplineSurface *splinePatch = new BSplineSurface();
//                splinePatch->setKnotsU(knots);
//                splinePatch->setKnotsV(knots);
//                splinePatch->setOrder(m_splineOrder);
//                splinePatch->setControlPoints(patchControlPoints);
//                bSplineSurfaces.append(splinePatch);
//            }
//            m_stlWriter->writeBSplineSurfacesToStl(bSplineSurfaces, QString("vertexfit_%1.stl").arg(cvId), 51);
//            for (int cfNumber = 0; cfNumber < numberOfIncidentFaces; cfNumber++)
//                delete bSplineSurfaces[cfNumber];
//            bSplineSurfaces.clear();
            ///===========///
            /// DEBUG END ///
            ///===========///
        } else {
            qDebug() << "no G1 edge -> skip!";
        }

        vertexFittingDone[cvId] = true;
        int finishedVertexFits = 0;
        for (int i = 0; i < numberOfCellVertices; i++)
            if (vertexFittingDone[i])
                finishedVertexFits++;
        std::cout << "Vertex fits done: " << finishedVertexFits << "/" << numberOfCellVertices << std::endl;
    }

    //Next step: Do G1 fitting along the edges
    #pragma omp parallel for
    for (int ceId = 0; ceId < numberOfCellEdges; ceId++) {
        qDebug() << "CellEdge:" << ceId;

        CellEdge *ce = m_cellMesh->getEdge(ceId);

        if (ce->isG1Edge()) {
            const int incidentCfIds[2] = {ce->getIncidentFaceId(0), ce->getIncidentFaceId(1)};
            CellFace *incidentCellFaces[2] = {m_cellMesh->getFace(incidentCfIds[0]), m_cellMesh->getFace(incidentCfIds[1])};

            const int numberOfLocalControlPoints = 2 * numberOfControlPointsPerPatch;

            //count the total number of datapoints in cellfaces around the cellvertex
            const int numberOfLocalDataPoints = incidentCellFaces[0]->getMeshVertices()->size() + incidentCellFaces[1]->getMeshVertices()->size();

            //Construct coefficient matrix and data point vector for approximation (least square system)
            Eigen::SparseMatrix<double, Eigen::RowMajor> matrixLS(numberOfLocalDataPoints, numberOfLocalControlPoints);
            const int nonZeroBasisFunctions = m_properties->getOrder() + 1;
            matrixLS.reserve(Eigen::VectorXi::Constant(numberOfLocalDataPoints, nonZeroBasisFunctions * nonZeroBasisFunctions));
            Eigen::VectorXd dataPointsX(numberOfLocalDataPoints);
            Eigen::VectorXd dataPointsY(numberOfLocalDataPoints);
            Eigen::VectorXd dataPointsZ(numberOfLocalDataPoints);

            int lineOffset = 0;
            for (int cfNumber = 0; cfNumber < 2; cfNumber++) {
                CellFace *cf = incidentCellFaces[cfNumber];
                Parameterization *cfParam = m_cellMesh->getQuadCellParameterization(incidentCfIds[cfNumber]);
                QVector<Vertex *> *data = cf->getMeshVertices();

                const int numberOfPatchDataPoints = data->size();
                const int columnOffset = numberOfControlPointsPerPatch * cfNumber;

                for (int k = 0; k < numberOfPatchDataPoints; k++) {
                    Vertex *v = data->at(k);
                    QVector2D vParam = cfParam->getParameter(v);

                    dataPointsX[lineOffset + k] = v->getPosition().x();
                    dataPointsY[lineOffset + k] = v->getPosition().y();
                    dataPointsZ[lineOffset + k] = v->getPosition().z();

                    for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                        const double nv = m_properties->N(vParam.y(), j);
                        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                            const double nu = m_properties->N(vParam.x(), i);
                            const double product = nu * nv;
                            if (product != 0)
                                matrixLS.insert(lineOffset + k, columnOffset + BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow)) = product;
                        }
                    }
                }
                lineOffset += numberOfPatchDataPoints;
            }

            //construct matrix for G0 continuity constraints
            Eigen::SparseMatrix<double> matrixG0(numberOfControlPointsPerRow - 4, numberOfLocalControlPoints);
            Eigen::SparseMatrix<double> matrixG1(numberOfControlPointsPerRow - 4, numberOfLocalControlPoints);

            const int internalCeId[2] = {incidentCellFaces[0]->getInternalEdgeId(ce), incidentCellFaces[1]->getInternalEdgeId(ce)};
            const int offset[2] = {0, numberOfControlPointsPerPatch};

            int i[2];
            int j[2];
            int *var[2];
            int step[2];
            int iInner[2];
            int jInner[2];
            for (int k = 0; k < 2; k++) {
                switch (internalCeId[k]) {
                case 0:
                    i[k] = 0;
                    j[k] = 0;
                    var[k] = &i[k];
                    step[k] = 1;
                    iInner[k] = 0;
                    jInner[k] = 1;
                    break;
                case 1:
                    i[k] = numberOfControlPointsPerRow - 1;
                    j[k] = 0;
                    var[k] = &j[k];
                    step[k] = 1;
                    iInner[k] = -1;
                    jInner[k] = 0;
                    break;
                case 2:
                    i[k] = numberOfControlPointsPerRow - 1;
                    j[k] = numberOfControlPointsPerRow - 1;
                    var[k] = &i[k];
                    step[k] = -1;
                    iInner[k] = 0;
                    jInner[k] = -1;
                    break;
                case 3:
                    i[k] = 0;
                    j[k] = numberOfControlPointsPerRow - 1;
                    var[k] = &j[k];
                    step[k] = -1;
                    iInner[k] = 1;
                    jInner[k] = 0;
                    break;
                default: qDebug() << "ERROR in index computation along celledge: Internal Index of Edge is not in [0,3] ->" << internalCeId;
                }

                if (incidentCellFaces[k]->edgeIsInverted(internalCeId[k])) {
                    *(var[k]) = numberOfControlPointsPerRow - 1 - *(var[k]);
                    step[k] *= -1;
                }

                *(var[k]) += 2 * step[k];   //Otherwise we would start with the control points around the vertex (which we did in step 1)
            }


            QVector<QVector<int> > importantCtrPointIndices(2, QVector<int>(2 * (numberOfControlPointsPerRow - 4)));
            for (int iCPCounter = 2; iCPCounter < numberOfControlPointsPerRow - 2; iCPCounter++) {
                int globalCPIndex[2];
                int globalCPInnerIndex[2];
                for (int k = 0; k < 2; k++) {
                    globalCPIndex[k] = offset[k] + BSplineSurface::matToVecIndexLocal(i[k], j[k], numberOfControlPointsPerRow);
                    globalCPInnerIndex[k] = offset[k] + BSplineSurface::matToVecIndexLocal(i[k] + iInner[k], j[k] + jInner[k], numberOfControlPointsPerRow);
                    *(var[k]) += step[k];

                    importantCtrPointIndices[k][2 * (iCPCounter - 2)    ] = globalCPIndex[k];
                    importantCtrPointIndices[k][2 * (iCPCounter - 2) + 1] = globalCPInnerIndex[k];
                }

                const int indexRow = iCPCounter - 2;
                matrixG0.insert(indexRow, globalCPIndex[0]) =  1;
                matrixG0.insert(indexRow, globalCPIndex[1]) = -1;

                matrixG1.insert(indexRow, globalCPIndex[0]) = -1;
                matrixG1.insert(indexRow, globalCPIndex[1]) = -1;
                matrixG1.insert(indexRow, globalCPInnerIndex[0]) = 1;
                matrixG1.insert(indexRow, globalCPInnerIndex[1]) = 1;
            }

            //construct tripples of tangent vectors that need to be coplanar
            QVector<TangentVectorCoeffientsNormal> tangentCoefficients;
            for (int iEval = 0; iEval < m_numberOfG1EdgePoints; iEval++) {
                double t = (double) (iEval + 1) / (m_numberOfG1EdgePoints + 1);
                double u[2];
                double v[2];

                Eigen::SparseVector<double> NTu0(numberOfLocalControlPoints);
                Eigen::SparseVector<double> NTu1(numberOfLocalControlPoints);
                Eigen::SparseVector<double> NTv0(numberOfLocalControlPoints);
                Eigen::SparseVector<double> NTv1(numberOfLocalControlPoints);

                for (int k = 0; k < 2; k++) {
                    double tOriented = t;
                    if (incidentCellFaces[k]->edgeIsInverted(internalCeId[k]))
                        tOriented = 1 - t;
                    if (internalCeId[k] == 0) {
                        u[k] = tOriented;
                        v[k] = 0;
                    } else if (internalCeId[k] == 1) {
                        u[k] = 1;
                        v[k] = tOriented;
                    } else if (internalCeId[k] == 2) {
                        u[k] = 1 - tOriented;
                        v[k] = 1;
                    } else if (internalCeId[k] == 3) {
                        u[k] = 0;
                        v[k] = 1 - tOriented;
                    }
                }

                for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                    const double nu0 = m_properties->N(u[0], i);
                    const double nu0d = m_properties->NDerivative(u[0], i, 1);
                    const double nu1d = m_properties->NDerivative(u[1], i, 1);
                    const double nu1 = m_properties->N(u[1], i);
                    for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                        const int ctrPointId = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);

                        NTu0.insert(ctrPointId) = nu0d * m_properties->N(v[0], j);
                        NTv0.insert(ctrPointId) = nu0 * m_properties->NDerivative(v[0], j, 1);
                        NTu1.insert(offset[1] + ctrPointId) = nu1d * m_properties->N(v[1], j);
                        NTv1.insert(offset[1] + ctrPointId) = nu1 * m_properties->NDerivative(v[1], j, 1);
                    }
                }

                TangentVectorCoeffientsNormal tanCoefficients;
                tanCoefficients.Nu0 = NTu0;
                tanCoefficients.Nv0 = NTv0;
                tanCoefficients.Nu1 = NTu1;
                tanCoefficients.Nv1 = NTv1;

                tangentCoefficients.append(tanCoefficients);
            }
            const int numberOfNormalConstraints = tangentCoefficients.size();

            double weightLS = m_weightLSTerm;
            double weightEnergy = m_weightEnergyTerm;
            if (!m_useFixedWeights) {
                const double totalSurfaceArea = incidentCellFaces[0]->calculateAreaByCorners() + incidentCellFaces[1]->calculateAreaByCorners();
                weightLS = (double) 1 / numberOfLocalDataPoints;
                weightEnergy = 1 / totalSurfaceArea;
            }

            //Calculate A^T * A + diag(M, ... ,M)
            Eigen::SparseMatrix<double> matrixLStLS = weightLS * matrixLS.transpose() * matrixLS;
            for (int k = 0; k < 2; k++)
                for (int i = 0; i < numberOfControlPointsPerPatch; i++)
                    for (int j = 0; j < numberOfControlPointsPerPatch; j++)
                        matrixLStLS.coeffRef(k * numberOfControlPointsPerPatch + i, k * numberOfControlPointsPerPatch + j) += weightEnergy * m_energyMatrix.coeff(i, j);

            Eigen::VectorXd rightSideX = weightLS * matrixLS.transpose() * dataPointsX;
            Eigen::VectorXd rightSideY = weightLS * matrixLS.transpose() * dataPointsY;
            Eigen::VectorXd rightSideZ = weightLS * matrixLS.transpose() * dataPointsZ;

            const int matG0rows = matrixG0.rows();

            Eigen::VectorXd initialGuessLSG0G1X(numberOfLocalControlPoints);
            Eigen::VectorXd initialGuessLSG0G1Y(numberOfLocalControlPoints);
            Eigen::VectorXd initialGuessLSG0G1Z(numberOfLocalControlPoints);
            if (!m_initialGuess) {
                //calculate initial guess as solution of LS-fitting with G0 and simplified G1 constraints
                Eigen::SparseMatrix<double> matLSSystemG01(numberOfLocalControlPoints + 2 * matG0rows, numberOfLocalControlPoints + 2 * matG0rows);
                for (int i = 0; i < numberOfLocalControlPoints; i++)
                    for (int j = 0; j < numberOfLocalControlPoints; j++)
                        matLSSystemG01.insert(i, j) = matrixLStLS.coeff(i, j);

                for (int i = 0; i < matG0rows; i++) {
                    for (int j = 0; j < numberOfLocalControlPoints; j++) {
                        matLSSystemG01.insert(numberOfLocalControlPoints + i, j) = matrixG0.coeff(i, j);
                        matLSSystemG01.insert(j, numberOfLocalControlPoints + i) = matrixG0.coeff(i, j);
                        matLSSystemG01.insert(numberOfLocalControlPoints + matG0rows + i, j) = matrixG1.coeff(i, j);
                        matLSSystemG01.insert(j, numberOfLocalControlPoints + matG0rows + i) = matrixG1.coeff(i, j);
                    }
                }
                Eigen::VectorXd rightSideLSSystemG01X(numberOfLocalControlPoints + 2 * matG0rows);
                Eigen::VectorXd rightSideLSSystemG01Y(numberOfLocalControlPoints + 2 * matG0rows);
                Eigen::VectorXd rightSideLSSystemG01Z(numberOfLocalControlPoints + 2 * matG0rows);
                for (int i = 0; i < numberOfLocalControlPoints; i++) {
                    rightSideLSSystemG01X[i] = rightSideX[i];
                    rightSideLSSystemG01Y[i] = rightSideY[i];
                    rightSideLSSystemG01Z[i] = rightSideZ[i];
                }
                for (int i = numberOfLocalControlPoints; i < numberOfLocalControlPoints + 2 * matG0rows; i++) {
                    rightSideLSSystemG01X[i] = 0;
                    rightSideLSSystemG01Y[i] = 0;
                    rightSideLSSystemG01Z[i] = 0;
                }

                Eigen::SparseLU<Eigen::SparseMatrix<double> > initialSolver;
                initialSolver.compute(matLSSystemG01);
                initialGuessLSG0G1X = initialSolver.solve(rightSideLSSystemG01X);
                initialGuessLSG0G1Y = initialSolver.solve(rightSideLSSystemG01Y);
                initialGuessLSG0G1Z = initialSolver.solve(rightSideLSSystemG01Z);
            } else {
                for (int cfNumber = 0; cfNumber < 2; cfNumber++) {
                    for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                            const int ctrPointId = cfNumber * numberOfControlPointsPerPatch + BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
                            const QVector3D ctrPointInitialGuess = m_initialGuess->getBSplineSurface(incidentCfIds[cfNumber])->getControlPointsP()->at(i).at(j);
                            initialGuessLSG0G1X[ctrPointId] = ctrPointInitialGuess.x();
                            initialGuessLSG0G1Y[ctrPointId] = ctrPointInitialGuess.y();
                            initialGuessLSG0G1Z[ctrPointId] = ctrPointInitialGuess.z();
                        }
                    }
                }
            }

            //Create energy matrix for full system (contains individual energy matrices on its diagonals)
            Eigen::SparseMatrix<double> matrixEnergyAllPatches(numberOfLocalControlPoints, numberOfLocalControlPoints);
            for (int k = 0; k < 2; k++)
                for (int i = 0; i < numberOfControlPointsPerPatch; i++)
                    for (int j = 0; j < numberOfControlPointsPerPatch; j++)
                        matrixEnergyAllPatches.insert(k * numberOfControlPointsPerPatch + i, k * numberOfControlPointsPerPatch + j) = m_energyMatrix.coeff(i, j);

            //Calculate constant matrix and vector products here to speed up optimization
            Eigen::SparseMatrix<double> constMatGradient = weightLS * matrixLS.transpose() * matrixLS + weightEnergy * matrixEnergyAllPatches;
            Eigen::VectorXd constVecGradientX = weightLS * matrixLS.transpose() * dataPointsX;
            Eigen::VectorXd constVecGradientY = weightLS * matrixLS.transpose() * dataPointsY;
            Eigen::VectorXd constVecGradientZ = weightLS * matrixLS.transpose() * dataPointsZ;

            //Setup parameters for nlopt
            NonlinearOptParameters optParameters;
            optParameters.M = &matrixEnergyAllPatches;
            optParameters.N = &matrixLS;
            optParameters.matGradientFirstHalf = &constMatGradient;
            optParameters.dataX = &dataPointsX;
            optParameters.dataY = &dataPointsY;
            optParameters.dataZ = &dataPointsZ;
            optParameters.NDataX = &constVecGradientX;
            optParameters.NDataY = &constVecGradientY;
            optParameters.NDataZ = &constVecGradientZ;
            optParameters.weightLS = weightLS;
            optParameters.weightEnergy = weightEnergy;
            optParameters.weightNormalConstraints = m_g1EdgeConstraintsOptFactor;

            //setup optimizer
            nlopt_opt optimizer = nlopt_create(NLOPT_LD_SLSQP, 3 * numberOfLocalControlPoints);
            nlopt_set_min_objective(optimizer, nonlinearOptFunction, &optParameters);
            nlopt_set_xtol_rel(optimizer, m_nloXPrec);
            nlopt_set_ftol_rel(optimizer, m_nloFPrec);
            nlopt_set_maxeval(optimizer, m_nloMaxIterations);

            //setup G0 constraints
            QVector<G0ConstraintParameters *> g0Constraints(3 * matG0rows);
            for (int i = 0; i < matG0rows; i++) {
                int index0 = 0;
                while (matrixG0.coeff(i, index0) == 0)
                    index0++;
                int index1 = index0 + 1;
                while (matrixG0.coeff(i, index1) == 0)
                    index1++;

                for (int k = 0; k < 3; k++) {
                    G0ConstraintParameters *g0Params = new G0ConstraintParameters;
                    g0Params->index0 = index0 + (k * numberOfLocalControlPoints);
                    g0Params->index1 = index1 + (k * numberOfLocalControlPoints);
                    g0Params->coef0 = 1;
                    g0Params->coef1 = -1;
                    nlopt_add_equality_constraint(optimizer, g0Constraint, g0Params, m_nloG0constPrec);
                    g0Constraints[3 * i + k] = g0Params;
                }
            }

            //setup normal constraints (either as real constraints or as optimization values)
            QVector<NormalConstraintParameters *> normalConstraints;
            for (int i = 0; i < numberOfNormalConstraints; i++) {
                NormalConstraintParameters *params = new NormalConstraintParameters;
                params->numberOfControlPoints = numberOfLocalControlPoints;
                params->tangents = tangentCoefficients[i];
                params->id = i;

                normalConstraints << params;
                if (!m_g1EdgeConstraintsAsOpt)
                    nlopt_add_equality_constraint(optimizer, normalConstraint_normalized, params, m_nloG1constPrec);
            }
            if (m_g1EdgeConstraintsAsOpt)
                optParameters.normalConstraints = normalConstraints;

            //Setup equality constraints
            QVector<QPair<int, int> > relevantCPIndices(16);
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < 4; i++)
                    relevantCPIndices[k * 8 + i] = CP_around_V[internalCeId[k]][i];
                for (int i = 0; i < 4; i++)
                    relevantCPIndices[k * 8 + 4 + i] = CP_around_V[(internalCeId[k] + 1) % 4][i];
            }

            const int numberOfEqualConstraints = 48;
            QVector<EqualConstraintParameters *> equalConstraints(numberOfEqualConstraints);
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < 8; i++) {
                    const QPair<int, int> ij = relevantCPIndices[k * 8 + i];
                    const QVector3D &ctrPoint = globalControlPoints[incidentCfIds[k]][ij.first][ij.second];
                    const double ctrPointArray[3] = {ctrPoint.x(), ctrPoint.y(), ctrPoint.z()};
                    const int index = offset[k] + BSplineSurface::matToVecIndexLocal(ij, numberOfControlPointsPerRow);

                    for (int j = 0; j < 3; j++) {
                        EqualConstraintParameters *equalParams = new EqualConstraintParameters;
                        equalParams->index = index + j * numberOfLocalControlPoints;
                        equalParams->value = ctrPointArray[j];
                        nlopt_add_equality_constraint(optimizer, equalConstraint, equalParams, m_nloEqualConstPrec);
                        equalConstraints[(k * 8 + i) * 3 + j] = equalParams;
                    }
                }
            }

            //set initial guess
            double optVariables[3 * numberOfLocalControlPoints];
            for (int i = 0; i < numberOfLocalControlPoints; i++) {
                optVariables[i                                 ] = initialGuessLSG0G1X[i];
                optVariables[i +     numberOfLocalControlPoints] = initialGuessLSG0G1Y[i];
                optVariables[i + 2 * numberOfLocalControlPoints] = initialGuessLSG0G1Z[i];
            }

            //do optimization
            double optValue;
            int status = nlopt_optimize(optimizer, optVariables, &optValue);
            if (status < 0) {
                qDebug() << "ERROR" << status << "Value:" << optValue;
            } else {
                qDebug() << "Optimization done! Edge:" << ceId << "Status:" << status << "Value:" << optValue;
            }

            //clean up optimizer
            nlopt_destroy(optimizer);
            for (int i = 0; i < 3 * matG0rows; i++)
                delete g0Constraints[i];
            for (int i = 0; i < numberOfNormalConstraints; i++)
                delete normalConstraints[i];
            for (int i = 0; i < numberOfEqualConstraints; i++)
                delete equalConstraints[i];

            //extract the result (only the control points)
            Eigen::VectorXd ctrX(numberOfLocalControlPoints);
            Eigen::VectorXd ctrY(numberOfLocalControlPoints);
            Eigen::VectorXd ctrZ(numberOfLocalControlPoints);
            for (int i = 0; i < numberOfLocalControlPoints; i++) {
                ctrX[i] = optVariables[i];
                ctrY[i] = optVariables[i + numberOfLocalControlPoints];
                ctrZ[i] = optVariables[i + 2 * numberOfLocalControlPoints];
            }

            ///
            ///
            ///

            //Remove numerical errors in G0 constraints by using the average of the G0 points
            for (int i = 0; i < numberOfControlPointsPerRow - 4; i++) {
                int index0 = 0;
                while (matrixG0.coeff(i, index0) == 0)
                    index0++;
                int index1 = index0 + 1;
                while (matrixG0.coeff(i, index1) == 0)
                    index1++;
                const double averageX = (ctrX[index0] + ctrX[index1])/2;
                ctrX[index0] = averageX;
                ctrX[index1] = averageX;
                const double averageY = (ctrY[index0] + ctrY[index1])/2;
                ctrY[index0] = averageY;
                ctrY[index1] = averageY;
                const double averageZ = (ctrZ[index0] + ctrZ[index1])/2;
                ctrZ[index0] = averageZ;
                ctrZ[index1] = averageZ;
            }

            //Store the G1-edge-relevant control points in global control point matrix
            for (int cfNumber = 0; cfNumber < 2; cfNumber++) {
                for (int i = 0; i < 2 * (numberOfControlPointsPerRow - 4); i++) {
                    const int ctrPointIndex = importantCtrPointIndices[cfNumber][i];
                    const int localIndex = ctrPointIndex - (cfNumber * numberOfControlPointsPerPatch);
                    int index_i, index_j;
                    BSplineSurface::vecToMatIndexLocal(localIndex, index_i, index_j, numberOfControlPointsPerRow);
                    const QVector3D ctrPoint = QVector3D(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);

                    globalControlPoints[incidentCellFaces[cfNumber]->getId()][index_i][index_j] = ctrPoint;
                    computedControlPoints[incidentCellFaces[cfNumber]->getId()][index_i][index_j] = true;
                    debug_edgeFitPoints << ctrPoint;
                }
            }

            //==============//
            // OUTPUT Start //
            //==============//
            const double errorX = (matrixLS * ctrX - dataPointsX).norm();
            const double errorY = (matrixLS * ctrY - dataPointsY).norm();
            const double errorZ = (matrixLS * ctrZ - dataPointsZ).norm();
            const double lsError = errorX * errorX + errorY * errorY + errorZ * errorZ;
            const double energyX = ctrX.transpose() * matrixEnergyAllPatches * ctrX;
            const double energyY = ctrY.transpose() * matrixEnergyAllPatches * ctrY;
            const double energyZ = ctrZ.transpose() * matrixEnergyAllPatches * ctrZ;
            std::cout << "Squared error: " << lsError << " | per datapoint: " << lsError / numberOfLocalDataPoints << std::endl;
            std::cout << "Energy: " << energyX + energyY + energyZ << " x, y, z: " << energyX << " " << energyY << " " << energyZ << std::endl;


            //G0 errors will be 0 (due to postprocessing)
            //std::cout << "G0 errors [x]: " << (matrixG0 * ctrX).transpose() << " [y]: " << (matrixG0 * ctrY).transpose() << " [z]: " << (matrixG0 * ctrZ).transpose() << std::endl;

            QVector<double> normalErrors(numberOfNormalConstraints);
            double maxNormalError = -1;
            double averageNormalError = 0;
            for (int i = 0; i < numberOfNormalConstraints; i++) {
                const QVector3D Tu0(tangentCoefficients.at(i).Nu0.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nu0.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nu0.transpose() * ctrZ);
                const QVector3D Tv0(tangentCoefficients.at(i).Nv0.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nv0.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nv0.transpose() * ctrZ);
                const QVector3D Tu1(tangentCoefficients.at(i).Nu1.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nu1.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nu1.transpose() * ctrZ);
                const QVector3D Tv1(tangentCoefficients.at(i).Nv1.transpose() * ctrX,
                                    tangentCoefficients.at(i).Nv1.transpose() * ctrY,
                                    tangentCoefficients.at(i).Nv1.transpose() * ctrZ);

                const double normalError = (QVector3D::crossProduct(Tu0, Tv0).normalized() - QVector3D::crossProduct(Tu1, Tv1).normalized()).length();

                normalErrors[i] = normalError;
                averageNormalError += normalError;
                if (normalError > maxNormalError)
                    maxNormalError = normalError;
            }
            averageNormalError /= (double) numberOfNormalConstraints;
            std::cout << "Normal errors (" << numberOfNormalConstraints << " constraints): ";
            if (m_detailedOutput) {
                for (int i = 0; i < numberOfNormalConstraints; i++) {
                    //std::cout << "[" << normalError.x() << " " << normalError.y() << " " << normalError.z() << "](" << normalError.length() << ") ";
                    std::cout << "(" << normalErrors[i] << ") ";
                }
            }
            std::cout << "Average: " << averageNormalError << " Max: " << maxNormalError << std::endl;

//            std::cout << "Equality errors: ";
//            for (int i = 0; i < equalConstraints.size(); i++) {
//                double diff = optVariables[equalConstraints[i].first] - equalConstraints[i].second;
//                std::cout << diff << " ";
//            }
//            std::cout << std::endl;;
            //============//
            // OUTPUT END //
            //============//

            ///=============///
            /// DEBUG START ///
            ///=============///
//            if (m_debugingMode >= 2) {
//                bSplineSurfaces.clear();
//                for (int cfNumber = 0; cfNumber < 2; cfNumber++) {
//                    QVector<QVector<QVector3D> > patchControlPoints(m_numberOfControlPointsPerRow, QVector<QVector3D>(m_numberOfControlPointsPerRow));
//                    for (int i = 0; i < m_numberOfControlPointsPerRow; i++) {
//                        for (int j = 0; j < m_numberOfControlPointsPerRow; j++) {
//                            const int ctrPointIndex = cfNumber * m_numberOfControlPointsPerPatch + matToVecIndexLocal(i, j, m_numberOfControlPointsPerRow);
//                            patchControlPoints[i][j] = QVector3D(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);
//                        }
//                    }

//                    BSplineSurface *splinePatch = new BSplineSurface();
//                    splinePatch->setKnotsU(knots);
//                    splinePatch->setKnotsV(knots);
//                    splinePatch->setOrder(m_splineOrder);
//                    splinePatch->setControlPoints(patchControlPoints);
//                    bSplineSurfaces.append(splinePatch);
//                }

//        //        Vertex v0(m_bSplineSurfaces[0]->evaluate(0.8, 1), -1);
//        //        Vertex v1(m_bSplineSurfaces[1]->evaluate(0, 0.2), -1);
//        //        Vertex tv0(m_bSplineSurfaces[0]->evaluate(0.8, 1) + m_bSplineSurfaces[0]->evaluateTv(0.8, 1), -1);
//        //        Vertex tv1(m_bSplineSurfaces[1]->evaluate(0, 0.2) + m_bSplineSurfaces[1]->evaluateTv(0, 0.2), -1);
//        //        Vertex tu0(m_bSplineSurfaces[0]->evaluate(0.8, 1) + m_bSplineSurfaces[0]->evaluateTu(0.8, 1), -1);
//        //        Vertex tu1(m_bSplineSurfaces[1]->evaluate(0, 0.2) + m_bSplineSurfaces[1]->evaluateTu(0, 0.2), -1);

//        //        QVector<Edge *> tangents0;
//        //        tangents0 << new Edge(&v0, &tv0, -1) << new Edge(&v0, &tu0, -1);
//        //        QVector<Edge *> tangents1;
//        //        tangents1 << new Edge(&v1, &tv1, -1) << new Edge(&v1, &tu1, -1);
//        //        m_stlWriter->writeEdgesToStl(&tangents0, QString("../../tangents0.stl"),0.3);
//        //        m_stlWriter->writeEdgesToStl(&tangents1, QString("../../tangents1.stl"),0.3);

//                QTime t;
//                t.start();
//                this->bSplinesToStl(51, QString("edgefit_%1").arg(ceId), &t);
//                bSplineSurfaces.clear();
//                m_bSplineCurves.clear();
//            }
            ///===========///
            /// DEBUG END ///
            ///===========///
        } else {
            qDebug() << "no G1 edge -> need to do G0 fit!";
            //TODO need to do G0 fit
        }

        edgeFittingDone[ceId] = true;
        int finishedEdgeFits = 0;
        for (int i = 0; i < numberOfCellEdges; i++)
            if (edgeFittingDone[i])
                finishedEdgeFits++;
        std::cout << "Edge fits done: " << finishedEdgeFits << "/" << numberOfCellEdges << std::endl;
    }

    //const int numberOfInnerCPs = (numberOfControlPointsPerRow - 4) * (numberOfControlPointsPerRow - 4);

    //Final step: Do LS fitting for the rest
    #pragma omp parallel for
    for (int cfId = 0; cfId < numberOfCellFaces; cfId++) {
        qDebug() << "CellFace" << cfId;

        CellFace *cf = m_cellMesh->getFace(cfId);

        //count the total number of datapoints in cellfaces around the cellvertex
        const int numberOfLocalDataPoints = cf->getMeshVertices()->size();

        Parameterization *cfParam = m_cellMesh->getQuadCellParameterization(cfId);
        QVector<Vertex *> *data = cf->getMeshVertices();

        //Construct coefficient matrix and data point vector for approximation (least square system)
        Eigen::SparseMatrix<double, Eigen::RowMajor> matrixLS(numberOfLocalDataPoints, numberOfControlPointsPerPatch);
        const int nonZeroBasisFunctions = m_properties->getOrder() + 1;
        matrixLS.reserve(Eigen::VectorXi::Constant(numberOfLocalDataPoints, nonZeroBasisFunctions * nonZeroBasisFunctions));
        Eigen::VectorXd dataPointsX(numberOfLocalDataPoints);
        Eigen::VectorXd dataPointsY(numberOfLocalDataPoints);
        Eigen::VectorXd dataPointsZ(numberOfLocalDataPoints);

        for (int k = 0; k < numberOfLocalDataPoints; k++) {
            Vertex *v = data->at(k);
            QVector2D vParam = cfParam->getParameter(v);

            dataPointsX[k] = v->getPosition().x();
            dataPointsY[k] = v->getPosition().y();
            dataPointsZ[k] = v->getPosition().z();

            for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                const double nu = m_properties->N(vParam.x(), i);
                for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                    const double nv = m_properties->N(vParam.y(), j);
                    const double product = nu*nv;
                    if (product != 0)
                        matrixLS.insert(k, BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow)) = product;
                }
            }
        }

        double weightLS = m_weightLSTerm;
        double weightEnergy = m_weightEnergyTerm;
        if (!m_useFixedWeights) {
            weightLS = (double) 1 / numberOfLocalDataPoints;
            weightEnergy = 1 / cf->calculateAreaByCorners();
        }

        Eigen::VectorXd initialGuessX(numberOfControlPointsPerPatch);
        Eigen::VectorXd initialGuessY(numberOfControlPointsPerPatch);
        Eigen::VectorXd initialGuessZ(numberOfControlPointsPerPatch);
        if (!m_initialGuess) {
            //calculate initial guess
            Eigen::SparseMatrix<double> matrixLStLS = weightLS * matrixLS.transpose() * matrixLS + weightEnergy * m_energyMatrix;

            Eigen::VectorXd rightSideX = weightLS * matrixLS.transpose() * dataPointsX;
            Eigen::VectorXd rightSideY = weightLS * matrixLS.transpose() * dataPointsY;
            Eigen::VectorXd rightSideZ = weightLS * matrixLS.transpose() * dataPointsZ;

            Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
            solver.compute(matrixLStLS);
            initialGuessX = solver.solve(rightSideX);
            initialGuessY = solver.solve(rightSideY);
            initialGuessZ = solver.solve(rightSideZ);
        } else {
            for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                    const int ctrPointId = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
                    const QVector3D ctrPointInitialGuess = m_initialGuess->getBSplineSurface(cfId)->getControlPointsP()->at(i).at(j);
                    initialGuessX[ctrPointId] = ctrPointInitialGuess.x();
                    initialGuessY[ctrPointId] = ctrPointInitialGuess.y();
                    initialGuessZ[ctrPointId] = ctrPointInitialGuess.z();
                }
            }
        }

        //Calculate constant matrix and vector products here to speed up optimization
        Eigen::SparseMatrix<double> constMatGradient = weightLS * matrixLS.transpose() * matrixLS + weightEnergy * m_energyMatrix;
        Eigen::VectorXd constVecGradientX = weightLS * matrixLS.transpose() * dataPointsX;
        Eigen::VectorXd constVecGradientY = weightLS * matrixLS.transpose() * dataPointsY;
        Eigen::VectorXd constVecGradientZ = weightLS * matrixLS.transpose() * dataPointsZ;

        //Setup parameters for nlopt
        NonlinearOptParameters optParameters;
        optParameters.M = &m_energyMatrix;
        optParameters.N = &matrixLS;
        optParameters.matGradientFirstHalf = &constMatGradient;
        optParameters.dataX = &dataPointsX;
        optParameters.dataY = &dataPointsY;
        optParameters.dataZ = &dataPointsZ;
        optParameters.NDataX = &constVecGradientX;
        optParameters.NDataY = &constVecGradientY;
        optParameters.NDataZ = &constVecGradientZ;
        optParameters.weightLS = weightLS;
        optParameters.weightEnergy = weightEnergy;

        //setup optimizer
        nlopt_opt optimizer = nlopt_create(NLOPT_LD_SLSQP, 3 * numberOfControlPointsPerPatch);
        nlopt_set_min_objective(optimizer, nonlinearOptFunction, &optParameters);
        nlopt_set_xtol_rel(optimizer, m_nloXPrec);
        nlopt_set_ftol_rel(optimizer, m_nloFPrec);
        nlopt_set_maxeval(optimizer, m_nloMaxIterations);

        //setup equality constraints
        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
            for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                if (computedControlPoints[cfId][i][j]) {
                    QVector3D &ctrPoint = globalControlPoints[cfId][i][j];
                    const int index = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
                    EqualConstraintParameters *paramsX = new EqualConstraintParameters;
                    paramsX->index = index;
                    paramsX->value = ctrPoint.x();
                    nlopt_add_equality_constraint(optimizer, equalConstraint, paramsX, m_nloEqualConstPrec);

                    EqualConstraintParameters *paramsY = new EqualConstraintParameters;
                    paramsY->index = index + numberOfControlPointsPerPatch;
                    paramsY->value = ctrPoint.y();
                    nlopt_add_equality_constraint(optimizer, equalConstraint, paramsY, m_nloEqualConstPrec);

                    EqualConstraintParameters *paramsZ = new EqualConstraintParameters;
                    paramsZ->index = index + 2 * numberOfControlPointsPerPatch;
                    paramsZ->value = ctrPoint.z();
                    nlopt_add_equality_constraint(optimizer, equalConstraint, paramsZ, m_nloEqualConstPrec);
                }
            }
        }

        //set initial guess
        double optVariables[3 * numberOfControlPointsPerPatch];
        for (int i = 0; i < numberOfControlPointsPerPatch; i++) {
            optVariables[i] = initialGuessX[i];
            optVariables[i + numberOfControlPointsPerPatch] = initialGuessY[i];
            optVariables[i + 2 * numberOfControlPointsPerPatch] = initialGuessZ[i];
        }

        //do optimization
        double optValue;
        int status = nlopt_optimize(optimizer, optVariables, &optValue);
        if (status < 0) {
            qDebug() << "ERROR" << status << "Value:" << optValue;
        } else {
            qDebug() << "Optimization done! Face:" << cfId << "Status:" << status << "Value:" << optValue;
        }
        nlopt_destroy(optimizer);

        //extract the result (only the control points)
        Eigen::VectorXd ctrX(numberOfControlPointsPerPatch);
        Eigen::VectorXd ctrY(numberOfControlPointsPerPatch);
        Eigen::VectorXd ctrZ(numberOfControlPointsPerPatch);
        for (int i = 0; i < numberOfControlPointsPerPatch; i++) {
            ctrX[i] = optVariables[i];
            ctrY[i] = optVariables[i + numberOfControlPointsPerPatch];
            ctrZ[i] = optVariables[i + 2 * numberOfControlPointsPerPatch];
        }

        //Store the control points
        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
            for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                if (!computedControlPoints[cfId][i][j]) {
                    const int ctrPointIndex = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
                    const QVector3D ctrPoint(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);

                    globalControlPoints[cfId][i][j] = ctrPoint;
                    debug_innerFitPoints << ctrPoint;
                }
            }
        }

        //==============//
        // OUTPUT Start //
        //==============//
        const double errorX = (matrixLS * ctrX - dataPointsX).norm();
        const double errorY = (matrixLS * ctrY - dataPointsY).norm();
        const double errorZ = (matrixLS * ctrZ - dataPointsZ).norm();
        const double lsError = errorX * errorX + errorY * errorY + errorZ * errorZ;
        const double energyX = ctrX.transpose() * m_energyMatrix * ctrX;
        const double energyY = ctrY.transpose() * m_energyMatrix * ctrY;
        const double energyZ = ctrZ.transpose() * m_energyMatrix * ctrZ;
        std::cout << "Squared error: " << lsError << " | per datapoint: " << lsError / numberOfLocalDataPoints << std::endl;
        std::cout << "Energy: " << energyX + energyY + energyZ << " x, y, z: " << energyX << " " << energyY << " " << energyZ << std::endl;
        //============//
        // OUTPUT END //
        //============//

        faceFittingDone[cfId] = true;
        int finishedFaceFits = 0;
        for (int i = 0; i < numberOfCellFaces; i++)
            if (faceFittingDone[i])
                finishedFaceFits++;
        std::cout << "Face fits done: " << finishedFaceFits << "/" << numberOfCellFaces << std::endl;
    }

//    m_stlWriter->writePoints(&debug_vertexFitPoints, 0.4, QString("../out/debug_vertexfit_points.stl"));
//    m_stlWriter->writePoints(&debug_edgeFitPoints, 0.4, QString("../out/debug_edgefit_points.stl"));
//    m_stlWriter->writePoints(&debug_innerFitPoints, 0.4, QString("../out/debug_innerfit_points.stl"));

    //Setup B-Splines
    bSplineSurfaces.clear();
    for (int i = 0; i < numberOfCellFaces; i++) {
        BSplineSurface *splinePatch = new BSplineSurface(m_properties, m_properties, globalControlPoints[i]);
        bSplineSurfaces.append(splinePatch);
    }

    BSplinePatchNetwork *patchNetwork = new BSplinePatchNetwork(m_cellMesh, m_properties, bSplineSurfaces);
    if (m_useFixedWeights)
        patchNetwork->addInfo(QString("B-spline fitting done with weights %1 LS %2 Energy.\n").arg(m_weightLSTerm).arg(m_weightEnergyTerm));
    else
        patchNetwork->addInfo(QString("B-spline fitting done with adaptive weights\n"));

    return patchNetwork;
}

BSplinePatchNetwork *Splinefitter::doApproximationGlobalC0G1Constraints(bool onlyC0)
{
    qDebug() << "global method probably needs fixing!";

    const int numberOfCells = m_cellMesh->getNumberOfFaces();
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();
    const int numberOfControlPointsPerPatch = numberOfControlPointsPerRow * numberOfControlPointsPerRow;
    const int totalNumberOfControlPoints = numberOfControlPointsPerPatch * numberOfCells;

    int totalNumberOfDataPoints = 0;
    for (int fId = 0; fId < numberOfCells; fId++) {
        CellFace *cf = m_cellMesh->getFace(fId);
        totalNumberOfDataPoints += cf->getMeshVertices()->size();
        for (int i = 0; i < 4; i++) //TODO once, double endpoints are removed, consider reduce the total size accordingly
            totalNumberOfDataPoints += ((CellEdge *) cf->getEdge(i))->getPolylineVertices()->size();
    }

    //Construct coefficient matrix for data point approximation (least square system)
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrixLS(totalNumberOfDataPoints, totalNumberOfControlPoints);
    const int nonZeroBasisFunctions = m_properties->getOrder() + 1;
    matrixLS.reserve(Eigen::VectorXi::Constant(totalNumberOfDataPoints, nonZeroBasisFunctions * nonZeroBasisFunctions));
    Eigen::VectorXd dataPointsX(totalNumberOfDataPoints);
    Eigen::VectorXd dataPointsY(totalNumberOfDataPoints);
    Eigen::VectorXd dataPointsZ(totalNumberOfDataPoints);

    int lineOffset = 0;
    for (int fId = 0; fId < numberOfCells; fId++) {
        CellFace *cf = m_cellMesh->getFace(fId);
        //Parameterization *cfParam = Parameterization::harmonicMap(cf, m_cellMesh->getOriginalMesh(), Parameterization::WT_length, Parameterization::BT_equidistantUnitSquare);
        Parameterization *cfParam = m_cellMesh->getQuadCellParameterization(fId);
        QVector<Vertex *> data = *cf->getMeshVertices();
        for (int eId = 0; eId < 4; eId++) {
            QVector<PolylineVertex *> *polyVertices = ((CellEdge *) cf->getEdge(eId))->getPolylineVertices();
            for (int k = 0; k < polyVertices->size(); k++) {
                data.append(polyVertices->at(k));
                //TODO remove double endpoints
            }
        }
        const int numberOfDataPoints = data.size();
        const int columnOffset = numberOfControlPointsPerPatch * fId;

        for (int k = 0; k < numberOfDataPoints; k++) {
            Vertex *v = data.at(k);
            QVector2D vParam = cfParam->getParameter(v);

            dataPointsX[lineOffset + k] = v->getPosition().x();
            dataPointsY[lineOffset + k] = v->getPosition().y();
            dataPointsZ[lineOffset + k] = v->getPosition().z();

            for (int i = 0; i < numberOfControlPointsPerRow; i++) {
                const double nu = m_properties->N(vParam.x(), i);
                for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                    const double nv = m_properties->N(vParam.y(), j);
                    const double product = nu * nv;
                    if (product != 0)
                        matrixLS.insert(lineOffset + k, columnOffset + (i * numberOfControlPointsPerRow) + j) = product;
                }
            }
        }

        lineOffset += numberOfDataPoints;
    }

    //construct matrices for continuity constraints
    int numberOfG1Edges = 0;
    for (int i = 0; i < numberOfCellEdges; i++)
        if (m_cellMesh->getEdge(i)->isG1Edge())
            numberOfG1Edges++;

    if (onlyC0)
        numberOfG1Edges = 0;

    Eigen::SparseMatrix<double> matrixG0(numberOfControlPointsPerRow * numberOfCellEdges, totalNumberOfControlPoints);
    Eigen::SparseMatrix<double> matrixG1(numberOfControlPointsPerRow * numberOfG1Edges, totalNumberOfControlPoints);
    int iG1 = 0;
    for (int i = 0; i < numberOfCellEdges; i++) {
        CellEdge *ce = m_cellMesh->getEdge(i);

        const int fId[2] = {ce->getIncidentFaceId(0), ce->getIncidentFaceId(1)};
        CellFace *cf[2] = {m_cellMesh->getFace(fId[0]), m_cellMesh->getFace(fId[1])};
        const int internalCeId[2] = {cf[0]->getInternalEdgeId(ce), cf[1]->getInternalEdgeId(ce)};
        const int offset[2] = {fId[0] * numberOfControlPointsPerPatch, fId[1] * numberOfControlPointsPerPatch};

        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
            int globalCPIndex[2];
            int globalCPInnerIndex[2];
            for (int k = 0; k < 2; k++) {
                int controlPointCounter = j;
                if (cf[k]->edgeIsInverted(internalCeId[k]))
                    controlPointCounter = numberOfControlPointsPerRow - j - 1;
                if (internalCeId[k] == 0) {
                    globalCPIndex[k] = offset[k] + (controlPointCounter * numberOfControlPointsPerRow);
                    globalCPInnerIndex[k] = offset[k] + (controlPointCounter * numberOfControlPointsPerRow) + 1;
                } else if (internalCeId[k] == 1) {
                    globalCPIndex[k] = offset[k] + (numberOfControlPointsPerRow * (numberOfControlPointsPerRow - 1)) + controlPointCounter;
                    globalCPInnerIndex[k] = offset[k] + (numberOfControlPointsPerRow * (numberOfControlPointsPerRow - 2)) + controlPointCounter;
                } else if (internalCeId[k] == 2) {
                    globalCPIndex[k] = offset[k] + numberOfControlPointsPerRow - 1 + (numberOfControlPointsPerRow * (numberOfControlPointsPerRow - 1 - controlPointCounter));
                    globalCPInnerIndex[k] = offset[k] + numberOfControlPointsPerRow - 2 + (numberOfControlPointsPerRow * (numberOfControlPointsPerRow - 1 - controlPointCounter));
                } else if (internalCeId[k] == 3) {
                    globalCPIndex[k] = offset[k] + (numberOfControlPointsPerRow - 1 - controlPointCounter);
                    globalCPInnerIndex[k] = offset[k] + numberOfControlPointsPerRow + (numberOfControlPointsPerRow - 1 - controlPointCounter);
                } else {
                    qDebug() << "ERROR in B-spline approximation: Internal Index of edge is not in [0,3] ->" << internalCeId[k];
                }
            }

            matrixG0.insert(i * numberOfControlPointsPerRow + j, globalCPIndex[0]) = 1;
            matrixG0.insert(i * numberOfControlPointsPerRow + j, globalCPIndex[1]) = -1;

            if (ce->isG1Edge() && !onlyC0) {
                matrixG1.insert(iG1 * numberOfControlPointsPerRow + j, globalCPIndex[0]) = -1;
                matrixG1.insert(iG1 * numberOfControlPointsPerRow + j, globalCPInnerIndex[0]) = 1;
                matrixG1.insert(iG1 * numberOfControlPointsPerRow + j, globalCPIndex[1]) = -1;
                matrixG1.insert(iG1 * numberOfControlPointsPerRow + j, globalCPInnerIndex[1]) = 1;
            }
        }
        if (ce->isG1Edge())
            iG1++;
    }

    //Solve least squares with side constraints
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrixToInvert = m_weightLSTerm * matrixLS.transpose() * matrixLS;
    matrixToInvert.resize(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow,
                          totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow);

    Eigen::VectorXi nonZeros = Eigen::VectorXi::Constant(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow, totalNumberOfControlPoints + 6);
    for (int i = totalNumberOfControlPoints; i < totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow; i++)
        nonZeros[i] = 6;
    matrixToInvert.reserve(nonZeros);

    //Normalize
//    matrixLStLS *= 1.0/(double) totalNumberOfDataPoints;
//    dataPointsX *= 1.0/(double) totalNumberOfDataPoints;
//    dataPointsY *= 1.0/(double) totalNumberOfDataPoints;
//    dataPointsZ *= 1.0/(double) totalNumberOfDataPoints;
//    matrixEnergyGradCoef *= 1.0/m_isosurface->calculateTotalSurfaceArea();

    for (int k = 0; k < numberOfCells; k++) {
        for (int i = 0; i < numberOfControlPointsPerPatch; i++) {
            for (int j = 0; j < numberOfControlPointsPerPatch; j++) {
                matrixToInvert.coeffRef(k * numberOfControlPointsPerPatch + i, k * numberOfControlPointsPerPatch + j) += m_weightEnergyTerm * m_energyMatrix.coeff(i, j);
            }
        }
    }

    for (int i = 0; i < numberOfCellEdges * numberOfControlPointsPerRow; i++) {
        for (int j = 0; j < totalNumberOfControlPoints; j++) {
            const double coef = matrixG0.coeff(i, j);
            if (coef) {
                matrixToInvert.insert(totalNumberOfControlPoints + i, j) = coef;
                matrixToInvert.insert(j, totalNumberOfControlPoints + i) = coef;
            }
        }
    }

    for (int i = 0; i < numberOfG1Edges * numberOfControlPointsPerRow; i++) {
        for (int j = 0; j < totalNumberOfControlPoints; j++) {
            const double coef = matrixG1.coeff(i, j);
            if (coef) {
                matrixToInvert.insert(totalNumberOfControlPoints + numberOfCellEdges * numberOfControlPointsPerRow + i, j) = coef;
                matrixToInvert.insert(j, totalNumberOfControlPoints + numberOfCellEdges * numberOfControlPointsPerRow + i) = coef;
            }
        }
    }

    Eigen::VectorXd rightSideX(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow);
    Eigen::VectorXd rightSideY(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow);
    Eigen::VectorXd rightSideZ(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow);

    Eigen::VectorXd tmpX = matrixLS.transpose() * dataPointsX;
    Eigen::VectorXd tmpY = matrixLS.transpose() * dataPointsY;
    Eigen::VectorXd tmpZ = matrixLS.transpose() * dataPointsZ;
    for (int i = 0; i < totalNumberOfControlPoints; i++) {
        rightSideX[i] = tmpX[i];
        rightSideY[i] = tmpY[i];
        rightSideZ[i] = tmpZ[i];
    }
    for (int i = totalNumberOfControlPoints; i <  totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow; i++) {
        rightSideX[i] = 0;
        rightSideY[i] = 0;
        rightSideZ[i] = 0;
    }

    rightSideX *= m_weightLSTerm;
    rightSideY *= m_weightLSTerm;
    rightSideZ *= m_weightLSTerm;

    //TODO test other solvers (i.e. compare running time and accuracy)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > inverter;
    inverter.compute(matrixToInvert);
//    Eigen::SparseMatrix<double> matInverted(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow,
//                                            totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow);

//    for (int i = 0; i < totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow; i++) {
//        qDebug() << i << totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow;
//        Eigen::VectorXd tmpVec = inverter.solve(Eigen::VectorXd::Unit(totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow, i));
//        for (int j = 0; j < totalNumberOfControlPoints + (numberOfCellEdges + numberOfG1Edges) * numberOfControlPointsPerRow; j++) {
//            if (tmpVec[j] != 0)
//                matInverted.insert(i, j) = tmpVec[j];
//        }
//    }

//    //TODO Bx, By, Bz also contain the lagrange multipliers
//    Eigen::VectorXd Bx = matInverted * rightSideX;
//    Eigen::VectorXd By = matInverted * rightSideY;
//    Eigen::VectorXd Bz = matInverted * rightSideZ;
    Eigen::VectorXd Bx = inverter.solve(rightSideX);
    Eigen::VectorXd By = inverter.solve(rightSideY);
    Eigen::VectorXd Bz = inverter.solve(rightSideZ);

    //Extract control points
    QVector<BSplineSurface *> bSplineSurfaces(numberOfCells);
    for (int fId = 0; fId < numberOfCells; fId++) {
        const int offset = (numberOfControlPointsPerRow * numberOfControlPointsPerRow) * fId;

        QVector<QVector<QVector3D> > controlPoints(numberOfControlPointsPerRow, QVector<QVector3D>(numberOfControlPointsPerRow));
        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
            for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                controlPoints[i][j] = QVector3D(Bx[offset + (i * numberOfControlPointsPerRow) + j],
                                                By[offset + (i * numberOfControlPointsPerRow) + j],
                                                Bz[offset + (i * numberOfControlPointsPerRow) + j]);
            }
        }

        bSplineSurfaces[fId] = new BSplineSurface(m_properties, m_properties, controlPoints);
    }

    return new BSplinePatchNetwork(m_cellMesh, m_properties, bSplineSurfaces);
}

#include "energymatrixgenerator.h"
BSplinePatchNetwork *Splinefitter::doApproximationLocalC0WithBoundaryInterpolation(bool minimizeEnergyOnBoundaries, double weightLSBoundary, double weightEnergyBoundary)
{
    qDebug() << "move energy matrix generation outside of the spline fitter";

    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();
    Eigen::SparseMatrix<double> curveEnergyMatrix(numberOfControlPointsPerRow, numberOfControlPointsPerRow);
    if (minimizeEnergyOnBoundaries)
        curveEnergyMatrix = EnergyMatrixGenerator().getEnergyMatrixCurve(m_properties, 2,
                                                                                                m_properties->getInnerKnotsStart(),
                                                                                                m_properties->getInnerKnotsEnd());
    //Build curve network
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    QVector<BSplineCurve *> bSplineCurves(numberOfCellEdges);

    #pragma omp parallel for
    for (int ceId = 0; ceId < numberOfCellEdges; ceId++)
        bSplineCurves[ceId] = this->fitSingleCurve_endPointInterpol(weightLSBoundary, weightEnergyBoundary,
                                                                    m_cellMesh->getEdge(ceId)->getPolylineVertexPositions(),
                                                                    m_cellMesh->getBoundaryParameterizationFromIncidentQuadParam(ceId),
                                                                    curveEnergyMatrix);

    //Build surface patches
    const int numberOfCells = m_cellMesh->getNumberOfFaces();
    QVector<BSplineSurface *> bSplineSurfaces(numberOfCells);

//    qDebug() << "manual energy weights";
    #pragma omp parallel for
    for (int cfId = 0; cfId < numberOfCells; cfId++) {
//        if (cfId <= 1)
//            this->setWeightEnergyTerm(0.99, true);
//        else
//            this->setWeightEnergyTerm(0.05, true);
        CellFace *cf = m_cellMesh->getFace(cfId);
        QVector<BSplineCurve *> boundaryCurves(4);
        for (int i = 0; i < 4; i++)
            boundaryCurves[i] = bSplineCurves[cf->getEdge(i)->getId()];

        bSplineSurfaces[cfId] = this->fitSingleSurface_boundaryInterpolation(cfId, boundaryCurves, m_weightLSTerm, m_weightEnergyTerm);
    }

    return new BSplinePatchNetwork(m_cellMesh, m_properties, bSplineSurfaces, bSplineCurves);
}

BSplinePatchNetwork *Splinefitter::doIterativeCurveAndSurfaceFittingC0BoundaryInterpolation(const QVector<double> &energyWeightsCurve, const QVector<double> &energyWeightsSurf, const bool useNonlinearSolverForOpt)
{
    qDebug() << "move energy matrix generation outside of the spline fitter";
    Eigen::SparseMatrix<double> curveEnergyMatrix = EnergyMatrixGenerator().getEnergyMatrixCurve(m_properties, 2,
                                                                                                                        m_properties->getInnerKnotsStart(),
                                                                                                                        m_properties->getInnerKnotsEnd());

    const int numberOfCurveFittings = energyWeightsCurve.size();
    const int numberOfCellEdges = m_cellMesh->getNumberOfEdges();
    QVector<BSplineCurve *> bSplineCurves(numberOfCellEdges);

    //Build curve network
    #pragma omp parallel for
    for (int ceId = 0; ceId < numberOfCellEdges; ceId++) {
        QVector<double> ceParam = m_cellMesh->getBoundaryParameterizationFromIncidentQuadParam(ceId);
        QVector<QVector3D> dataPoints = m_cellMesh->getEdge(ceId)->getPolylineVertexPositions();
        for (int i = 0; i < numberOfCurveFittings; i++) {
            const double wEnergy = energyWeightsCurve[i];
            const double wLS = 1 - wEnergy;
            bSplineCurves[ceId] = this->fitSingleCurve_endPointInterpol(wLS, wEnergy, dataPoints, ceParam, curveEnergyMatrix);
            ceParam = ParameterOptimizer::optimizeParameterizationSingleCurve(bSplineCurves[ceId], dataPoints, ceParam, true, useNonlinearSolverForOpt);
        }
        m_cellMesh->updateBoundaryParameterization(ceId, ceParam, true);
    }
    qDebug() << "Curves done";

    //Build surface patches
    const int numberOfSurfaceFittings = energyWeightsSurf.size();
    const int numberOfCells = m_cellMesh->getNumberOfFaces();
    QVector<BSplineSurface *> bSplineSurfaces(numberOfCells);

    for (int i = 0; i < numberOfSurfaceFittings; i++) {
        const double wEnergy = energyWeightsSurf[i];
        const double wLS = 1 - wEnergy;
        qDebug() << wEnergy << wLS;

        //TODO check this! (running it in parallel causes segmentation faults)
        //#pragma omp parallel for
        for (int cfId = 0; cfId < numberOfCells; cfId++) {
            CellFace *cf = m_cellMesh->getFace(cfId);
            QVector<BSplineCurve *> boundaryCurves(4);
            for (int i = 0; i < 4; i++)
                boundaryCurves[i] = bSplineCurves[cf->getEdge(i)->getId()];

            bSplineSurfaces[cfId] = this->fitSingleSurface_boundaryInterpolation(cfId, boundaryCurves, wLS, wEnergy);
        }

        //Not pretty, but destructor destroys surfaces and curves as well
        qDebug() << "pointers to patch notwork are still there";
        BSplinePatchNetwork *tmpSurfaceNetwork = new BSplinePatchNetwork(m_cellMesh, m_properties, bSplineSurfaces, bSplineCurves);

        const double optValue = ParameterOptimizer::optimizeParameterizationFullSurfaceNetwork(tmpSurfaceNetwork, false, useNonlinearSolverForOpt);

        qDebug() << "Surfaces optimized:" << optValue;
    }

    return new BSplinePatchNetwork(m_cellMesh, m_properties, bSplineSurfaces, bSplineCurves);
}

BSplineCurve *Splinefitter::fitSingleCurve_endPointInterpol(const double weightLS, const double weightEnergy, const QVector<QVector3D> &dataPoints, const QVector<double> &parameterization, Eigen::SparseMatrix<double> &curveEnergyMatrix)
{
    const int numberOfDataPoints = dataPoints.size();
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();

    //Set up least squares matrix
    Eigen::SparseMatrix<double> N(numberOfDataPoints, numberOfControlPointsPerRow);
    Eigen::MatrixXd dataPointMatrix(numberOfDataPoints, 3);

    for (int j = 0; j < numberOfDataPoints; j++) {
        dataPointMatrix(j, 0) = dataPoints[j].x();
        dataPointMatrix(j, 1) = dataPoints[j].y();
        dataPointMatrix(j, 2) = dataPoints[j].z();
        for (int i = 0; i < numberOfControlPointsPerRow; i++)
            N.insert(j, i) = m_properties->N(parameterization.at(j), i);
    }

    //endpoints that are to be interpolated
    Eigen::MatrixXd endPoints(2, 3);
    endPoints(0, 0) = dataPoints.first().x();
    endPoints(0, 1) = dataPoints.first().y();
    endPoints(0, 2) = dataPoints.first().z();
    endPoints(1, 0) = dataPoints.last().x();
    endPoints(1, 1) = dataPoints.last().y();
    endPoints(1, 2) = dataPoints.last().z();

    //Derivative of objective function (LS and energy term)
    Eigen::SparseMatrix<double> NtN_M = weightLS * N.transpose() * N + weightEnergy * curveEnergyMatrix;

    //Construct block matrix that needs to be inverted later
    const int blockMatDimension = numberOfControlPointsPerRow + 2;
    Eigen::SparseMatrix<double> blockMatrix(blockMatDimension, blockMatDimension);

    //TODO alternative way to do it: create from triples (Eigen documentation)
    for (int i = 0; i < numberOfControlPointsPerRow; i++)
        for (int j = 0; j < numberOfControlPointsPerRow; j++)
            if (NtN_M.coeff(i, j))
                blockMatrix.insert(i, j) = NtN_M.coeff(i, j);

    const double numericStabilityFactor = 1000;
    blockMatrix.insert(numberOfControlPointsPerRow    , 0                                ) = numericStabilityFactor;
    blockMatrix.insert(numberOfControlPointsPerRow + 1, numberOfControlPointsPerRow - 1) = numericStabilityFactor;
    blockMatrix.insert(0                                , numberOfControlPointsPerRow    ) = numericStabilityFactor;
    blockMatrix.insert(numberOfControlPointsPerRow - 1, numberOfControlPointsPerRow + 1) = numericStabilityFactor;

    //Compute result of blockMatrix * [b lambda] = [dataPoints boundary]
    Eigen::MatrixXd result(blockMatDimension, 3);
    for (int k = 0; k < 3; k++) {
        Eigen::VectorXd constVector(blockMatDimension);
        constVector << weightLS * N.transpose() * dataPointMatrix.col(k), numericStabilityFactor * endPoints.col(k);

        //TODO solution should have zero error (since we are solving derivative = 0), maybe use LR solver instead?
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver(blockMatrix);
        result.col(k) << solver.solve(constVector);
    }

    //Extract control points from result
    QVector<QVector3D> controlPoints;
    for (int j = 0; j < numberOfControlPointsPerRow; j++)
        controlPoints << QVector3D(result(j, 0), result(j, 1), result(j, 2));

    return new BSplineCurve(m_properties, controlPoints);
}

BSplineCurve *Splinefitter::fitSingleCurve_noConstraints(const double weightLS, const double weightEnergy, const QVector<QVector3D> &dataPoints, const QVector<double> &parameterization, Eigen::SparseMatrix<double> &curveEnergyMatrix)
{
    const int numberOfDataPoints = dataPoints.size();
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();

    //Set up least squares matrix
    Eigen::SparseMatrix<double> N(numberOfDataPoints, numberOfControlPointsPerRow);
    Eigen::MatrixXd dataPointMatrix(numberOfDataPoints, 3);

    for (int j = 0; j < numberOfDataPoints; j++) {
        dataPointMatrix(j, 0) = dataPoints[j].x();
        dataPointMatrix(j, 1) = dataPoints[j].y();
        dataPointMatrix(j, 2) = dataPoints[j].z();
        for (int i = 0; i < numberOfControlPointsPerRow; i++)
            N.insert(j, i) = m_properties->N(parameterization.at(j), i);
    }

    //Derivative of objective function (LS and energy term)
    Eigen::SparseMatrix<double> NtN_M = weightLS * N.transpose() * N + weightEnergy * curveEnergyMatrix;

    //Compute result of NtN_M * b = dataPoints
    Eigen::MatrixXd result(numberOfControlPointsPerRow, 3);
    for (int k = 0; k < 3; k++) {
        Eigen::VectorXd constVector = weightLS * N.transpose() * dataPointMatrix.col(k);

        //TODO solution should have zero error (since we are solving derivative = 0), maybe use LR solver instead?
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver(NtN_M);
        result.col(k) << solver.solve(constVector);
    }

    //Extract control points from result
    QVector<QVector3D> controlPoints;
    for (int j = 0; j < numberOfControlPointsPerRow; j++)
        controlPoints << QVector3D(result(j, 0), result(j, 1), result(j, 2));

    return new BSplineCurve(m_properties, controlPoints);
}

BSplineSurface *Splinefitter::fitSingleSurface_noConstraints(int cfId, const double weightLS, const double weightEnergy)
{
    CellFace *cf = m_cellMesh->getFace(cfId);
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();
    const int numberOfControlPointsPerSurface = numberOfControlPointsPerRow * numberOfControlPointsPerRow;

    QVector<Vertex *> *data = cf->getMeshVertices();
    const int numberOfLocalDataPoints = data->size();

    Parameterization *cfParam = m_cellMesh->getQuadCellParameterization(cfId);

    //Construct coefficient matrix and data point vector for approximation (least square system)
    Eigen::SparseMatrix<double> matrixLS(numberOfLocalDataPoints, numberOfControlPointsPerSurface);
//    const int nonZeroBasisFunctions = m_properties->getOrder() + 1;
//    matrixLS.reserve(Eigen::VectorXi::Constant(numberOfLocalDataPoints, nonZeroBasisFunctions * nonZeroBasisFunctions));
    Eigen::VectorXd dataPointsX(numberOfLocalDataPoints);
    Eigen::VectorXd dataPointsY(numberOfLocalDataPoints);
    Eigen::VectorXd dataPointsZ(numberOfLocalDataPoints);

    for (int k = 0; k < numberOfLocalDataPoints; k++) {
        Vertex *v = data->at(k);
        QVector2D vParam = cfParam->getParameter(v);

        dataPointsX[k] = v->getPosition().x();
        dataPointsY[k] = v->getPosition().y();
        dataPointsZ[k] = v->getPosition().z();

        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
            const double nu = m_properties->N(vParam.x(), i);
            for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                const double nv = m_properties->N(vParam.y(), j);
                const double product = nu*nv;
                if (product != 0)
                    matrixLS.insert(k, BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow)) = product;
            }
        }
    }

    double wLS = weightLS;
    double wEnergy = weightEnergy;
    if (!m_useFixedWeights) {
        wLS = (double) 1 / numberOfLocalDataPoints;
        wEnergy = 1 / cf->calculateAreaByCorners();
    }

    //calculate initial guess
    Eigen::SparseMatrix<double> matrixLStLS = wLS * matrixLS.transpose() * matrixLS + wEnergy * m_energyMatrix;

    Eigen::VectorXd rightSideX = wLS * matrixLS.transpose() * dataPointsX;
    Eigen::VectorXd rightSideY = wLS * matrixLS.transpose() * dataPointsY;
    Eigen::VectorXd rightSideZ = wLS * matrixLS.transpose() * dataPointsZ;

    //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.compute(matrixLStLS);
    Eigen::VectorXd ctrX = solver.solve(rightSideX);
    Eigen::VectorXd ctrY = solver.solve(rightSideY);
    Eigen::VectorXd ctrZ = solver.solve(rightSideZ);

    std::cout << matrixLStLS * ctrX - rightSideX << std::endl << "Norm: " << (matrixLStLS * ctrX - rightSideX).norm() << std::endl;
    std::cout << matrixLStLS * ctrY - rightSideY << std::endl << "Norm: " << (matrixLStLS * ctrY - rightSideY).norm() << std::endl;
    std::cout << matrixLStLS * ctrZ - rightSideZ << std::endl << "Norm: " << (matrixLStLS * ctrZ - rightSideZ).norm() << std::endl;



    //debug stuff
    const double lsErrorSqrtX = (matrixLS * ctrX - dataPointsX).norm();
    const double lsErrorSqrtY = (matrixLS * ctrY - dataPointsY).norm();
    const double lsErrorSqrtZ = (matrixLS * ctrZ - dataPointsZ).norm();

    const double lsErrorX = lsErrorSqrtX * lsErrorSqrtX;
    const double lsErrorY = lsErrorSqrtY * lsErrorSqrtY;
    const double lsErrorZ = lsErrorSqrtZ * lsErrorSqrtZ;

    const double lsErrorTotal = lsErrorX + lsErrorY + lsErrorZ;

    qDebug() << "LS Error:" << lsErrorTotal << "per datapoint:" << lsErrorTotal/numberOfLocalDataPoints << "datapoints:" << numberOfLocalDataPoints;
    qDebug() << "Individual Errors" << lsErrorX << lsErrorY << lsErrorZ;

    const double energyX = ctrX.transpose() * m_energyMatrix * ctrX;
    const double energyY = ctrY.transpose() * m_energyMatrix * ctrY;
    const double energyZ = ctrZ.transpose() * m_energyMatrix * ctrZ;

    const double energyTotal = energyX + energyY + energyZ;
    //TODO energy per area
    qDebug() << "Energy:" << energyTotal;
    qDebug() << "Individual energy terms:" << energyX << energyY << energyZ;
    qDebug() << "Weighted sum:" << 0.5 * wLS * lsErrorTotal + 0.5 * wEnergy * energyTotal;

    //Convert into 2D array of 3D points
    QVector<QVector<QVector3D> > controlPoints(numberOfControlPointsPerRow, QVector<QVector3D>(numberOfControlPointsPerRow));
    for (int i = 0; i < numberOfControlPointsPerRow; i++) {
        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
            const int ctrPointIndex = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
            const QVector3D ctrPoint(ctrX[ctrPointIndex], ctrY[ctrPointIndex], ctrZ[ctrPointIndex]);

            controlPoints[i][j] = ctrPoint;
        }
    }

    return new BSplineSurface(m_properties, m_properties, controlPoints);
}

BSplineSurface *Splinefitter::fitSingleSurface_boundaryInterpolation(int cfId, QVector<BSplineCurve *> boundaryCurves, const double weightLS, const double weightEnergy)
{
    CellFace *cf = m_cellMesh->getFace(cfId);
    const int numberOfControlPointsPerRow = m_properties->getNumberOfControlPoints();
    const int numberOfControlPointsPerSurface = numberOfControlPointsPerRow * numberOfControlPointsPerRow;

    //Minimize || Nb - d ||^2 + b^T M b (for each coordinate)
    // s.t. C b - b^(c) = 0     (boundary control points equal to control points of boundary curves)
    //Lagrange Multipliers:
    // /N^T N + M | C^T\ /    b   \ = / N^T d \  (C contains exactly one "1" in each row)
    // \    C     |  0 / \ lambda /   \ b^(c) /

    //Construct C (leave out the last CP per row to avoid duplicates)
    const double numericStabilityFactor = 1000;
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(4 * (numberOfControlPointsPerRow - 1), numberOfControlPointsPerSurface);

    for (int i = 0; i < numberOfControlPointsPerRow - 1; i++) {
        C(i                                  ,
          BSplineSurface::matToVecIndexLocal(i,
                                             0,
                                             numberOfControlPointsPerRow)
          ) = numericStabilityFactor;
        C(i +   (numberOfControlPointsPerRow - 1),
          BSplineSurface::matToVecIndexLocal(numberOfControlPointsPerRow - 1,
                                             i,
                                             numberOfControlPointsPerRow)
          ) = numericStabilityFactor;
        C(i + 2*(numberOfControlPointsPerRow - 1),
          BSplineSurface::matToVecIndexLocal(numberOfControlPointsPerRow - 1 - i,
                                             numberOfControlPointsPerRow - 1,
                                             numberOfControlPointsPerRow)
          ) = numericStabilityFactor;
        C(i + 3*(numberOfControlPointsPerRow - 1),
          BSplineSurface::matToVecIndexLocal(0,
                                             numberOfControlPointsPerRow - 1 - i,
                                             numberOfControlPointsPerRow)
          ) = numericStabilityFactor;
    }

    Eigen::MatrixXd borderCP(4 * (numberOfControlPointsPerRow - 1), 3);
    for (int iCe = 0; iCe < 4; iCe++) {
        const QVector<QVector3D> *contrPoints = boundaryCurves[iCe]->getControlPointsP();

        if (!cf->edgeIsInverted(iCe)) {
            for (int j = 0; j < numberOfControlPointsPerRow-1; j++) {
                borderCP(iCe * (numberOfControlPointsPerRow - 1) + j, 0) = contrPoints->at(j).x();
                borderCP(iCe * (numberOfControlPointsPerRow - 1) + j, 1) = contrPoints->at(j).y();
                borderCP(iCe * (numberOfControlPointsPerRow - 1) + j, 2) = contrPoints->at(j).z();
            }
        } else {
            for (int j = 0; j < numberOfControlPointsPerRow-1; j++) {
                borderCP(iCe * (numberOfControlPointsPerRow - 1) + j, 0) = contrPoints->at(numberOfControlPointsPerRow-1-j).x();
                borderCP(iCe * (numberOfControlPointsPerRow - 1) + j, 1) = contrPoints->at(numberOfControlPointsPerRow-1-j).y();
                borderCP(iCe * (numberOfControlPointsPerRow - 1) + j, 2) = contrPoints->at(numberOfControlPointsPerRow-1-j).z();
            }
        }
    }
    borderCP *= numericStabilityFactor;

    Eigen::SparseMatrix<double, Eigen::RowMajor> N = this->constructCoefficientMatrixForOneCell(cfId, m_properties);
    Eigen::MatrixXd dataPoints = this->constructDataPointVector(m_cellMesh->getFace(cfId)->getMeshVertices());

    const int blockMatDimension = numberOfControlPointsPerSurface + 4 * (numberOfControlPointsPerRow - 1);

//    Eigen::SparseMatrix<double, Eigen::ColMajor> blockMatrix = weightLS * N.transpose() * N + weightEnergy * m_energyMatrix;
//    blockMatrix.resize(blockMatDimension, blockMatDimension);
//    blockMatrix.reserve(Eigen::VectorXi::Constant(blockMatDimension, numberOfControlPointsPerSurface + 2));

    Eigen::SparseMatrix<double> NtN_M = weightLS * N.transpose() * N + weightEnergy * m_energyMatrix;
    //Construct block matrix
    Eigen::SparseMatrix<double, Eigen::ColMajor> blockMatrix(blockMatDimension, blockMatDimension);
    blockMatrix.reserve(Eigen::VectorXi::Constant(blockMatDimension, numberOfControlPointsPerSurface + 2));


    //TODO alternative way to do it: create from triples (Eigen documentation)
    for (int i = 0; i < numberOfControlPointsPerSurface; i++)
        for (int j = 0; j < numberOfControlPointsPerSurface; j++)
            if (NtN_M.coeff(i, j))
                blockMatrix.insert(i, j) = NtN_M.coeff(i, j);

    for (int i = 0; i < 4 * (numberOfControlPointsPerRow - 1); i++)
        for (int j = 0; j < numberOfControlPointsPerSurface; j++) {
            const double coef = C(i, j);
            if (coef) {
                blockMatrix.insert(numberOfControlPointsPerSurface + i, j) = coef;
                blockMatrix.insert(j, numberOfControlPointsPerSurface + i) = coef;
            }
        }

    Eigen::MatrixXd result(blockMatDimension, 3);
    blockMatrix.makeCompressed();
    for (int k = 0; k < 3; k++) {
        Eigen::VectorXd constVector(blockMatDimension);
        constVector << weightLS * N.transpose() * dataPoints.col(k), borderCP.col(k);

        //TODO solution should have zero error (since we are solving derivative = 0), maybe use LR solver instead?
        //Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver(blockMatrix);
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solver(blockMatrix);
        solver.analyzePattern(blockMatrix);
        solver.factorize(blockMatrix);
        result.col(k) << solver.solve(constVector);
    }

    //Convert into 2D array of 3D points
    QVector<QVector<QVector3D> > controlPoints(numberOfControlPointsPerRow, QVector<QVector3D>(numberOfControlPointsPerRow));
    for (int i = 0; i < numberOfControlPointsPerRow; i++) {
        for (int j = 0; j < numberOfControlPointsPerRow; j++) {
            const int ctrPointIndex = BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow);
            const QVector3D ctrPoint(result(ctrPointIndex, 0), result(ctrPointIndex, 1), result(ctrPointIndex, 2));

            controlPoints[i][j] = ctrPoint;
        }
    }

    return new BSplineSurface(m_properties, m_properties, controlPoints);
}

Eigen::SparseMatrix<double, Eigen::RowMajor> Splinefitter::constructCoefficientMatrixForOneCell(int cfId, BSplineProperties *properties)
{
    CellFace *cf = m_cellMesh->getFace(cfId);
    Parameterization *cfParam = m_cellMesh->getQuadCellParameterization(cfId);

    const int numberOfControlPointsPerRow = properties->getNumberOfControlPoints();
    const int numberOfControlPointsPerSurface = numberOfControlPointsPerRow * numberOfControlPointsPerRow;
    const int numberOfLocalDataPoints = cf->getMeshVertices()->size();
    Eigen::SparseMatrix<double, Eigen::RowMajor> matrixLS(numberOfLocalDataPoints, numberOfControlPointsPerSurface);
    const int nonZeroBasisFunctions = properties->getOrder() + 1;
    matrixLS.reserve(Eigen::VectorXi::Constant(numberOfLocalDataPoints, nonZeroBasisFunctions * nonZeroBasisFunctions));
    for (int k = 0; k < numberOfLocalDataPoints; k++) {
        Vertex *v = cf->getMeshVertices()->at(k);
        QVector2D vParam = cfParam->getParameter(v);

        for (int i = 0; i < numberOfControlPointsPerRow; i++) {
            const double nu = properties->N(vParam.x(), i);
            if (nu != 0) {
                for (int j = 0; j < numberOfControlPointsPerRow; j++) {
                    const double nv = properties->N(vParam.y(), j);
                    const double product = nu * nv;
                    if (product != 0)
                        matrixLS.insert(k, BSplineSurface::matToVecIndexLocal(i, j, numberOfControlPointsPerRow)) = product;
                }
            }
        }
    }

    return matrixLS;
}

Eigen::MatrixXd Splinefitter::constructDataPointVector(QVector<Vertex *> *dataPoints)
{
    const int numberOfLocalDataPoints = dataPoints->size();
    Eigen::MatrixXd dataPointMatrix(numberOfLocalDataPoints, 3);

    for (int k = 0; k < numberOfLocalDataPoints; k++) {
        Vertex *v = dataPoints->at(k);

        dataPointMatrix(k,0) = v->getPosition().x();
        dataPointMatrix(k,1) = v->getPosition().y();
        dataPointMatrix(k,2) = v->getPosition().z();
    }

    return dataPointMatrix;
}

void Splinefitter::setInitialGuess(BSplinePatchNetwork *initialGuess)
{
    m_initialGuess = initialGuess;
}

BSplineProperties *Splinefitter::getBSplineProperties() const
{
    return m_properties;
}

void Splinefitter::setBSplineProperties(BSplineProperties *properties)
{
    m_properties = properties;
}

Eigen::SparseMatrix<double> Splinefitter::getEnergyMatrix() const
{
    return m_energyMatrix;
}

void Splinefitter::setEnergyMatrix(const Eigen::SparseMatrix<double> &energyMatrix)
{
    m_energyMatrix = energyMatrix;
}

double Splinefitter::getWeightLSTerm() const
{
    return m_weightLSTerm;
}

void Splinefitter::setWeightLSTerm(double weightLSTerm, bool convexCombination)
{
    m_weightLSTerm = weightLSTerm;
    if (convexCombination)
        m_weightEnergyTerm = 1 - weightLSTerm;
}
double Splinefitter::getWeightEnergyTerm() const
{
    return m_weightEnergyTerm;
}

void Splinefitter::setWeightEnergyTerm(double weightEnergyTerm, bool convexCombination)
{
    m_weightEnergyTerm = weightEnergyTerm;
    if (convexCombination)
        m_weightLSTerm = 1 - weightEnergyTerm;
}

bool Splinefitter::getUseFixedWeights() const
{
    return m_useFixedWeights;
}

void Splinefitter::setUseFixedWeights(bool useFixedWeights)
{
    m_useFixedWeights = useFixedWeights;
}

int Splinefitter::getNumberOfG1VertexPoints() const
{
    return m_numberOfG1VertexPoints;
}

void Splinefitter::setNumberOfG1VertexPoints(int numberOfG1VertexPoints)
{
    m_numberOfG1VertexPoints = numberOfG1VertexPoints;
}
double Splinefitter::getSizeOfG1VertexDomain() const
{
    return m_sizeOfG1VertexDomain;
}

void Splinefitter::setSizeOfG1VertexDomain(double sizeOfG1VertexDomain)
{
    m_sizeOfG1VertexDomain = sizeOfG1VertexDomain;
}
int Splinefitter::getNumberOfG1EdgePoints() const
{
    return m_numberOfG1EdgePoints;
}

void Splinefitter::setNumberOfG1EdgePoints(int numberOfG1EdgePoints)
{
    m_numberOfG1EdgePoints = numberOfG1EdgePoints;
}

bool Splinefitter::getG1VertexConstraintsAsOpt() const
{
    return m_g1VertexConstraintsAsOpt;
}

void Splinefitter::setG1VertexConstraintsAsOpt(bool g1VertexConstraintsAsOpt)
{
    m_g1VertexConstraintsAsOpt = g1VertexConstraintsAsOpt;
}

double Splinefitter::getG1VertexConstraintsOptFactor() const
{
    return m_g1VertexConstraintsOptFactor;
}

void Splinefitter::setG1VertexConstraintsOptFactor(double g1VertexConstraintsOptFactor)
{
    m_g1VertexConstraintsOptFactor = g1VertexConstraintsOptFactor;
}

bool Splinefitter::getG1EdgeConstraintsAsOpt() const
{
    return m_g1EdgeConstraintsAsOpt;
}

void Splinefitter::setG1EdgeConstraintsAsOpt(bool g1EdgeConstraintsAsOpt)
{
    m_g1EdgeConstraintsAsOpt = g1EdgeConstraintsAsOpt;
}

double Splinefitter::getG1EdgeConstraintsOptFactor() const
{
    return m_g1EdgeConstraintsOptFactor;
}

void Splinefitter::setG1EdgeConstraintsOptFactor(double g1EdgeConstraintsOptFactor)
{
    m_g1EdgeConstraintsOptFactor = g1EdgeConstraintsOptFactor;
}

int Splinefitter::getNloMaxIterations() const
{
    return m_nloMaxIterations;
}

void Splinefitter::setNloMaxIterations(int nloIterations)
{
    m_nloMaxIterations = nloIterations;
}
double Splinefitter::getNloFPrec() const
{
    return m_nloFPrec;
}

void Splinefitter::setNloFPrec(double nloFPrec)
{
    m_nloFPrec = nloFPrec;
}
double Splinefitter::getNloXPrec() const
{
    return m_nloXPrec;
}

void Splinefitter::setNloXPrec(double nloXPrec)
{
    m_nloXPrec = nloXPrec;
}
double Splinefitter::getNloG0constPrec() const
{
    return m_nloG0constPrec;
}

void Splinefitter::setNloG0constPrec(double nloG0constPrec)
{
    m_nloG0constPrec = nloG0constPrec;
}
double Splinefitter::getNloEqualConstPrec() const
{
    return m_nloEqualConstPrec;
}

void Splinefitter::setNloEqualConstPrec(double nloEqualConstPrec)
{
    m_nloEqualConstPrec = nloEqualConstPrec;
}
double Splinefitter::getNloG1constPrec() const
{
    return m_nloG1constPrec;
}

void Splinefitter::setNloG1constPrec(double nloG1constPrec)
{
    m_nloG1constPrec = nloG1constPrec;
}

bool Splinefitter::getDetailedOutput() const
{
    return m_detailedOutput;
}

void Splinefitter::setDetailedOutput(bool detailedOutput)
{
    m_detailedOutput = detailedOutput;
}













