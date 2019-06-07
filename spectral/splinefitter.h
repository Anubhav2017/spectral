#ifndef SPLINEFITTER_H
#define SPLINEFITTER_H

#include "BSpline/bsplinepatchnetwork.h"
#include "CellMesh/cellmesh.h"

// Eigen
#include <Sparse>
#include <Dense>

//Qt
#include <QVector>

class Splinefitter
{
public:
    Splinefitter(CellMesh *cm, BSplineProperties *properties);

    BSplinePatchNetwork *doApproximationLocalG0ConstraintG1Sample();
    BSplinePatchNetwork *doApproximationGlobalC0G1Constraints(bool onlyC0 = false);
    BSplinePatchNetwork *doApproximationLocalC0WithBoundaryInterpolation(bool minimizeEnergyOnBoundaries, double weightLSBoundary, double weightEBoundary);

    BSplinePatchNetwork *doIterativeCurveAndSurfaceFittingC0BoundaryInterpolation(const QVector<double> &energyWeightsCurve, const QVector<double> &energyWeightsSurf, const bool useNonlinearSolverForOpt);

    BSplineCurve *fitSingleCurve_endPointInterpol(const double weightLS, const double weightEnergy, const QVector<QVector3D> &dataPoints, const QVector<double> &parameterization, Eigen::SparseMatrix<double> &curveEnergyMatrix);
    BSplineCurve *fitSingleCurve_noConstraints(const double weightLS, const double weightEnergy, const QVector<QVector3D> &dataPoints, const QVector<double> &parameterization, Eigen::SparseMatrix<double> &curveEnergyMatrix);
    BSplineSurface *fitSingleSurface_noConstraints(int cfId, const double weightLS, const double weightEnergy);
    BSplineSurface *fitSingleSurface_boundaryInterpolation(int cfId, QVector<BSplineCurve *> boundaryCurves, const double weightLS, const double weightEnergy);

    Eigen::SparseMatrix<double, Eigen::RowMajor> constructCoefficientMatrixForOneCell(int cfId, BSplineProperties *properties);
    Eigen::MatrixXd constructDataPointVector(QVector<Vertex *> *dataPoints);


    //getter and setter methods
    void setInitialGuess(BSplinePatchNetwork *initialGuess);

    BSplineProperties *getBSplineProperties() const;
    void setBSplineProperties(BSplineProperties *properties);

    Eigen::SparseMatrix<double> getEnergyMatrix() const;
    void setEnergyMatrix(const Eigen::SparseMatrix<double> &energyMatrix);

    double getWeightLSTerm() const;
    void setWeightLSTerm(double weightLSTerm, bool convexCombination = true);

    double getWeightEnergyTerm() const;
    void setWeightEnergyTerm(double weightEnergyTerm, bool convexCombination = true);

    bool getUseFixedWeights() const;
    void setUseFixedWeights(bool useFixedWeights);

    int getNumberOfG1VertexPoints() const;
    void setNumberOfG1VertexPoints(int numberOfG1VertexPoints);

    double getSizeOfG1VertexDomain() const;
    void setSizeOfG1VertexDomain(double sizeOfG1VertexDomain);

    int getNumberOfG1EdgePoints() const;
    void setNumberOfG1EdgePoints(int numberOfG1EdgePoints);

    bool getG1VertexConstraintsAsOpt() const;
    void setG1VertexConstraintsAsOpt(bool g1VertexConstraintsAsOpt);

    double getG1VertexConstraintsOptFactor() const;
    void setG1VertexConstraintsOptFactor(double g1VertexConstraintsOptFactor);

    bool getG1EdgeConstraintsAsOpt() const;
    void setG1EdgeConstraintsAsOpt(bool g1EdgeConstraintsAsOpt);

    double getG1EdgeConstraintsOptFactor() const;
    void setG1EdgeConstraintsOptFactor(double g1EdgeConstraintsOptFactor);

    int getNloMaxIterations() const;
    void setNloMaxIterations(int nloIterations);

    double getNloFPrec() const;
    void setNloFPrec(double nloFPrec);

    double getNloXPrec() const;
    void setNloXPrec(double nloXPrec);

    double getNloG0constPrec() const;
    void setNloG0constPrec(double nloG0constPrec);

    double getNloEqualConstPrec() const;
    void setNloEqualConstPrec(double nloEqualConstPrec);

    double getNloG1constPrec() const;
    void setNloG1constPrec(double nloG1constPrec);

    bool getDetailedOutput() const;
    void setDetailedOutput(bool detailedOutput);

private:
    //Data
    CellMesh *const m_cellMesh;
    BSplinePatchNetwork *m_initialGuess;

    //General B-Spline parameters
    BSplineProperties *m_properties;

    //Energy matrix and weights
    Eigen::SparseMatrix<double> m_energyMatrix;
    double m_weightLSTerm;
    double m_weightEnergyTerm;
    bool m_useFixedWeights;

    //Sample G1-Algorithm parameters
    int m_numberOfG1VertexPoints;
    double m_sizeOfG1VertexDomain;
    int m_numberOfG1EdgePoints;
    bool m_g1VertexConstraintsAsOpt;
    double m_g1VertexConstraintsOptFactor;
    bool m_g1EdgeConstraintsAsOpt;
    double m_g1EdgeConstraintsOptFactor;


    //NLO parameters
    int m_nloMaxIterations;
    double m_nloFPrec;
    double m_nloXPrec;
    double m_nloG0constPrec;
    double m_nloEqualConstPrec;
    double m_nloG1constPrec;

    //Debug paramewters
    bool m_detailedOutput;
};

#endif // SPLINEFITTER_H
