#include "parameteroptimizer.h"

#include "CellMesh/parameterization.h"

//Eigen
#include <Dense>

//Omp
#include <omp.h>

//NLOpt
#include <nlopt.h>

//Qt
#include <QDebug>

struct nlo_distPointSurfaceParam {
    QVector3D point;
    BSplineSurface *surface;
};

double nlo_distPointSurface(unsigned int n, const double *x, double *grad, void *params) {
    nlo_distPointSurfaceParam *parameters = (nlo_distPointSurfaceParam *) params;
    BSplineSurface *surface = parameters->surface;

    const double u = x[0];
    const double v = x[1];

    //TODO evaluate N(u,v), Nu(u,v), Nv(u,v) as a vector once and multiply the control point vector
    const QVector3D S = surface->evaluate(u, v);
    const QVector3D Tu = surface->evaluateTu(u, v);
    const QVector3D Tv = surface->evaluateTv(u, v);

    const QVector3D errorVector = S - parameters->point;

    if (grad) {
        grad[0] = QVector3D::dotProduct(Tu, errorVector);
        grad[1] = QVector3D::dotProduct(Tv, errorVector);

        for (unsigned int i = 2; i < n; i++)
            grad[i] = 0;
    }

    return 0.5 * errorVector.lengthSquared();
}

struct nlo_distPointCurveParam {
    QVector3D point;
    BSplineCurve *curve;
};

double nlo_distPointCurve(unsigned int n, const double *x, double *grad, void *params) {
    nlo_distPointCurveParam *parameters = (nlo_distPointCurveParam *) params;
    BSplineCurve *curve = parameters->curve;

    const double t = *x;

    //TODO evaluate N(u,v), Nu(u,v), Nv(u,v) as a vector once and multiply the control point vector
    const QVector3D C = curve->evaluate(t);
    const QVector3D tangent = curve->evaluateTangent(t);

    const QVector3D errorVector = C - parameters->point;

    if (grad) {
        grad[0] = QVector3D::dotProduct(tangent, errorVector);

        for (unsigned int i = 1; i < n; i++)
            grad[i] = 0;
    }

    return 0.5 * errorVector.lengthSquared();
}

double ParameterOptimizer::optimizeParameterizationFullSurfaceNetwork(BSplinePatchNetwork *bSplineSurfaceNetwork, const bool alsoOptimizeBoundaries, const bool useNonlinearSolver)
{
    CellMesh *cm = bSplineSurfaceNetwork->getCellMesh();

    const int numberOfSurfaces = bSplineSurfaceNetwork->getNumberOfBSplineSurfaces();
    const int numberOfDataPoints = cm->getOriginalMesh()->getNumberOfVertices();

    QVector<double> dataErrorsBefore(numberOfDataPoints, 0);
    QVector<double> dataErrorsAfter(numberOfDataPoints, 0);

    if (alsoOptimizeBoundaries)
        ParameterOptimizer::optimizeParameterizationFullCurveNetwork(bSplineSurfaceNetwork, useNonlinearSolver);

    #pragma omp parallel for
    for (int iCell = 0; iCell < numberOfSurfaces; iCell++) {
        CellFace *cf = cm->getFace(iCell);
        Parameterization *cfParam = cm->getQuadCellParameterization(iCell);
        BSplineSurface *bs = bSplineSurfaceNetwork->getBSplineSurface(iCell);

        const int numberOfVerticesInCellFace = cf->getMeshVertices()->size();
        double cfErrorBefore = 0;
        double cfErrorAfter = 0;
        for (int k = 0; k < numberOfVerticesInCellFace; k++) {
            Vertex *v = cf->getMeshVertices()->at(k);
            QVector2D param = cfParam->getParameter(v);

            const QVector3D vertexPos = v->getPosition();
            const double errorSquaredBefore = (bs->evaluate(param.x(), param.y()) - vertexPos).lengthSquared();

            if (useNonlinearSolver)
                param = ParameterOptimizer::optimizeParameter_surfacePoint_nlo_bounded(bs, vertexPos, param);
            else
                param = ParameterOptimizer::optimizeParameter_surfacePoint_projection_bounded(bs, vertexPos, param);

            const double errorSquaredAfter = (bs->evaluate(param.x(), param.y()) - vertexPos).lengthSquared();

            cfErrorBefore += errorSquaredBefore;
            dataErrorsBefore[v->getId()] = errorSquaredBefore;
            if (errorSquaredAfter < errorSquaredBefore) {
                cfParam->updateParameter(v, param.x(), param.y());
                cfErrorAfter += errorSquaredAfter;
                dataErrorsAfter[v->getId()] = errorSquaredAfter;
            } else {
                cfErrorAfter += errorSquaredBefore;
                dataErrorsAfter[v->getId()] = errorSquaredBefore;
            }
        }

        if (useNonlinearSolver)
            cfParam->addInfo(QString("Optimized by finding nearest surface points (NLO). Sum of squared errors before: %1, after: %2\n").arg(cfErrorBefore).arg(cfErrorAfter));
        else
            cfParam->addInfo(QString("Optimized by projection onto tangent plane of given B-Splines. Sum of squared errors before: %1, after: %2\n").arg(cfErrorBefore).arg(cfErrorAfter));
    }

    double lsErrorBefore = 0;
    double lsErrorAfter = 0;
    for (int i = 0; i < numberOfDataPoints; i++) {
        lsErrorBefore += dataErrorsBefore[i];
        lsErrorAfter += dataErrorsAfter[i];
    }

    return lsErrorAfter / lsErrorBefore;
}

void ParameterOptimizer::optimizeParameterizationFullCurveNetwork(BSplinePatchNetwork *bSplineSurfaceNetwork, const bool useNonlinearSolver)
{
    const int numberOfCurves = bSplineSurfaceNetwork->getNumberOfBSplineCurves();

    CellMesh *cm = bSplineSurfaceNetwork->getCellMesh();

    #pragma omp parallel for
    for (int iCurv = 0; iCurv < numberOfCurves; iCurv++) {
        CellEdge *ce = cm->getEdge(iCurv);
        BSplineCurve *curve = bSplineSurfaceNetwork->getBSplineCurve(iCurv);

        const int incidentSurfaceId = ce->getIncidentFaceId(0);
        Parameterization *cfParam = cm->getQuadCellParameterization(incidentSurfaceId);

        CellFace *incidentCf = cm->getFace(incidentSurfaceId);   //should be the same as if extracted from cfParam[1]
        QVector<double> ceOldParam = cfParam->getParameterizationOfBoundaryCurve(incidentCf->getInternalEdgeId(ce));
        const int numPLV = ce->getPolylineVertices()->size();
        QVector<QVector3D> dataPoints(numPLV);
        for (int i = 0; i < numPLV; i++)
            dataPoints[i] = ce->getPolylineVertices()->at(i)->getPosition();

        QVector<double> ceOptParam = ParameterOptimizer::optimizeParameterizationSingleCurve(curve, dataPoints, ceOldParam, true, useNonlinearSolver);

        cm->updateBoundaryParameterization(iCurv, ceOptParam, true);
    }
}

QVector<double> ParameterOptimizer::optimizeParameterizationSingleCurve(BSplineCurve * const curve, const QVector<QVector3D> &dataPoints, const QVector<double> &parameters, const bool keepEndpointsFixed, const bool useNonlinearSolver)
{
    const int numDataPoints = dataPoints.size();
    if (numDataPoints != parameters.size())
        return QVector<double>();

    int skipEnds = 0;
    if (keepEndpointsFixed)
        skipEnds = 1;

    QVector<double> result(numDataPoints, 0);

    for (int i = skipEnds; i < numDataPoints - skipEnds; i++) {
        const double tInitial = parameters[i];
        const QVector3D point = dataPoints[i];

        const double errorSquaredBefore = (curve->evaluate(tInitial) - point).lengthSquared();

        double tNew = -1;
        if (useNonlinearSolver)
            tNew = ParameterOptimizer::optimizeParameter_curvePoint_nlo_bounded(curve, point, tInitial);
        else
            tNew = ParameterOptimizer::optimizeParameter_curvePoint_projection_bounded(curve, point, tInitial);

        const double errorSquaredAfter = (curve->evaluate(tNew) - point).lengthSquared();

        if (errorSquaredAfter < errorSquaredBefore)
            result[i] = tNew;
        else
            result[i] = tInitial;
    }

    return result;
}

QVector2D ParameterOptimizer::optimizeParameter_surfacePoint_nlo_bounded(BSplineSurface *surface, const QVector3D &point, const QVector2D &initialParameter, const int nloptMaxIterations, const double nloptPrecision)
{
    nlo_distPointSurfaceParam parameters;
    parameters.surface = surface;
    parameters.point = point;

    //NLOPT_LD_MMA -> return value 5, 1000 iterations yield good results and take 20 ms, no errors
    //NLOPT_LD_SLSQP -> also good results but for some reasons often errors, same time as MMA
    //NLOPT_LD_LBFGS -> very fast, good results, almost always error -1

    nlopt_opt optimizer = nlopt_create(NLOPT_LD_LBFGS, 2);
    nlopt_set_min_objective(optimizer, nlo_distPointSurface, &parameters);
    nlopt_set_ftol_rel(optimizer, nloptPrecision);
    nlopt_set_xtol_rel(optimizer, nloptPrecision);
    nlopt_set_maxeval(optimizer, nloptMaxIterations);

    nlopt_set_lower_bounds1(optimizer, 0);
    nlopt_set_upper_bounds1(optimizer, 1);

    double optVariables[2] = {initialParameter.x(), initialParameter.y()};

    double optValue;
    nlopt_optimize(optimizer, optVariables, &optValue);
//    const int status = nlopt_optimize(optimizer, optVariables, &optValue);
//    if (status < 0)
//        qDebug() << "ERROR while parameter optimization! Status:" << status << "Value:" << optValue << "Parameter: (" << optVariables[0] << "|" << optVariables[1] << ")";

    //qDebug() << "Status:" << status << "Value:" << optValue << "Parameter: (" << optVariables[0] << "|" << optVariables[1] << ")" << "time:" << timer.elapsed();

    nlopt_destroy(optimizer);

    //const double initialError = (surface->evaluate(initialGuess.x(), initialGuess.y()) - point).length();
    //const double newError = (surface->evaluate(optVariables[0], optVariables[1]) - point).length();
    //qDebug() << "Error before:" << initialError << "after:" << newError;

    return QVector2D(optVariables[0], optVariables[1]);
}

QVector2D ParameterOptimizer::optimizeParameter_surfacePoint_projection_unbounded(BSplineSurface *surface, const QVector3D &point, const QVector2D &initialParameter)
{
    const QVector3D bSplinePosQVector = surface->evaluate(initialParameter.x(), initialParameter.y());
    const QVector3D Tu = surface->evaluateTu(initialParameter.x(), initialParameter.y());
    const QVector3D Tv = surface->evaluateTv(initialParameter.x(), initialParameter.y());
    const QVector3D normal = QVector3D::crossProduct(Tu, Tv);

    const Eigen::Vector3d pointPos(point.x(), point.y(), point.z());
    const Eigen::Vector3d bSplinePos(bSplinePosQVector.x(), bSplinePosQVector.y(), bSplinePosQVector.z());
    Eigen::Matrix3d coefficients;
    coefficients << Tu.x(), Tv.x(), normal.x(),
                    Tu.y(), Tv.y(), normal.y(),
                    Tu.z(), Tv.z(), normal.z();
    const Eigen::Vector3d factors = coefficients.inverse() * (pointPos - bSplinePos);

    return initialParameter + QVector2D(factors[0], factors[1]);
}

QVector2D ParameterOptimizer::optimizeParameter_surfacePoint_projection_bounded(BSplineSurface *surface, const QVector3D &point, const QVector2D &initialParameter)
{
    const QVector3D bSplinePosQVector = surface->evaluate(initialParameter.x(), initialParameter.y());
    const QVector3D Tu = surface->evaluateTu(initialParameter.x(), initialParameter.y());
    const QVector3D Tv = surface->evaluateTv(initialParameter.x(), initialParameter.y());
    const QVector3D normal = QVector3D::crossProduct(Tu, Tv);

    const Eigen::Vector3d pointPos(point.x(), point.y(), point.z());
    const Eigen::Vector3d bSplinePos(bSplinePosQVector.x(), bSplinePosQVector.y(), bSplinePosQVector.z());
    Eigen::Matrix3d coefficients;
    coefficients << Tu.x(), Tv.x(), normal.x(),
                    Tu.y(), Tv.y(), normal.y(),
                    Tu.z(), Tv.z(), normal.z();
    const Eigen::Vector3d factors = coefficients.inverse() * (pointPos - bSplinePos);
    const QVector2D paramDirection = QVector2D(factors[0], factors[1]);

    QVector2D newParam = initialParameter + paramDirection;
    double stepsize = 1.0;
    while (newParam.x() < 0 || newParam.x() > 1 || newParam.y() < 0 || newParam.y() > 1) {
        stepsize /= 2;
        newParam = initialParameter + (stepsize * paramDirection);
    }

    return newParam;
}

double ParameterOptimizer::optimizeParameter_curvePoint_nlo_bounded(BSplineCurve *curve, const QVector3D &point, const double initialParameter, const int nloptMaxIterations, const double nloptPrecision)
{
    nlo_distPointCurveParam parameters;
    parameters.curve = curve;
    parameters.point = point;

    nlopt_opt optimizer = nlopt_create(NLOPT_LD_LBFGS, 1);
    nlopt_set_min_objective(optimizer, nlo_distPointCurve, &parameters);
    nlopt_set_ftol_rel(optimizer, nloptPrecision);
    nlopt_set_xtol_rel(optimizer, nloptPrecision);
    nlopt_set_maxeval(optimizer, nloptMaxIterations);

    nlopt_set_lower_bounds1(optimizer, 0);
    nlopt_set_upper_bounds1(optimizer, 1);

    double optVariable = initialParameter;

    double optValue;
    nlopt_optimize(optimizer, &optVariable, &optValue);
//    const int status = nlopt_optimize(optimizer, &optVariable, &optValue);
//    if (status < 0)
//        qDebug() << "ERROR while parameter optimization! Status:" << status << "Value:" << optValue << "Parameter:" << optVariable;

//    qDebug() << "Status:" << status << "Value:" << optValue << "Parameter:" << optVariable;

    nlopt_destroy(optimizer);

    return optVariable;
}

double ParameterOptimizer::optimizeParameter_curvePoint_projection_unbounded(BSplineCurve *curve, const QVector3D &point, const double initialParameter)
{
    const QVector3D curvePos = curve->evaluate(initialParameter);
    const QVector3D tangent = curve->evaluateTangent(initialParameter);

    const QVector3D errorVec = point - curvePos;

    //QVector3D projectedPoint = QVector3D::dotProduct(errorVec, tangent)/tangent.lengthSquared() * tangent;
    const double factor = QVector3D::dotProduct(errorVec, tangent)/tangent.lengthSquared();    //we only need the factor of the projection
    return initialParameter + factor;
}

double ParameterOptimizer::optimizeParameter_curvePoint_projection_bounded(BSplineCurve *curve, const QVector3D &point, const double initialParameter)
{
    const QVector3D curvePos = curve->evaluate(initialParameter);
    const QVector3D tangent = curve->evaluateTangent(initialParameter);

    const QVector3D errorVec = point - curvePos;

    //QVector3D projectedPoint = QVector3D::dotProduct(errorVec, tangent)/tangent.lengthSquared() * tangent;
    const double factor = QVector3D::dotProduct(errorVec, tangent)/tangent.lengthSquared();    //we only need the factor of the projection
    const double newParameter = initialParameter + factor;
    if (newParameter < 0)
        return 0;
    else if (newParameter > 1)
        return 1;
    else
        return newParameter;
}
